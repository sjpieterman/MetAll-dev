#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <zlib.h>
#include "common.h"
#include "index.h"
#include "classifier.h"
#include "fastq_reader.h"
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <memory>
#include <atomic>
#include <map>

namespace fs = std::filesystem;

// Kraken rank code mapping
std::string get_rank_code(const char* rank_name, const char* scientific_name) {
    std::string r = rank_name;
    std::string n = scientific_name;
    
    if (n == "root") return "R";
    if (n == "unclassified") return "U";
    
    if (r == "superkingdom" || r == "domain") return "D";
    if (r == "kingdom") return "K";
    if (r == "phylum") return "P";
    if (r == "class") return "C";
    if (r == "order") return "O";
    if (r == "family") return "F";
    if (r == "genus") return "G";
    if (r == "species") return "S";
    if (r == "subspecies") return "S1";
    if (r == "strain") return "S2";
    
    // Prefix matches for sub-ranks
    if (r.find("sub") == 0) {
        if (r.find("kingdom") != std::string::npos) return "K1";
        if (r.find("phylum") != std::string::npos) return "P1";
        if (r.find("class") != std::string::npos) return "C1";
        if (r.find("order") != std::string::npos) return "O1";
        if (r.find("family") != std::string::npos) return "F1";
        if (r.find("genus") != std::string::npos) return "G1";
        if (r.find("species") != std::string::npos) return "S1";
    }
    
    // Root level organisms often get R1, R2 etc.
    if (r == "cellular root" || r == "acellular root" || n == "cellular organisms") return "R1";

    return "-";
}

void print_report_node(uint64_t idx, int depth, uint64_t total_reads, 
                       const std::vector<uint64_t>& direct_counts, 
                       const std::vector<uint64_t>& clade_counts,
                       const KmerIndex& index, std::ostream& out) {
    const auto& nodes = index.getNodes();
    if (clade_counts[idx] == 0) return;

    const char* name_data = index.getNameData().data();
    const char* rank_data = index.getRankData().data();
    const char* scientific_name = &name_data[nodes[idx].name_offset];
    const char* rank_name = &rank_data[nodes[idx].rank_offset];
    
    double pct = (total_reads > 0) ? (clade_counts[idx] * 100.0 / total_reads) : 0.0;
    
    out << std::fixed << std::setprecision(2) << std::setw(6) << pct << "\t"
        << clade_counts[idx] << "\t"
        << direct_counts[idx] << "\t"
        << get_rank_code(rank_name, scientific_name) << "\t"
        << nodes[idx].external_id << "\t";
    
    for (int i = 0; i < depth; ++i) out << "  ";
    out << scientific_name << "\n";

    for (uint64_t i = 0; i < nodes[idx].child_count; ++i) {
        print_report_node(nodes[idx].first_child + i, depth + 1, total_reads, 
                          direct_counts, clade_counts, index, out);
    }
}


// Host-side depths for fast LCA
// Struct to hold a batch of reads for processing
struct BatchData {
    std::vector<PackedRead> reads;
    std::vector<std::string> headers1;
    std::vector<std::string> headers2;
    size_t n1;
    std::string base_name;
    bool is_first_in_sample;
    bool is_last_in_sample;
};

// Thread-safe queue for pre-fetching batches
template <typename T>
class ThreadSafeQueue {
    std::queue<T> queue;
    std::mutex mutex;
    std::condition_variable cv_push;
    std::condition_variable cv_pop;
    size_t max_size;
    bool finished = false;

public:
    ThreadSafeQueue(size_t max_sz = 20) : max_size(max_sz) {}

    void push(T item) {
        std::unique_lock<std::mutex> lock(mutex);
        cv_push.wait(lock, [this] { return queue.size() < max_size || finished; });
        if (finished) return;
        queue.push(std::move(item));
        cv_pop.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex);
        cv_pop.wait(lock, [this] { return !queue.empty() || finished; });
        if (queue.empty() && finished) return false;
        item = std::move(queue.front());
        queue.pop();
        cv_push.notify_one();
        return true;
    }

    void finish() {
        std::unique_lock<std::mutex> lock(mutex);
        finished = true;
        cv_pop.notify_all();
        cv_push.notify_all();
    }
};

std::vector<uint32_t> host_depths;

std::string clean_base_name(std::string base_name, bool paired_mode) {
    std::vector<std::string> exts = {".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"};
    for (const auto& ext : exts) {
        if (base_name.size() > ext.size() && 
            base_name.compare(base_name.size() - ext.size(), ext.size(), ext) == 0) {
            base_name.erase(base_name.size() - ext.size());
            break;
        }
    }

    if (paired_mode) {
        std::vector<std::string> r1_markers = {"_R1_001", "_R1", "_1", ".mate1"};
        for (const auto& marker : r1_markers) {
            if (base_name.size() >= marker.size()) {
                size_t pos = base_name.size() - marker.size();
                if (base_name.compare(pos, marker.size(), marker) == 0) {
                    base_name.erase(pos);
                    break;
                }
            }
        }
    }
    return base_name;
}

std::atomic<size_t> global_f_idx{0};

void producer_func(std::vector<std::string> input_files, bool paired_mode, ThreadSafeQueue<std::shared_ptr<BatchData>>* queue, std::atomic<int>* active_count) {
    size_t step = paired_mode ? 2 : 1;
    while (true) {
        size_t f_idx = global_f_idx.fetch_add(step);
        if (f_idx >= input_files.size()) break;

        FastqReader reader1(input_files[f_idx]);
        FastqReader* reader2 = (paired_mode && f_idx + 1 < input_files.size()) ? new FastqReader(input_files[f_idx + 1]) : nullptr;

        if (!reader1.is_open() || (paired_mode && !reader2->is_open())) {
            std::cerr << "Error: Could not open input files: " << input_files[f_idx] << std::endl;
            if (reader2) delete reader2;
            continue;
        }

        std::string base_name = clean_base_name(fs::path(input_files[f_idx]).filename().string(), paired_mode);

        bool is_first = true;
        while (true) {
            auto batch = std::make_shared<BatchData>();
            batch->base_name = base_name;
            batch->is_first_in_sample = is_first;
            batch->is_last_in_sample = false;
            is_first = false;

            batch->n1 = reader1.read_batch(batch->reads, batch->headers1, 50000, 0);
            if (batch->n1 == 0) break;

            if (paired_mode) {
                size_t n2 = reader2->read_batch(batch->reads, batch->headers2, batch->n1, 0);
                if (n2 != batch->n1) {
                    // This is still a bit problematic with parallel producers but good enough for warnings
                }
            }
            queue->push(batch);
        }
        
        auto last_batch = std::make_shared<BatchData>();
        last_batch->base_name = base_name;
        last_batch->is_first_in_sample = false;
        last_batch->is_last_in_sample = true;
        last_batch->n1 = 0;
        queue->push(last_batch);

        if (reader2) delete reader2;
    }
    if (active_count->fetch_sub(1) == 1) {
        queue->finish();
    }
}

// Host-side LCA for merging pass results
uint32_t get_lca_host(uint32_t a, uint32_t b, const KmerIndex& index) {
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == b) return a;

    const auto& nodes = index.getNodes();
    uint32_t da = host_depths[a];
    uint32_t db = host_depths[b];

    while (da > db) {
        a = nodes[a].parent_id;
        da--;
    }
    while (db > da) {
        b = nodes[b].parent_id;
        db--;
    }
    while (a != b) {
        a = nodes[a].parent_id;
        b = nodes[b].parent_id;
    }
    return a;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input1.fastq> [input2.fastq ...] [-d <db_dir>] [-o <output_dir>] [-k <k_length>] [-l true/false] [--paired]" << std::endl;
        return 1;
    }

    std::vector<std::string> input_files;
    std::string output_dir = "output";
    std::string db_dir = "/media/baseripper/volume_28TB_2/k2_standard_08_GB_20260226";
    int k_override = 0;
    bool use_managed = false;
    bool paired_mode = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-o" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "-d" && i + 1 < argc) {
            db_dir = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            k_override = std::stoi(argv[++i]);
        } else if (arg == "-l" && i + 1 < argc) {
            std::string val = argv[++i];
            use_managed = (val == "true");
        } else if (arg == "--paired") {
            paired_mode = true;
        } else {
            input_files.push_back(arg);
        }
    }

    if (!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
    }

    std::cout << "Starting Kraken-GPU Prototype..." << std::endl;

    // 1. Initialize Index
    std::string hash_path = (fs::path(db_dir) / "hash.k2d").string();
    KmerIndex index(0);
    index.loadFromFiles({hash_path}, use_managed); // Load database

    if (index.getHostPtr() == nullptr) {
        std::cerr << "Error: Database could not be loaded. Exiting." << std::endl;
        return 1;
    }
    if (index.getNodes().empty()) {
        std::cerr << "Error: Taxonomy could not be loaded. Exiting." << std::endl;
        return 1;
    }

    index.uploadToDevice(); // Uploads taxonomy structure
    
    // Pre-calculate host depths for LCA merging (if needed)
    const auto& nodes = index.getNodes();
    host_depths.assign(nodes.size(), 0);
    if (!nodes.empty()) {
        host_depths[1] = 1; // Root
        for (size_t i = 2; i < nodes.size(); ++i) {
            uint32_t p = (uint32_t)nodes[i].parent_id;
            if (p < i) host_depths[i] = host_depths[p] + 1;
        }
    }

    Classifier classifier(&index, k_override, 3);
    const std::vector<uint32_t>& taxo_map = index.getTaxonomyMap();

    ThreadSafeQueue<std::shared_ptr<BatchData>> queue(20);
    int num_producers = std::min((int)std::thread::hardware_concurrency(), 16);
    if (input_files.size() < (size_t)num_producers) num_producers = (int)input_files.size() / (paired_mode ? 2 : 1);
    if (num_producers < 1) num_producers = 1;

    std::atomic<int> active_producers{num_producers};
    std::vector<std::thread> producers;
    for (int i = 0; i < num_producers; ++i) {
        producers.emplace_back(producer_func, input_files, paired_mode, &queue, &active_producers);
    }

    struct SampleState {
        std::string base_name;
        std::string out_path;
        std::ofstream out_file;
        uint64_t total_pairs = 0;
        uint64_t total_classified = 0;
        std::vector<uint64_t> direct_counts;
        std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    };
    std::map<std::string, std::shared_ptr<SampleState>> active_samples;

    std::shared_ptr<BatchData> batch;
    while (queue.pop(batch)) {
        std::shared_ptr<SampleState> state;
        if (batch->is_first_in_sample) {
            state = std::make_shared<SampleState>();
            state->base_name = batch->base_name;
            state->out_path = (fs::path(output_dir) / (state->base_name + ".output.txt")).string();
            state->out_file.open(state->out_path);
            state->direct_counts.assign(index.getNodes().size(), 0);
            state->start_time = std::chrono::high_resolution_clock::now();
            active_samples[batch->base_name] = state;
            
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "Starting Sample: " << state->base_name << std::endl;
        } else {
            state = active_samples[batch->base_name];
        }

        if (batch->n1 > 0 && state) {
            size_t n1 = batch->n1;
            std::vector<ClassificationResult> batch_results(batch->reads.size());
            Batch gpu_batch;
            gpu_batch.reads = std::move(batch->reads);
            
            classifier.classifyBatch(gpu_batch, batch_results.data());
            classifier.waitAll();

            for (size_t i = 0; i < n1; ++i) {
                uint32_t internal_id;
                uint32_t read_len;
                const std::string& header = batch->headers1[i];

                if (paired_mode) {
                    uint32_t id1 = batch_results[i].taxid;
                    uint32_t id2 = (i + n1 < batch_results.size()) ? batch_results[i + n1].taxid : 0;
                    internal_id = get_lca_host(id1, id2, index);
                    read_len = gpu_batch.reads[i].length + ((i + n1 < gpu_batch.reads.size()) ? gpu_batch.reads[i + n1].length : 0);
                } else {
                    internal_id = batch_results[i].taxid;
                    read_len = gpu_batch.reads[i].length;
                }

                uint32_t external_id = (internal_id < taxo_map.size()) ? taxo_map[internal_id] : 0;
                
                if (external_id > 0) {
                    state->out_file << "C\t" << header << "\t" << external_id << "\t" << read_len << "\t" << "GPU:1" << "\n";
                    state->total_classified++;
                    if (internal_id < state->direct_counts.size()) state->direct_counts[internal_id]++;
                } else {
                    state->out_file << "U\t" << header << "\t0\t" << read_len << "\t0" << "\n";
                }
            }
            state->total_pairs += n1;
            // Only print progress for one of the active samples to avoid flicker, or just print total
            std::cout << "    [" << state->base_name << "] Processed " << state->total_pairs << " reads...\r" << std::flush;
        }

        if (batch->is_last_in_sample && state) {
            state->out_file.close();

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end_time - state->start_time;

            std::cout << "\nSummary for " << state->base_name << (paired_mode ? " (paired)" : "") << ":" << std::endl;
            std::cout << "  Total Pairs/Reads: " << state->total_pairs << " Classified: " << state->total_classified 
                      << " (" << (state->total_pairs > 0 ? (state->total_classified * 100.0 / state->total_pairs) : 0.0) << "%)" << std::endl;
            std::cout << "  Time: " << diff.count() << "s" << std::endl;

            // Generate Report
            std::string report_path = (fs::path(output_dir) / (state->base_name + ".report.txt")).string();
            std::ofstream report_file(report_path);
            
            size_t num_nodes = nodes.size();
            std::vector<uint64_t> clade_counts = state->direct_counts;
            for (int i = num_nodes - 1; i > 0; --i) {
                uint64_t parent = nodes[i].parent_id;
                if (parent < num_nodes && parent != (uint64_t)i) {
                    clade_counts[parent] += clade_counts[i];
                }
            }

            uint64_t unclassified_count = state->total_pairs - state->total_classified;
            double unclass_pct = (state->total_pairs > 0) ? (unclassified_count * 100.0 / state->total_pairs) : 0.0;
            report_file << std::fixed << std::setprecision(2) << std::setw(6) << unclass_pct << "\t"
                        << unclassified_count << "\t" << unclassified_count << "\tU\t0\tunclassified\n";
            
            uint32_t root_idx = 1;
            for (size_t i = 0; i < num_nodes; ++i) if (nodes[i].external_id == 1) { root_idx = (uint32_t)i; break; }
            print_report_node(root_idx, 0, state->total_pairs, state->direct_counts, clade_counts, index, report_file);
            report_file.close();
            std::cout << "  Report written to " << report_path << std::endl;
            
            active_samples.erase(state->base_name);
        }
    }

    for (auto& t : producers) t.join();
    return 0;
}
