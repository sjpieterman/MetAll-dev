#include "index.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

KmerIndex::KmerIndex(size_t cap) : capacity(cap), k(35), l(31), toggle_mask(0), spaced_seed_mask(0), tag_bits(0), value_bits(0), revcom_version(1), h_entries(nullptr), d_entries(nullptr), d_parents(nullptr), d_depths(nullptr) {
}

KmerIndex::~KmerIndex() {
    if (h_entries) {
        cudaFree(h_entries);
    }
    if (d_parents) cudaFree(d_parents);
    if (d_depths) cudaFree(d_depths);
}

void KmerIndex::uploadToDevice() {
    // Note: d_entries is now handled by uploadChunk for multi-pass
    // or manually if capacity fits in VRAM.
    
    if (!nodes.empty()) {
        std::vector<uint32_t> h_parents(nodes.size());
        std::vector<uint32_t> h_depths(nodes.size());
        
        // Calculate depths on host
        for (size_t i = 0; i < nodes.size(); ++i) {
            h_parents[i] = (uint32_t)nodes[i].parent_id;
        }

        // root is node 1 usually, node 0 is unclassified
        h_depths[0] = 0; // unclassified
        h_depths[1] = 1; // root
        for (size_t i = 2; i < nodes.size(); ++i) {
            uint32_t p = h_parents[i];
            if (p < i) {
                h_depths[i] = h_depths[p] + 1;
            } else {
                // Should not happen in BFS order, but for safety:
                uint32_t d = 0;
                uint32_t curr = i;
                while (curr != h_parents[curr] && d < 100) {
                    curr = h_parents[curr];
                    d++;
                }
                h_depths[i] = d;
            }
        }

        CUDA_CHECK(cudaMalloc(&d_parents, nodes.size() * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemcpy(d_parents, h_parents.data(), nodes.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
        
        CUDA_CHECK(cudaMalloc(&d_depths, nodes.size() * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemcpy(d_depths, h_depths.data(), nodes.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    }
}

void KmerIndex::uploadChunk(size_t start_idx, size_t count) {
    if (!h_entries) return;
    if (start_idx + count > capacity) count = capacity - start_idx;
    
    if (!d_entries) {
        CUDA_CHECK(cudaMalloc(&d_entries, count * sizeof(KmerEntry)));
    }
    CUDA_CHECK(cudaMemcpy(d_entries, h_entries + start_idx, count * sizeof(KmerEntry), cudaMemcpyHostToDevice));
}

void KmerIndex::loadFromFiles(const std::vector<std::string>& files, bool use_managed) {
    if (files.empty()) {
        std::cerr << "No files provided for database loading." << std::endl;
        return;
    }

    std::string hash_file = files[0];
    std::string base_dir = hash_file.substr(0, hash_file.find_last_of("/\\") + 1);
    std::string opts_file = base_dir + "opts.k2d";

    // 1. Load Options (k, l, Spaced Mask, Toggle Mask, dna_db, min_hash, revcom_version)
    std::ifstream opts_f(opts_file, std::ios::binary);
    if (opts_f) {
        uint64_t opts_header[8];
        opts_f.read((char*)opts_header, 64);
        k = (uint32_t)opts_header[0];
        l = (uint32_t)opts_header[1];
        spaced_seed_mask = opts_header[2];
        toggle_mask = opts_header[3];
        revcom_version = (int)opts_header[6];
        std::cout << "Loaded from opts.k2d: k=" << k << ", l=" << l 
                  << ", Toggle mask=0x" << std::hex << toggle_mask 
                  << ", Spaced mask=0x" << spaced_seed_mask << std::dec 
                  << ", Revcom v" << revcom_version << std::endl;
    } else {
        std::cerr << "Warning: Could not open opts.k2d. Using zero mask." << std::endl;
    }

    // 2. Load Hash Table
    std::cout << "Loading Kraken2 database from: " << hash_file << std::endl;

    int fd = open(hash_file.c_str(), O_RDONLY);
    if (fd == -1) {
        perror("Error opening hash.k2d");
        return;
    }

    struct stat st;
    if (fstat(fd, &st) == -1) {
        perror("Error stating hash.k2d");
        close(fd);
        return;
    }
    size_t file_size = st.st_size;

    // Read header (at least 32 bytes)
    uint64_t header[4];
    if (read(fd, header, 32) != 32) {
        std::cerr << "Error reading header from " << hash_file << std::endl;
        close(fd);
        return;
    }

    // Kraken2 hash.k2d structure:
    // [8 bytes capacity] [8 bytes count] [8 bytes tag_bits] [8 bytes value_bits]
    size_t db_capacity = header[0];
    tag_bits = (uint32_t)header[2];
    value_bits = (uint32_t)header[3];

    std::cout << "Database reported capacity: " << db_capacity << " slots" << std::endl;
    std::cout << "Bits: tag=" << tag_bits << ", value=" << value_bits << std::endl;

    // Adjust capacity based on what we can handle
    if (capacity == 0 || capacity > db_capacity) {
        capacity = db_capacity;
    }

    if (use_managed) {
        std::cout << "Allocating " << (capacity * sizeof(KmerEntry)) / (1024*1024) << " MB (Managed Memory)..." << std::endl;
        CUDA_CHECK(cudaMallocManaged(&h_entries, capacity * sizeof(KmerEntry)));
    } else {
        std::cout << "Allocating " << (capacity * sizeof(KmerEntry)) / (1024*1024) << " MB (VRAM)..." << std::endl;
        CUDA_CHECK(cudaMalloc(&h_entries, capacity * sizeof(KmerEntry)));
    }

    // Use mmap to read from the file efficiently
    void* mapped = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        perror("mmap failed, falling back to read()");
        lseek(fd, 32, SEEK_SET);
        if (use_managed) {
            size_t entries_read = read(fd, h_entries, capacity * sizeof(KmerEntry)) / sizeof(KmerEntry);
            capacity = entries_read;
        } else {
            // Need temporary host buffer if not managed
            KmerEntry* tmp_h = new KmerEntry[capacity];
            size_t entries_read = read(fd, tmp_h, capacity * sizeof(KmerEntry)) / sizeof(KmerEntry);
            capacity = entries_read;
            CUDA_CHECK(cudaMemcpy(h_entries, tmp_h, capacity * sizeof(KmerEntry), cudaMemcpyHostToDevice));
            delete[] tmp_h;
        }
    } else {
        // Kraken2 entries start at offset 32
        KmerEntry* entries_ptr = (KmerEntry*)((char*)mapped + 32);
        if (use_managed) {
            memcpy(h_entries, entries_ptr, capacity * sizeof(KmerEntry));
        } else {
            CUDA_CHECK(cudaMemcpy(h_entries, entries_ptr, capacity * sizeof(KmerEntry), cudaMemcpyHostToDevice));
        }
        munmap(mapped, file_size);
    }

    std::cout << "Successfully loaded " << capacity << " entries (" << (capacity * sizeof(KmerEntry)) / (1024*1024) << " MB)." << std::endl;
    close(fd);

    // 3. Load Taxonomy
    std::string taxo_file = base_dir + "taxo.k2d";
    std::ifstream taxo_f(taxo_file, std::ios::binary);
    if (taxo_f) {
        uint64_t taxo_header[4];
        taxo_f.read((char*)taxo_header, 32);
        size_t num_nodes = taxo_header[1];
        size_t name_data_len = taxo_header[2];
        size_t rank_data_len = taxo_header[3];
        
        nodes.resize(num_nodes);
        taxo_f.read((char*)nodes.data(), num_nodes * sizeof(TaxonomyNode));
        
        name_data.resize(name_data_len);
        taxo_f.read(name_data.data(), name_data_len);
        
        rank_data.resize(rank_data_len);
        taxo_f.read(rank_data.data(), rank_data_len);

        // Fill legacy taxonomy_map for compatibility
        taxonomy_map.resize(num_nodes, 0);
        for (size_t i = 0; i < num_nodes; ++i) {
            taxonomy_map[i] = (uint32_t)nodes[i].external_id;
        }
        
        std::cout << "Loaded taxonomy: " << num_nodes << " nodes, " 
                  << name_data_len << " bytes names, " << rank_data_len << " bytes ranks." << std::endl;
    } else {
        std::cerr << "Warning: Could not open taxo.k2d." << std::endl;
    }
}
