#include "classifier.h"
#include <device_launch_parameters.h>
#include <iostream>

// MurmurHash3 fmix64 (Kraken2's hash function)
__device__ uint64_t hash_kmer(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

// Optimized bit-reversal RC matching Kraken2's version 1
__device__ uint64_t reverse_complement_k2(uint64_t kmer, int n, int version) {
    kmer = ((kmer & 0xCCCCCCCCCCCCCCCCULL) >> 2)
         | ((kmer & 0x3333333333333333ULL) << 2);
    kmer = ((kmer & 0xF0F0F0F0F0F0F0F0ULL) >> 4)
         | ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4);
    kmer = ((kmer & 0xFF00FF00FF00FF00ULL) >> 8)
         | ((kmer & 0x00FF00FF00FF00FFULL) << 8);
    kmer = ((kmer & 0xFFFF0000FFFF0000ULL) >> 16)
         | ((kmer & 0x0000FFFF0000FFFFULL) << 16);
    kmer = ( kmer >> 32 ) | ( kmer << 32);
    
    if (version == 0)
        return (~kmer) & ((1ULL << (n * 2)) - 1);
    else
        return ((~kmer) >> (64 - n * 2)) & ((1ULL << (n * 2)) - 1);
}

// Extract minimizer matching Kraken2's logic
__device__ uint64_t get_minimizer_k2(const uint64_t* seq, int start_pos, int k, int l, uint64_t toggle_mask, uint64_t spaced_mask, int rev_version) {
    uint64_t min_candidate = 0xFFFFFFFFFFFFFFFFULL;
    uint64_t min_lmer = 0;
    
    for (int j = 0; j <= k - l; ++j) {
        uint64_t lmer = 0;
        for (int b = 0; b < l; ++b) {
            int pos = start_pos + j + b;
            int uint_idx = pos / 32;
            int bit_offset = (31 - (pos % 32)) * 2;
            uint64_t base = (seq[uint_idx] >> bit_offset) & 0x3;
            lmer = (lmer << 2) | base;
        }
        
        uint64_t rc = reverse_complement_k2(lmer, l, rev_version);
        uint64_t canonical = (lmer < rc) ? lmer : rc;
        
        if (spaced_mask) canonical &= spaced_mask;
        
        uint64_t candidate = canonical ^ toggle_mask;
        if (candidate < min_candidate) {
            min_candidate = candidate;
            min_lmer = canonical;
        }
    }
    return min_lmer;
}

// Compact Hash Table lookup matching Kraken2's Get()
__device__ taxid_t lookup_k2(uint64_t key, KmerEntry* index, size_t capacity, uint32_t tag_bits, uint32_t value_bits) {
    uint64_t hc = hash_kmer(key);
    uint32_t tag = (uint32_t)(hc >> (32 + value_bits));
    size_t idx = hc % capacity;
    size_t first_idx = idx;
    uint64_t step = 0;
    
    for (int i = 0; i < 512; ++i) { // Probe limit for GPU efficiency
        uint32_t cell = index[idx].data;
        if (cell == 0) return 0;
        
        uint32_t cell_tag = cell >> value_bits;
        if (cell_tag == tag) return cell & ((1ULL << value_bits) - 1);
        
        if (step == 0) step = (hc >> 8) | 1;
        idx = (idx + step) % capacity;
        if (idx == first_idx) break;
    }
    return 0;
}

__device__ taxid_t get_lca(taxid_t a, taxid_t b, const uint32_t* parents, const uint32_t* depths) {
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == b) return a;
    uint32_t da = depths[a];
    uint32_t db = depths[b];
    while (da > db) { a = parents[a]; da--; }
    while (db > da) { b = parents[b]; db--; }
    while (a != b) { a = parents[a]; b = parents[b]; }
    return a;
}

__global__ void classify_kernel(PackedRead* reads, int num_reads, 
                               KmerEntry* index, size_t capacity, 
                               uint64_t toggle_mask, uint64_t spaced_mask,
                               uint32_t tag_bits, uint32_t value_bits,
                               int k, int l, int rev_version,
                               const uint32_t* parents, const uint32_t* depths,
                               ClassificationResult* results) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;

    PackedRead read = reads[tid];
    int len = read.length;
    
    taxid_t best_taxid = 0;
    uint32_t max_depth = 0;
    uint32_t hit_count = 0;
    
    if (len >= k) {
        for (int i = 0; i <= len - k; ++i) {
            uint64_t minimizer = get_minimizer_k2(read.seq, i, k, l, toggle_mask, spaced_mask, rev_version);
            taxid_t t = lookup_k2(minimizer, index, capacity, tag_bits, value_bits);
            if (t > 0) {
                uint32_t d = depths[t];
                if (best_taxid == 0) {
                    best_taxid = t;
                    max_depth = d;
                } else {
                    if (d > max_depth) {
                        best_taxid = t;
                        max_depth = d;
                    } else if (d == max_depth) {
                        if (t != best_taxid) {
                            best_taxid = get_lca(t, best_taxid, parents, depths);
                            max_depth = depths[best_taxid];
                        }
                    }
                }
                hit_count++;
            }
        }
    }
    
    results[tid].taxid = best_taxid;
    results[tid].hit_count = hit_count;
}

Classifier::Classifier(KmerIndex* idx, int k_override, int ns) : index(idx), num_streams(ns), current_stream_idx(0) {
    k = (k_override > 0) ? k_override : index->getK();
    l = index->getL();
    toggle_mask = index->getToggleMask();
    spaced_seed_mask = index->getSpacedSeedMask();
    tag_bits = index->getTagBits();
    value_bits = index->getValueBits();
    revcom_version = index->getRevcomVersion();
    d_parents = index->getDeviceParents();
    d_depths = index->getDeviceDepths();
    
    std::cout << "Classifier initialized with k=" << k << ", l=" << l << std::endl;
    
    streams = new cudaStream_t[num_streams];
    for (int i = 0; i < num_streams; ++i) {
        cudaStreamCreate(&streams[i]);
        size_t max_reads = 1000000;
        cudaMalloc(&d_reads[i], max_reads * sizeof(PackedRead));
        cudaMalloc(&d_results[i], max_reads * sizeof(ClassificationResult));
    }
}

Classifier::~Classifier() {
    for (int i = 0; i < num_streams; ++i) {
        cudaStreamDestroy(streams[i]);
        cudaFree(d_reads[i]);
        cudaFree(d_results[i]);
    }
    delete[] streams;
}

void Classifier::classifyBatch(const Batch& batch, ClassificationResult* h_results) {
    int s_idx = current_stream_idx;
    int num_reads = batch.reads.size();
    
    cudaMemcpyAsync(d_reads[s_idx], batch.reads.data(), num_reads * sizeof(PackedRead),
                    cudaMemcpyHostToDevice, streams[s_idx]);
    
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_reads + threadsPerBlock - 1) / threadsPerBlock;
    classify_kernel<<<blocksPerGrid, threadsPerBlock, 0, streams[s_idx]>>>(
        d_reads[s_idx], num_reads, index->getHostPtr(), index->getCapacity(), 
        toggle_mask, spaced_seed_mask, tag_bits, value_bits, k, l, revcom_version, 
        d_parents, d_depths, d_results[s_idx]);
    
    cudaMemcpyAsync(h_results, d_results[s_idx], num_reads * sizeof(ClassificationResult),
                    cudaMemcpyDeviceToHost, streams[s_idx]);
    
    current_stream_idx = (current_stream_idx + 1) % num_streams;
}

void Classifier::waitAll() {
    for (int i = 0; i < num_streams; ++i) {
        cudaStreamSynchronize(streams[i]);
    }
}
