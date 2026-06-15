#pragma once

#include "common.h"
#include "index.h"

class Classifier {
public:
    Classifier(KmerIndex* index, int k_override = 0, int num_streams = 3);
    ~Classifier();

    void classifyBatch(const Batch& batch, ClassificationResult* h_results);
    void waitAll();

private:
    KmerIndex* index;
    int num_streams;
    int k, l;
    cudaStream_t* streams;
    uint64_t toggle_mask;
    uint64_t spaced_seed_mask;
    uint32_t tag_bits;
    uint32_t value_bits;
    int revcom_version;
    
    uint32_t* d_parents;
    uint32_t* d_depths;
    
    // Staging buffers for each stream
    PackedRead* d_reads[3];
    ClassificationResult* d_results[3];
    
    int current_stream_idx;
};
