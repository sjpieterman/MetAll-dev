#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <cuda_runtime.h>

typedef uint64_t kmer_t;
typedef uint32_t taxid_t;

const int MAX_READ_LEN = 256;
const int KMER_SIZE = 31;
const int MINIMIZER_WINDOW = 35;

struct PackedRead {
    uint64_t seq[8]; // 256 bp at 2-bits per base
    uint16_t length;
    uint16_t sample_id;
};

struct ClassificationResult {
    taxid_t taxid;
    uint32_t hit_count; 
};

struct Batch {
    std::vector<PackedRead> reads;
    int stream_id;
};

#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
