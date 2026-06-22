# MetaLL/GKraker - pipeline


#MetaLL

Metatranscriptomics pipeline with ultra fast GPU-accelerated speed for large metatranscriptomic datasets.

#GKraker
A GPU-accelerated k-mer classifier inspired by Kraken2/Bracken, designed for high-throughput multi-sample streaming on NVIDIA GPUs (targeted at RTX A5000 24GB).

## Key Features
- **GPU-Resident Index**: Compact hash table for fast k-mer lookups.
- **High-Performance Streaming Reader**: Uses a custom buffered FASTQ reader with native GZIP support to stream reads directly to the GPU.
- **Multi-Threaded I/O**: Parallelized sample loading with up to 16 concurrent producers to saturate GPU classification.
- **Async Streaming**: Uses multiple CUDA streams and a producer-consumer architecture to overlap I/O and computation.
- **2-bit Encoding**: Reads are packed into 2-bit representations to minimize VRAM footprint.
- **Paired-End Support**: Fully supports paired-end classification with host-side LCA merging.
- **Massive Database Support**: Supports databases up to 251GB using CUDA Unified Memory and multi-pass processing.

## Requirements
- Ubuntu 20.04+ (or similar Linux)
- NVIDIA GPU with Compute Capability 8.6+ (Ampere, e.g., RTX A5000)
- CUDA Toolkit 11.0+
- CMake 3.24+
- GCC/G++ 9+

## Building
```bash
mkdir build
cd build
cmake ..
make
```

## Usage
```bash
./kraken_gpu <input1.fastq> [input2.fastq ...] -d <db_dir> [-o <output_dir>] [-k <k_length>] [-l true/false] [--paired]
```
- `-d`: Path to the Kraken2 database directory (containing `hash.k2d`, `opts.k2d`, `taxo.k2d`).
- `-o`: Output directory (default: `output`).
- `-k`: Classification k-mer length (must be >= minimizer length `l`).
- `-l`: Use Unified Managed Memory (default: false). Use `-l true` for databases > 24GB.
- `--paired`: Process files in pairs.

## Implementation Details

### Streaming Pipeline
The system uses a producer-consumer architecture. Up to 16 background threads (producers) decompress and parse FASTQ files concurrently, populating a thread-safe queue with read batches. The main thread (consumer) retrieves these batches and submits them to the GPU via multiple CUDA streams. This ensures the GPU is constantly fed with data even when processing many small samples.

### VRAM Management
With 24GB VRAM, the index is designed to fit ~1.5B k-mer entries (at 12 bytes per entry: 8-byte k-mer + 4-byte TaxID).
Current prototype uses a simpler flat hash table with linear probing.

## To-Do / Future Optimizations
- [ ] Implement BGZF parallel decompression.
- [ ] Add Warp-level primitives for k-mer extraction to reduce divergence.
- [ ] CPU-side Bracken-style abundance estimation.
