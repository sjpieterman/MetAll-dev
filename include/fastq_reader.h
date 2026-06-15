#pragma once

#include <string>
#include <vector>
#include <zlib.h>
#include "common.h"

class FastqReader {
public:
    FastqReader(const std::string& path);
    ~FastqReader();

    bool is_open() const { return fp != nullptr; }
    
    /**
     * @brief Reads next sequence from FASTQ file.
     * @param header Output for read ID (everything before first space).
     * @param seq Output for sequence string.
     * @return true if successful, false on EOF or error.
     */
    bool next_read(std::string& header, std::string& seq);

    /**
     * @brief Reads a batch of sequences and packs them.
     * @param reads Output vector to append packed reads.
     * @param headers Output vector to append headers.
     * @param max_reads Maximum number of reads to append.
     * @param sample_id ID to assign to packed reads.
     * @return Number of reads actually appended.
     */
    size_t read_batch(std::vector<PackedRead>& reads, std::vector<std::string>& headers, size_t max_reads, uint16_t sample_id);

private:
    gzFile fp;
    char* buffer;
    static const size_t BUFFER_SIZE = 1024 * 1024; // 1MB buffer
    size_t buffer_pos;
    size_t bytes_in_buffer;
    bool eof_reached;

    void fill_buffer();
    bool read_line(std::string& line);
};

void pack_read_internal(const std::string& seq, PackedRead& read, uint16_t sample_id);
