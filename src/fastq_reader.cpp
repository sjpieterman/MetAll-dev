#include "fastq_reader.h"
#include <cstring>
#include <algorithm>
#include <iostream>

FastqReader::FastqReader(const std::string& path) 
    : buffer_pos(0), bytes_in_buffer(0), eof_reached(false) {
    fp = gzopen(path.c_str(), "rb");
    if (fp) {
        gzbuffer(fp, 128 * 1024); // Internal zlib buffer
        buffer = new char[BUFFER_SIZE];
    } else {
        buffer = nullptr;
    }
}

FastqReader::~FastqReader() {
    if (fp) gzclose(fp);
    if (buffer) delete[] buffer;
}

void FastqReader::fill_buffer() {
    if (eof_reached) return;
    
    int n = gzread(fp, buffer, BUFFER_SIZE);
    if (n <= 0) {
        eof_reached = true;
        bytes_in_buffer = 0;
    } else {
        bytes_in_buffer = (size_t)n;
    }
    buffer_pos = 0;
}

bool FastqReader::read_line(std::string& line) {
    line.clear();
    while (true) {
        if (buffer_pos >= bytes_in_buffer) {
            fill_buffer();
            if (bytes_in_buffer == 0) return !line.empty();
        }

        // Find newline in remaining buffer
        char* start = buffer + buffer_pos;
        char* end = (char*)memchr(start, '\n', bytes_in_buffer - buffer_pos);
        
        if (end) {
            size_t len = end - start;
            line.append(start, len);
            buffer_pos += len + 1; // skip newline
            break;
        } else {
            // No newline in current buffer
            line.append(start, bytes_in_buffer - buffer_pos);
            buffer_pos = bytes_in_buffer;
        }
    }

    // Trim \r
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return true;
}

bool FastqReader::next_read(std::string& header, std::string& seq) {
    if (!is_open()) return false;

    std::string line;
    // 1. Header
    if (!read_line(line)) return false;
    if (line.empty() || line[0] != '@') return false;
    
    header = line.substr(1);
    size_t space_pos = header.find(' ');
    if (space_pos != std::string::npos) {
        header = header.substr(0, space_pos);
    }

    // 2. Sequence
    if (!read_line(seq)) return false;

    // 3. +
    if (!read_line(line)) return false;

    // 4. Quality
    if (!read_line(line)) return false;

    return true;
}

void pack_read_internal(const std::string& seq, PackedRead& read, uint16_t sample_id) {
    memset(&read, 0, sizeof(PackedRead));
    read.length = (uint16_t)std::min((size_t)MAX_READ_LEN, seq.length());
    read.sample_id = sample_id;

    for (int i = 0; i < read.length; ++i) {
        uint64_t base;
        switch (seq[i]) {
            case 'A': case 'a': base = 0; break;
            case 'C': case 'c': base = 1; break;
            case 'G': case 'g': base = 2; break;
            case 'T': case 't': base = 3; break;
            default: base = 0; break;
        }
        int uint_idx = i / 32;
        int bit_offset = (31 - (i % 32)) * 2; 
        read.seq[uint_idx] |= (base << bit_offset);
    }
}

size_t FastqReader::read_batch(std::vector<PackedRead>& reads, std::vector<std::string>& headers, size_t max_reads, uint16_t sample_id) {
    size_t count = 0;
    std::string h, s;
    while (count < max_reads && next_read(h, s)) {
        headers.push_back(h);
        PackedRead pr;
        pack_read_internal(s, pr, sample_id);
        reads.push_back(pr);
        count++;
    }
    return count;
}
