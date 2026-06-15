#pragma once

#include "common.h"

struct KmerEntry {
    uint32_t data;
} __attribute__((packed));

struct TaxonomyNode {
    uint64_t parent_id;
    uint64_t first_child;
    uint64_t child_count;
    uint64_t name_offset;
    uint64_t rank_offset;
    uint64_t external_id;
    uint64_t godparent_id;
};

class KmerIndex {
public:
    KmerIndex(size_t capacity);
    ~KmerIndex();

    // Host methods
    void loadFromFiles(const std::vector<std::string>& files, bool use_managed = false);
    void uploadToDevice();
    void uploadChunk(size_t start_idx, size_t count);

    // Device access
    KmerEntry* getDevicePtr() { return d_entries; }
    KmerEntry* getHostPtr() { return h_entries; }
    uint32_t* getDeviceParents() { return d_parents; }
    uint32_t* getDeviceDepths() { return d_depths; }
    size_t getCapacity() const { return capacity; }
    uint64_t getToggleMask() const { return toggle_mask; }
    uint64_t getSpacedSeedMask() const { return spaced_seed_mask; }
    uint32_t getK() const { return k; }
    uint32_t getL() const { return l; }
    uint32_t getTagBits() const { return tag_bits; }
    uint32_t getValueBits() const { return value_bits; }
    int getRevcomVersion() const { return revcom_version; }
    
    // Taxonomy access
    const std::vector<TaxonomyNode>& getNodes() const { return nodes; }
    const std::vector<char>& getNameData() const { return name_data; }
    const std::vector<char>& getRankData() const { return rank_data; }
    const std::vector<uint32_t>& getTaxonomyMap() const { return taxonomy_map; }

private:
    size_t capacity;
    uint32_t k;
    uint32_t l;
    uint64_t toggle_mask;
    uint64_t spaced_seed_mask;
    uint32_t tag_bits;
    uint32_t value_bits;
    int revcom_version;
    
    std::vector<TaxonomyNode> nodes;
    std::vector<char> name_data;
    std::vector<char> rank_data;
    std::vector<uint32_t> taxonomy_map; // Mapping of internal ID to external ID (compatibility)

    KmerEntry* h_entries;
    KmerEntry* d_entries;
    uint32_t* d_parents;
    uint32_t* d_depths;
};
