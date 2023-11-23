//
// Created by ldd on 2023/3/7.
//

#ifndef MLCDEC_MLCHASHTABLE_H
#define MLCDEC_MLCHASHTABLE_H

#include "../Structures/DataBuf.h"
#include "../Structures/Trie.h"

#include "../Utils/ArrayUtils.h"

const uint HASH_BUCKET_SIZE = (1 << 10); // 1024

/*
 * MLCNode alignment: | === k === | length | === mlc_core === |
 */

struct MLCNode;

struct MLCHashNode {
    MLCNode *MLC_node{};
    MLCHashNode *next_node{nullptr};
};

class MLCHashTable {
public:
    explicit MLCHashTable(uint n, uint dim_, uint bucket_size_ = HASH_BUCKET_SIZE);

    ~MLCHashTable();

    [[nodiscard]] ll_uint GetNumOfMLC() const {
        return n_mlcs;
    }


    ll_uint CountDistinct();

    void Insert(uint *k, uint *mlc, uint length);

    uint Search(vector<uint> &k_vec, uint *mlc);

    void Flush(const string &file);

    static MLCHashTable *Load(const string &file);

    bool Equals(MLCHashTable &mlc_ht);

    [[nodiscard]] uint GetN() const {
        return n;
    }

    void PrintLayout();

    // ========== iterator ============
    struct iterator {
        uint bucket_id;
        uint *mlc;
        MLCHashNode *mlc_hash_node;
    };

    iterator GetIterator() {
        auto bid = GetNextValidBuk(0);
        auto mlc_hash_node = bucket[bid];

        return {bid, get_mlc_addr(mlc_hash_node->MLC_node), mlc_hash_node};
    }

    bool GetNext(iterator &iter);

    inline uint GetNextValidBuk(uint bid) {

        while (bid < bucket_size) {
            if (bucket[bid]) return bid;
            bid += 1;
        }
        return -1;
    }

    // ========== iterator ============

private:
    MDataBuf<char> hash_node_buf;
    DataBuf<uint> mlc_node_buf;

    MLCHashNode **bucket;

    uint dim;
    uint n;
    uint bucket_size;
    ll_uint n_mlcs;

    size_t k_size;
    size_t hash_node_size;

    // boost/functional/hash/hash.hpp
    inline uint range_hash_k(const uint *k) const {
        std::size_t seed = 0;
        hash<uint> hash;

        for (uint i = 0; i < dim; i++) {
            seed ^= hash(k[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        return (uint) (seed % bucket_size);
    }

    inline MLCHashNode *AllocateHashNode() {
        return reinterpret_cast<MLCHashNode *>(hash_node_buf.Allocate(hash_node_size));
    }

    inline MLCNode *AllocateMLCNode(uint length) {
        return reinterpret_cast<MLCNode *>(mlc_node_buf.Allocate(dim + length + 1));
    }

    static inline uint *get_k_vec_addr(MLCNode *node) {
        return reinterpret_cast<uint *> (node);
    }

    inline uint *get_mlc_addr(MLCNode *node) const {
        return reinterpret_cast<uint *> (reinterpret_cast<char *> (node) + k_size);
    }

    void Load(std::basic_ifstream<char> &f);
};


#endif //MLCDEC_MLCHASHTABLE_H
