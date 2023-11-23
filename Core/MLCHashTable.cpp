//
// Created by ldd on 2023/3/7.
//

#include "MLCHashTable.h"


MLCHashTable::MLCHashTable(uint n_, uint dim_, uint bucket_size_) : n(n_), dim(dim_), bucket_size(bucket_size_), n_mlcs(0) {
    k_size = dim * sizeof(uint);
    hash_node_size = sizeof(MLCHashNode);

    bucket = new MLCHashNode *[bucket_size];

    for (uint i = 0; i < bucket_size; i++) bucket[i] = nullptr;
}

MLCHashTable::~MLCHashTable() {
    delete[] bucket;

    mlc_node_buf.Release();
    hash_node_buf.Release();
}

ll_uint MLCHashTable::CountDistinct() {
    uint *mlc;
    Trie trie;

    MLCHashTable::iterator iter = GetIterator();

    while (true) {
        mlc = iter.mlc;
        sort(mlc + 1, mlc + 1 + *mlc);
        trie.Add(mlc + 1, *mlc);
        if (!GetNext(iter)) {
            break;
        }
    }

    return trie.GetCount();
}

void MLCHashTable::Insert(uint *k, uint *mlc, uint length) {
    uint buk_id, *ml_core;
    MLCNode *mlc_node;
    MLCHashNode *new_hash_node;

    mlc_node = AllocateMLCNode(length);
    memcpy(mlc_node, k, k_size);

    ml_core = get_mlc_addr(mlc_node);
    ml_core[0] = length;
    memcpy(ml_core + 1, mlc, length * sizeof(uint));

    buk_id = range_hash_k(k);
    new_hash_node = AllocateHashNode();
    new_hash_node->MLC_node = mlc_node;
    new_hash_node->next_node = bucket[buk_id];
    bucket[buk_id] = new_hash_node;

    n_mlcs++;
}

uint MLCHashTable::Search(vector<uint> &k_vec, uint *mlc) {
    uint *k, *core;
    MLCNode *mlc_node;
    MLCHashNode *hash_node;

    hash_node = bucket[range_hash_k(k_vec.data())];
    while (hash_node) {
        mlc_node = hash_node->MLC_node;
        k = get_k_vec_addr(mlc_node);
        if (!memcmp(k_vec.data(), k, k_size)) {
            core = get_mlc_addr(mlc_node);
            memcpy(mlc, core + 1, core[0] * sizeof(uint));
            return core[0];
        } else hash_node = hash_node->next_node;
    }
    return 0;
}

void MLCHashTable::Flush(const string &file) {
    uint length;
    MLCNode *mlc_node;
    MLCHashNode *hash_node;

    auto f = std::ofstream(file);

    f.write(reinterpret_cast<char *> (&n), sizeof(uint));
    f.write(reinterpret_cast<char *> (&dim), sizeof(uint));
    f.write(reinterpret_cast<char *> (&bucket_size), sizeof(uint));
    f.write(reinterpret_cast<char *> (&n_mlcs), sizeof(uint));

    for (uint i = 0; i < bucket_size; i++) {
        hash_node = bucket[i];
        while (hash_node) {
            mlc_node = hash_node->MLC_node;
            length = *(get_mlc_addr(mlc_node));
            f.write(reinterpret_cast<char *> (mlc_node), (long) ((dim + length + 1) * sizeof(uint)));
            hash_node = hash_node->next_node;
        }
    }

    f.close();
}

MLCHashTable *MLCHashTable::Load(const string &file) {
    uint num, dimension, buk_size;

    auto f = std::ifstream(file);

    f.read(reinterpret_cast<char *> (&num), sizeof(uint));
    f.read(reinterpret_cast<char *> (&dimension), sizeof(uint));
    f.read(reinterpret_cast<char *> (&buk_size), sizeof(uint));

    auto *ht = new MLCHashTable(num, dimension, buk_size);
    ht->Load(f);
    f.close();
    return ht;
}

void MLCHashTable::Load(std::basic_ifstream<char> &f) {
    uint k[dim], buk_id, *core, length;
    MLCNode *mlc_node;
    MLCHashNode *hash_node;

    f.read(reinterpret_cast<char *> (&n_mlcs), sizeof(uint));
    for (uint i = 0; i < n_mlcs; i++) {
        f.read(reinterpret_cast<char *> (&k[0]), long(k_size));
        f.read(reinterpret_cast<char *> (&length), sizeof(uint));

        mlc_node = AllocateMLCNode(length);
        memcpy(mlc_node, k, k_size);
        core = get_mlc_addr(mlc_node);
        core[0] = length;
        f.read(reinterpret_cast<char *> (core + 1), long(length * sizeof(uint)));

        buk_id = range_hash_k(k);
        hash_node = AllocateHashNode();
        hash_node->MLC_node = mlc_node;
        hash_node->next_node = bucket[buk_id];
        bucket[buk_id] = hash_node;
    }
}

bool MLCHashTable::Equals(MLCHashTable &mlc_ht) {
    uint size, *core, *mlc_ht_core;
    MLCHashNode *hash_node;

    auto MLC_cmp = [this](MLCNode *MLC1, MLCNode *MLC2) {
        uint *k1 = get_k_vec_addr(MLC1), *k2 = get_k_vec_addr(MLC2);
        for (uint i = 0; i < dim; i++) {
            if (k1[i] < k2[i]) return true;
            else if (k1[i] > k2[i]) return false;
        }
        return false;
    };

    // Equals basic elements
    if (dim != mlc_ht.dim || bucket_size != mlc_ht.bucket_size || n_mlcs != mlc_ht.n_mlcs) {
        return false;
    }

    for (uint i = 0; i < bucket_size; i++) {
        vector<MLCNode *> mlc_nodes, mlc_ht_mlc_nodes;

        hash_node = bucket[i];
        while (hash_node) {
            mlc_nodes.emplace_back(hash_node->MLC_node);
            hash_node = hash_node->next_node;
        }

        hash_node = mlc_ht.bucket[i];
        while (hash_node) {
            mlc_ht_mlc_nodes.emplace_back(hash_node->MLC_node);
            hash_node = hash_node->next_node;
        }

        if (mlc_nodes.size() != mlc_ht_mlc_nodes.size()) return false;
        if (mlc_nodes.empty()) continue;

        // sort and compare
        sort(mlc_nodes.begin(), mlc_nodes.end(), MLC_cmp);
        sort(mlc_ht_mlc_nodes.begin(), mlc_ht_mlc_nodes.end(), MLC_cmp);

        size = (uint) mlc_nodes.size();
        for (uint j = 0; j < size; j++) {
            if (memcmp(mlc_nodes[j], mlc_ht_mlc_nodes[j], k_size) != 0) return false;

            core = get_mlc_addr(mlc_nodes[j]);
            mlc_ht_core = get_mlc_addr(mlc_ht_mlc_nodes[j]);

            if (core[0] != mlc_ht_core[0]) return false;

            sort(core + 1, core + 1 + core[0]);
            sort(mlc_ht_core + 1, mlc_ht_core + 1 + core[0]);

            if (memcmp(core + 1, mlc_ht_core + 1, core[0] * sizeof(uint)) != 0) return false;
        }
    }
    return true;
}

bool MLCHashTable::GetNext(iterator &iter) {
    MLCHashNode *new_mlc_hash_node;

    if (iter.mlc_hash_node->next_node) {
        new_mlc_hash_node = iter.mlc_hash_node->next_node;
        iter = {iter.bucket_id, get_mlc_addr(new_mlc_hash_node->MLC_node), new_mlc_hash_node};
        return true;

    } else {
        auto new_bid = GetNextValidBuk(iter.bucket_id + 1);

        if (new_bid != -1) {
            new_mlc_hash_node = bucket[new_bid];
            iter = {new_bid, get_mlc_addr(new_mlc_hash_node->MLC_node), new_mlc_hash_node};
            return true;
        }
    }
    return false;
}

void MLCHashTable::PrintLayout() {
    MLCHashNode *hash_node;
    ll_uint n_nodes = 0;

    for (int i = 0; i < bucket_size; i++) {
        if (bucket[i]) {
            cout << "Bucket " << i << ": ";
            hash_node = bucket[i];
            while (hash_node) {
                cout << Arr2Str(get_k_vec_addr(hash_node->MLC_node), dim) << " ";
                n_nodes += 1;
                hash_node = hash_node->next_node;
            }
            cout << endl;
        }
    }

    assert(n_nodes == n_mlcs);
    cout << "Total " << n_nodes << " mlcs." << endl;
}