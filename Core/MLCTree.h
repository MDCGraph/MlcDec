//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_MLCTREE_H
#define MLCDEC_MLCTREE_H


#include "../Structures/DataBuf.h"
#include "../Utils/ArrayUtils.h"

/*
 *  node alignment: | === k === | IncD(k) | Diff |  chd_l |  chd_2 | ...
 */

struct Node;

struct Diff {
    uint num;
    uint *vtx_ptr;
};

class MLCTree {
public:
    explicit MLCTree(uint dim_, uint n_, bool is_isolate = true);
    ~MLCTree();

    [[nodiscard]] ll_uint GetNumOfNodes() const {
        return num_of_nodes;
    }

    void SetNumOfNodes(ll_uint num_of_nodes_) {
        num_of_nodes = num_of_nodes_;
    }

    void SetRoot(Node *r) {
        root = r;
    }

    Node *GetRoot() {
        return root;
    }

    [[nodiscard]] bool Isolate() const{
        return isolate;
    }

    void SetIsolate() {
        isolate = true;
    }

    [[nodiscard]] uint GetLayerNumber() const {
        return dim;
    }

    [[nodiscard]] uint GetN() const {
        return n;
    }

    inline Node *GetNewTreeNode() {
        Node *new_node;

        new_node = Allocate(0);
        memset(new_node, 0, i_off);
        *GetIncKAddr(new_node) = 0;

        num_of_nodes++;
        return new_node;
    }

    inline Node *GetNewTreeNode(uint *k, uint inc_k) {
        Node *new_node;

        new_node = Allocate(inc_k);
        memcpy(new_node, k, i_off);
        *GetIncKAddr(new_node) = inc_k;

        num_of_nodes++;
        return new_node;
    }

    uint* AllocateDiff(uint size) {
        return difference_buf.Allocate(size);
    }

    static inline uint *GetK(Node *node) {
        return reinterpret_cast<uint *>(node);
    }

    inline uint GetIncK(Node *node) const {
        return *reinterpret_cast<uint *>(reinterpret_cast<char *> (node) + i_off);
    }

    inline Node *GetRelChd(Node *node, uint i) const {
        return *(reinterpret_cast<Node **>(reinterpret_cast<char *> (node) + c_off + (i - GetIncK(node)) * sizeof(Node *)));
    }

    inline void SetRelChd(Node *node, uint i, Node *chd) const {
        *(reinterpret_cast<Node **>(reinterpret_cast<char *> (node) + c_off + (i - GetIncK(node)) * sizeof(Node *))) = chd;
    }

    inline Node *GetRmChd(Node *node) const {
        return *(reinterpret_cast<Node **>(reinterpret_cast<char *> (node) + c_off +
        (dim - GetIncK(node) - 1) * sizeof(Node *)));
    }

    inline Diff *GetDiff(Node *node) const {
        return reinterpret_cast<Diff *>(reinterpret_cast<char *> (node) + d_off);
    }

//    inline void SetDiff(Node *node, uint length, uint *vtx = nullptr) {
//        Diff *d = GetDiff(node);
//        d->num = length;
//        if (length) {
//            d->vtx_ptr = difference_buf.Allocate(length);
//            memcpy(d->vtx_ptr, vtx, length * sizeof(uint));
//        }
//    }
//
    inline void SetDiff(Node *node, uint length, uint *vtx, const uint * relabel = nullptr) {
        Diff *d = GetDiff(node);
        d->num = length;
        if (length) {
            d->vtx_ptr = difference_buf.Allocate(length);
            if (!relabel) memcpy(d->vtx_ptr, vtx, length * sizeof(uint));
            else {
                for (uint i = 0; i < length; i++) {
                    d->vtx_ptr[i] = relabel[vtx[i]];
                }
            }
        }
    }

    uint Search(const vector<uint> &k_vec, uint *core, uint &length);
    uint Retrieve(Node *node, uint *core) const;

    /* test parallel version */
    uint PSearch(const vector<uint> &k_vec, uint *core, uint &length);
    uint PRetrieve(Node *node, uint *core) const;
    uint Fd(Node* node, uint* core) const;

    void Flush(const string &file);
    void Flush(std::basic_ofstream<char> &f);

    static MLCTree *Load(const string &file);
    void Load(std::basic_ifstream<char> &f);

    bool Equals(MLCTree &mlc_t);
    void MergeBuf(MLCTree *mlc_t);

private:

    MDataBuf<char> tree_node_buf;
    DataBuf<uint> difference_buf;

    ll_uint num_of_nodes;

    Node *root;
    uint dim;
    uint n;

    size_t i_off;  // offset of the value of inc_k
    size_t d_off;  // offset of difference
    size_t c_off;  // offset of branches

    bool isolate;

    inline Node *Allocate(uint inc_k) {
        return reinterpret_cast<Node *>(tree_node_buf.Allocate(c_off + (dim - inc_k) * sizeof(Node *)));
    }

    inline uint *GetIncKAddr(Node *node) const {
        return reinterpret_cast<uint *>(reinterpret_cast<char *> (node) + i_off);
    }

    inline Node *GetLmChd(Node *node) const {
        return *(reinterpret_cast<Node **>(reinterpret_cast<char *> (node) + c_off));
    }

    bool IsEqual(Node *node, Node *mlc_tree_node);

    void FlushNode(Node *node, unordered_map<Node *, ll_uint> &node2id, std::basic_ofstream<char> &f) const;
    void LoadNode(Node *node, unordered_map<ll_uint, Node *> &id2node, std::basic_ifstream<char> &f);
};


#endif //MLCDEC_MLCTREE_H
