//
// Created by ldd on 2023/3/6.
//

#include "MLCTree.h"

MLCTree::MLCTree(uint dim_, uint n_, bool is_isolate) : root(nullptr), dim(dim_), n(n_), num_of_nodes(0),
                                                        isolate(is_isolate) {
    i_off = dim_ * sizeof(uint);
    d_off = i_off + sizeof(uint);
    c_off = d_off + sizeof(struct Diff);
}

MLCTree::~MLCTree() {
    if (isolate) {
        tree_node_buf.Release();
        difference_buf.Release();
    }
}

uint MLCTree::Search(const vector<uint> &k_vec, uint *core, uint &length) {
    uint i, num_of_visits;
    Node *node;

    i = 0;
    node = root;
    num_of_visits = 1;  // visiting root

    // locate (k,p)-node
    while (true) {
        while (node && GetK(node)[i] < k_vec[i]) {
            node = GetLmChd(node);
            num_of_visits++;
        }
        if (!node) break;
        else {
            i++;
            while (i < dim && GetK(node)[i] >= k_vec[i]) i++;
            if (i == dim) break;
            else {
                node = GetRelChd(node, i);
                num_of_visits++;
            }
        }
    }

    length = Retrieve(node, core);
    return num_of_visits;
}

uint MLCTree::Retrieve(Node *node, uint *core) const {
    uint length = 0;
    Diff *diff;

    while (node) {
        diff = GetDiff(node);  // collect vertices along the right-most path.
        memcpy(core + length, diff->vtx_ptr, diff->num * sizeof(uint));
        length += diff->num;
        node = GetRmChd(node);
    }

    return length;
}

uint MLCTree::PSearch(const vector<uint> &k_vec, uint *core, uint &length) {
    uint i, num_of_visits;
    Node *node;

    i = 0;
    node = root;
    num_of_visits = 1;  // visiting root

    // locate (k,p)-node
    while (true) {
        while (node && GetK(node)[i] < k_vec[i]) {
            node = GetLmChd(node);
            num_of_visits++;
        }
        if (!node) break;
        else {
            i++;
            while (i < dim && GetK(node)[i] >= k_vec[i]) i++;
            if (i == dim) break;
            else {
                node = GetRelChd(node, i);
                num_of_visits++;
            }
        }
    }

    length = node ? PRetrieve(node, core) : 0;
    return num_of_visits;
}

uint MLCTree::PRetrieve(Node *node, uint *core) const {
    uint length;
#pragma omp parallel
    {
#pragma omp master
        length = Fd(node, core);
    };

    return length;
}

uint MLCTree::Fd(Node *node, uint *core) const {
    auto diff = GetDiff(node);  // collect vertices along the right-most path.
    uint length = diff->num;

    auto chd = GetRmChd(node);

    if (chd) {
#pragma omp task shared (length)
        length += Fd(chd, core + diff->num);
    }

    if (diff->num) {
        memcpy(core, diff->vtx_ptr, diff->num * sizeof(uint));
    }
#pragma omp taskwait
    return length;
}

void MLCTree::Flush(const string &file) {
    auto f = std::ofstream(file);

    f.write(reinterpret_cast<char *> (&dim), sizeof(uint));
    f.write(reinterpret_cast<char *> (&n), sizeof(uint));
    Flush(f);
    f.close();
}

void MLCTree::Flush(std::basic_ofstream<char> &f) {
    unordered_map<Node *, ll_uint> node2id;

    node2id.emplace(root, 0);
    FlushNode(root, node2id, f);
}

void MLCTree::FlushNode(Node *node, unordered_map<Node *, ll_uint> &node2id, std::basic_ofstream<char> &f) const {
    Node *child;
    ll_uint child_id;
    Diff *difference;

    f.write(reinterpret_cast<char *> (node), (long) d_off);

    difference = GetDiff(node);
    f.write(reinterpret_cast<char *>(&difference->num), sizeof(uint));
    if (difference->num) {
        f.write(reinterpret_cast<char *> (difference->vtx_ptr), long(difference->num * sizeof(uint)));
    }

    for (int i = (int) dim - 1; i >= (int) GetIncK(node); i--) {
        child = GetRelChd(node, i);
        child_id = -1;

        if (child) {
            auto iter = node2id.find(child);
            if (iter == node2id.end()) {
                child_id = node2id.size();
                f.write(reinterpret_cast<char *>(&child_id), sizeof(ll_uint));

                node2id.emplace(child, child_id);
                FlushNode(child, node2id, f);
            } else {
                f.write(reinterpret_cast<char *>(&iter->second), sizeof(ll_uint));
            }
        } else f.write(reinterpret_cast<char *>(&child_id), sizeof(ll_uint));
    }
}

MLCTree *MLCTree::Load(const string &file) {
    uint dimension, n_;
    auto f = std::ifstream(file);

    f.read(reinterpret_cast<char *> (&dimension), sizeof(uint));
    f.read(reinterpret_cast<char *> (&n_), sizeof(uint));

    auto *tree = new MLCTree(dimension, n_);
    tree->Load(f);
    f.close();
    return tree;
}

void MLCTree::Load(std::basic_ifstream<char> &f) {
    unordered_map<ll_uint, Node *> id2node;

    root = Allocate(0);
    num_of_nodes++;
    id2node.emplace(0, root);

    LoadNode(root, id2node, f);
}

void MLCTree::LoadNode(Node *node, unordered_map<ll_uint, Node *> &id2node, std::basic_ifstream<char> &f) {
    Node *child;
    ll_uint child_id;
    Diff *difference;

    f.read(reinterpret_cast<char *> (node), (long) d_off);

    difference = GetDiff(node);
    f.read(reinterpret_cast<char *>(&difference->num), sizeof(uint));
    if (difference->num) {
        difference->vtx_ptr = difference_buf.Allocate(difference->num);
        f.read(reinterpret_cast<char *> (difference->vtx_ptr), long(difference->num * sizeof(uint)));
    }

    for (int i = (int) dim - 1; i >= (int) GetIncK(node); i--) {
        f.read(reinterpret_cast<char *>(&child_id), sizeof(ll_uint));

        if (child_id != -1) {
            auto iter = id2node.find(child_id);
            if (iter == id2node.end()) {
                child = Allocate(i);
                num_of_nodes++;

                SetRelChd(node, i, child);
                id2node.emplace(child_id, child);

                LoadNode(child, id2node, f);
            } else SetRelChd(node, i, iter->second);
        }
    }
}

bool MLCTree::Equals(MLCTree &mlc_t) {

    if (dim != mlc_t.dim || n != mlc_t.n || num_of_nodes != mlc_t.num_of_nodes) {
        return false;
    }
    return IsEqual(root, mlc_t.root);
}

bool MLCTree::IsEqual(Node *node, Node *mlc_tree_node) {
    Diff *d, *mlc_d;

    if (node == nullptr && mlc_tree_node == nullptr) return true;
    else if (node == nullptr || mlc_tree_node == nullptr) return false;

    // compare k
    if (!ArrEq(GetK(node), GetK(mlc_tree_node), dim)) {
        return false;
    }

    // compare inc_k
    if ((GetIncK(node)) != (GetIncK(mlc_tree_node))) {
        return false;
    }

    // compare diff
    d = GetDiff(node);
    mlc_d = GetDiff(mlc_tree_node);

    if (d->num != mlc_d->num || (memcmp(d->vtx_ptr, mlc_d->vtx_ptr, d->num * sizeof(uint)) != 0)) {
        return false;
    }

    for (uint i = GetIncK(node); i < dim; i++) {
        if (!IsEqual(GetRelChd(node, i), GetRelChd(mlc_tree_node, i))) {
            return false;
        }
    }
    return true;
}

void MLCTree::MergeBuf(MLCTree *mlc_t) {

    if (mlc_t->num_of_nodes) {
        tree_node_buf.Merge(mlc_t->tree_node_buf);
        difference_buf.Merge(mlc_t->difference_buf);
    } else {
        mlc_t->SetIsolate();
    }
}
