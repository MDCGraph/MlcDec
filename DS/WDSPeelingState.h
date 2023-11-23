//
// Created by ldd on 2023/10/31.
//

#ifndef MLCDEC_WDSPEELINGSTATE_H
#define MLCDEC_WDSPEELINGSTATE_H

#include "DSMLCTree.h"

struct WDS_state_simp {
    uint ln;
    float *bl;
    const float *w;

    DSMLCTree *mlc_tree;
    Node *max_den_node;
    float max_den;

    uint g_size;
    uint *n_edges;
    float *wd;

    WDS_state_simp(uint ln_, const float *w_, float *bl_, DSMLCTree *mlc_tree_) : ln(ln_), w(w_), bl(bl_),
                                                                                  mlc_tree(mlc_tree_), max_den(0),
                                                                                  g_size(0), max_den_node(nullptr) {
        wd = new float[ln];
        n_edges = new uint[ln];
    }

    ~WDS_state_simp() {
        delete[] wd;
        delete[] bl;
    }

    void Clear() {
        memset(n_edges, 0, ln * sizeof(uint));
        g_size = 0;
    }

    float GetDensity(uint delta_g_size) {
        float den, max_ww = 0;

        g_size += delta_g_size;
        for (uint l = 0; l < ln; l++) {
            wd[l] = (float) n_edges[l] / (float) g_size * w[l];
        }

        std::sort(wd, wd + ln);
        for (int i = (int) ln - 1; i >= 0; i--) {
            den = wd[i] * bl[i];

            if (den > max_ww) {
                max_ww = den;
            }
        }

        return max_ww;
    }

    void RecSearch(Node *r) {
        float max_ww;
        Node *child;
        DS_Diff *diff;

        /* find the first child representing different core with r*/
        child = r;
        while (child) {
            diff = mlc_tree->GetDiff(child);
            child = mlc_tree->GetRmChd(child);
            if (diff->num) break;
        }

        if (child) RecSearch(child);

        for (uint i = 0; i < ln; i++) {
            n_edges[i] += diff->n_edge[i];
        }
        max_ww = GetDensity(diff->num);

        if (max_ww > max_den) {
            max_den = max_ww;;
            max_den_node = r;
        }
    }
};

struct WDS_state {
//    MultilayerGraph *mg{nullptr};

    uint *buf;
    uint *n_edges;

    bool *sign;
    bool *counted;

    float *bl;
    float *wd;

    uint n{0};
    uint ln{0};

    const float *w;

    WDS_state(MultilayerGraph *mg, const float *w_, float beta) : n(mg->GetN()), ln(mg->GetLayerNumber()),
                                                                  w(w_) {
        size_t size = (n + ln) * sizeof(uint) + (n << 1) * sizeof(bool) + (ln << 1) * sizeof(float);

        buf = reinterpret_cast<uint *> (new char[size]);
        n_edges = buf + n;
        sign = reinterpret_cast<bool *>(n_edges + ln);
        counted = sign + n;

        bl = reinterpret_cast<float *>(counted + n);
        wd = bl + ln;

        for (uint i = 0; i < ln; i++) {
            bl[i] = (float) pow(ln - i, beta);
        }

    }

    ~WDS_state() {
        delete[] buf;
    }

    void Clear() const {
        memset(n_edges, 0, ln * sizeof(uint));
        memset(sign, false, n * sizeof(bool));
        memset(counted, false, n * sizeof(bool));
    }

    void SetSign(const uint *arr, uint length) const {
        for (uint i = 0; i < length; i++) {
            sign[arr[i]] = true;
        }
    }

    void SetCounted(const uint *arr, uint length) const {
        for (uint i = 0; i < length; i++) {
            counted[arr[i]] = true;
        }
    }

    void IncNEdges(uint **adj_lst, uint l, const uint *arr, uint length) const {
        uint v, u, inc = 0;

        for (uint i = 0; i < length; i++) {
            v = arr[i];
            for (uint j = 1; j <= adj_lst[v][0]; j++) {
                u = adj_lst[v][j];
                if (counted[u] || (sign[u] && u > v)) {
                    inc++;
                }
            }
        }
        n_edges[l] += inc;
    }

    float ComputeDensityByInc(MultilayerGraph *mg, const uint *arr, uint length, uint g_size) const {
        float den, max_ww = 0;

        for (uint l = 0; l < ln; l++) {
            IncNEdges(mg->GetGraph(l).GetAdjLst(), l, arr, length);
            wd[l] = (float) n_edges[l] / (float) g_size * w[l];
        }

        std::sort(wd, wd + ln);
        for (int i = (int) ln - 1; i >= 0; i--) {
            den = wd[i] * bl[i];

            if (den > max_ww) {
                max_ww = den;
            }
        }

        return max_ww;
    }

    float ComputeDensityByInc(Graph **graphs, const uint *arr, uint length, uint g_size) const {
        float den, max_ww = 0;

        for (uint l = 0; l < ln; l++) {
            IncNEdges(graphs[l]->GetAdjLst(), l, arr, length);
            wd[l] = (float) n_edges[l] / (float) g_size * w[l];
        }

        std::sort(wd, wd + ln);
        for (int i = (int) ln - 1; i >= 0; i--) {
            den = wd[i] * bl[i];

            if (den > max_ww) {
                max_ww = den;
            }
        }

        return max_ww;
    }
};

struct WDS_search_state : WDS_state {

    MultilayerGraph *mg{};
    Node *max_den_node{nullptr};
    float max_den{0};

    MLCTree *mlc_t{nullptr};

    WDS_search_state(MultilayerGraph *mg_, const float *w_, float beta_, MLCTree *mlc_t_) : mg(mg_),
                                                                                            WDS_state(mg_, w_, beta_),
                                                                                            mlc_t(mlc_t_) {}

    uint SearchRMPath(Node *r, uint offset = 0) {
        uint size, *buf_start;
        float max_ww;
        Node *child;
        Diff *diff;

        /* find the first child representing different core with r*/
        child = r;
        while (child) {
            diff = mlc_t->GetDiff(child);
            child = mlc_t->GetRmChd(child);
            if (diff->num) break;
        }

        buf_start = buf + offset;
        memcpy(buf_start, diff->vtx_ptr, diff->num * sizeof(uint));
        size = diff->num;

        if (child) {
            size += SearchRMPath(child, offset + diff->num);
        }

        SetSign(buf_start, diff->num);
        max_ww = ComputeDensityByInc(mg, buf_start, diff->num, size);
        SetCounted(buf_start, diff->num);

        if (max_ww > max_den) {
            max_den = max_ww;;
            max_den_node = r;
        }

        return size;
    }

};

struct WDS_peeling_state : WDS_state {

    MDataBuf<uint> *k_buf;
    Graph **g_layers;
    uint ld;
//    uint n_nodes = 0;
    uint *relabel;

    uint *k_vec;
    uint **degs;

    CoreIndex *kci;

    uint *max_den_k{};
    uint *core{};
    uint length{0};
    uint *max_den_relabel;
    float max_den{0};

    char *buf{nullptr};

    WDS_peeling_state(MultilayerGraph *mg_, MDataBuf<uint> *k_buf_, const float *w_, float beta_, float curr_max_den,
                      uint *buff_, uint small_bnd) : WDS_state(mg_, w_, beta_), k_buf(k_buf_), ld(ln - 1),
                                                     max_den(curr_max_den) {

        n = small_bnd;
        uint off = 0, re_use = 0;
        size_t size, req_size[7 + ln], reuse_size = 0, avai_reuse_size = mg_->GetN() * sizeof(uint), new_size = 0;

        req_size[off++] = ln * sizeof(Graph *);                              // graph_layers
        req_size[off++] = ln * sizeof(uint);                                // k_vec
        req_size[off++] = ln * sizeof(uint *);                              // degs
        for (uint i = 0; i < ln; i++) req_size[off++] = n * sizeof(uint);   // degs[i]
        req_size[off++] = sizeof(CoreIndex);                                // Coreindex
        req_size[off++] = (n << 1) * sizeof(uint);                          // Coreindex's vert, pos
        req_size[off++] = ln * sizeof(uint);
        req_size[off++] = n * sizeof(uint);

        while (re_use < off) {
            reuse_size += req_size[re_use];
            if (reuse_size > avai_reuse_size) {
                break;
            } else re_use++;
        }

        if (re_use < off) {
            for (uint i = re_use; i < off; i++) {
                new_size += req_size[i];
            }
        }

        if (new_size) buf = new char[new_size];


        size = 0;
        auto avai_buf = reinterpret_cast<char *>(buff_);

        auto alloc_next = [&](uint id) {
            if (id == re_use) {
                size = 0;
                avai_buf = buf;
            }

            auto ptr = avai_buf + size;
            size += req_size[id];
            return ptr;
        };

        off = 0;

        // allocate space for g_layers
        g_layers = reinterpret_cast<Graph **> (alloc_next(off++));

        // allocate space for k_vec
        k_vec = reinterpret_cast<uint *> (alloc_next(off++));

        // allocate space for degs
        degs = reinterpret_cast<uint **> (alloc_next(off++));
        for (uint i = 0; i < ln; i++) {
            degs[i] = reinterpret_cast<uint *>(alloc_next(off++));
        }

        // allocate space for kci
        kci = reinterpret_cast<CoreIndex *>(alloc_next(off++));
        kci->Init(n, reinterpret_cast<uint *>(alloc_next(off++)));

        // allocate space for max_den_k
        max_den_k = reinterpret_cast<uint *> (alloc_next(off++));
        core = reinterpret_cast<uint *> (alloc_next(off++));
    }

    void Set(G_small *ref_g) {
        uint **adj_lst, *deg;

        n = ref_g->GetN();
        kci->Set(n);
        relabel = ref_g->GetRelabel();

        for (uint i = 0; i < ln; i++) {

            g_layers[i] = &(ref_g->GetGraph(i));
            adj_lst = g_layers[i]->GetAdjLst();
            deg = degs[i];

            for (uint j = 0; j < n; j++) {
                deg[j] = adj_lst[j][0];
            }
        }
    }

    inline void SetKVec(uint *k) const {
        memcpy(k_vec, k, ln * sizeof(uint));
    }

    /* no wds-related operations */
    [[nodiscard]] uint PeelCore() const {
        uint k, old_e, v, u, v_neighbor, **adj_lst;

        auto &s = kci->s;
        auto &e = kci->e;

        for (uint i = 0; i < ln; i++) {
            k = k_vec[i];
            if (k) {
                for (uint j = e; j < n; j++) {
                    v = kci->vert[j];
                    if (kci->pos[v] >= e && degs[i][v] < k) {
                        kci->Remove(v);
                    }
                }
            }
        }

        // ======= Peeling =========
        while (s + n < (e << 1)) {  // deletes more
            s = e;
            for (uint i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = g_layers[i]->GetAdjLst();
                for (uint j = e; j < n; j++) {
                    v = kci->vert[j];
                    v_neighbor = 0;
                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        if (kci->pos[adj_lst[v][l]] >= s) v_neighbor++;
                    }
                    if (v_neighbor < k) kci->Remove(v);
                    else degs[i][v] = v_neighbor;
                }
            }
        }

        while (s < e) {
            old_e = e;
            for (uint i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = g_layers[i]->GetAdjLst();
                for (uint j = s; j < old_e; j++) {
                    v = kci->vert[j];
                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        u = adj_lst[v][l];
                        if (kci->pos[u] >= e) {
                            degs[i][u]--;
                            if (degs[i][u] < k) kci->Remove(u);
                        }
                    }
                }
            }
            s = old_e;
        }

        return n - e;
    }

    uint *GetDupK() {
        auto new_k = k_buf->Allocate(ln);
        memcpy(new_k, k_vec, ln * sizeof(uint));
        return new_k;
    }

    void Peel() const {
        uint k, old_e;
        uint **adj_lst, v, u;

        auto &s = kci->s;
        auto &e = kci->e;

        while (s < e) {
            old_e = e;
            for (uint i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = g_layers[i]->GetAdjLst();

                for (uint j = s; j < old_e; j++) {
                    v = kci->vert[j];
                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        u = adj_lst[v][l];
                        degs[i][u]--;
                        if (kci->pos[u] >= e && degs[i][u] < k) {
                            kci->Remove(u);
                        }
                    }
                }
            }
            s = old_e;
        }
    }

    void PeelLast() const {
        uint i, k, old_e, v_neighbor;
        uint **adj_lst, v, u;

        auto &s = kci->s;
        auto &e = kci->e;

        while (s + n < (e << 1)) {  // deletes more
            s = e;
            for (i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = g_layers[i]->GetAdjLst();

                for (uint j = e; j < n; j++) {
                    v = kci->vert[j];
                    v_neighbor = 0;

                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        if (kci->pos[adj_lst[v][l]] >= s) v_neighbor++;
                    }

                    degs[i][v] = v_neighbor;
                    if (v_neighbor < k) {
                        kci->Remove(v);
                    }
                }
            }
        }

        while (s < e) { // removed more
            old_e = e;
            for (i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = g_layers[i]->GetAdjLst();

                for (uint j = s; j < old_e; j++) {
                    v = kci->vert[j];
                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        u = adj_lst[v][l];
                        if (kci->pos[u] >= e) {
                            degs[i][u]--;
                            if (degs[i][u] < k) {
                                kci->Remove(u);
                            }
                        }
                    }
                }
            }

            s = old_e;
        }
    }

    void Restore(uint old_e) const {
        uint v, u, *deg, **adj_lst;

        for (uint i = 0; i < ln; i++) {
            deg = degs[i];
            adj_lst = g_layers[i]->GetAdjLst();
            for (uint j = old_e; j < kci->e; j++) {
                v = kci->vert[j];
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    deg[adj_lst[v][l]]++;
                }
            }
        }

        kci->e = old_e;
        kci->s = old_e;
    }

    [[nodiscard]] bool BuildBranch(uint i) const {
        uint v, old_e;

        k_vec[i] += 1;
        old_e = kci->e;

        for (uint j = old_e; j < n; j++) {
            v = kci->vert[j];
            if (degs[i][v] < k_vec[i]) {
                kci->Remove(v);
            }
        }

        if (old_e != kci->e) Peel();

        return kci->e < n;
    }

    [[nodiscard]] uint BuildRMBranch() const {
        uint v, old_e;

        k_vec[ld] += 1;
        old_e = kci->e;

        for (uint j = old_e; j < n; j++) {
            v = kci->vert[j];
            if (degs[ld][v] < k_vec[ld]) {
                kci->Remove(v);
            }
        }

        if (old_e != kci->e) PeelLast();

        return kci->e;
    }

    void SearchRMBranch() {
        uint e = kci->e, size, *buf_start, child_e;
        float max_ww;

        buf_start = kci->vert + kci->s;
        child_e = BuildRMBranch();

        fflush(stdout);

        if (child_e < n) {
//            n_nodes++;
            SearchRMBranch();
        }

        k_vec[ld]--;

        size = child_e - e;

        if (size) {
            SetSign(buf_start, size);
            max_ww = ComputeDensityByInc(g_layers, buf_start, size, n - e);
            SetCounted(buf_start, size);

            if (max_ww > max_den) {
                max_den = max_ww;
                length = n - e;
                memcpy(max_den_k, k_vec, ln * sizeof(uint));
                memcpy(core, buf_start, length * sizeof(uint));
                max_den_relabel = relabel;
            }
        }
    }

    uint RecursiveSearchChild(uint inc_k) {

        bool has_child;
        uint old_e = kci->e, size, *buf_start, child_e = n;
        float max_ww;

        for (int i = (int) ld; i >= (int) inc_k; i--) {
            has_child = BuildBranch(i);

            if (has_child) {
//                n_nodes++;

                if (i != ld) Clear();
                child_e = RecursiveSearchChild(i);
            }

            k_vec[i]--;

            if (i == ld) {

                size = child_e - old_e;

                buf_start = kci->vert + old_e;
                SetSign(buf_start, size);
                max_ww = ComputeDensityByInc(g_layers, buf_start, size, n - old_e);
                SetCounted(buf_start, size);

                if (max_ww > max_den) {
                    max_den = max_ww;
                    length = n - old_e;
                    memcpy(max_den_k, k_vec, ln * sizeof(uint));
                    memcpy(core, buf_start, length * sizeof(uint));
                    max_den_relabel = relabel;
                }
            }

            Restore(old_e);

        }

        return old_e;
    }
};


#endif //MLCDEC_WDSPEELINGSTATE_H
