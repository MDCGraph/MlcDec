//
// Created by ldd on 2023/3/10.
//

#ifndef MLCDEC_MLCSTATE_H
#define MLCDEC_MLCSTATE_H

#include "../Graphs/MultilayerGraph.h"
#include "../Structures/IntLinearHeap.h"
#include "CoreIndex.h"
#include "MLCTree.h"
#include "MLCGSmall.h"

struct MLCState {
//
//    MultilayerGraph &mg;
//    G_small *g_small;
//
    Graph **g_layers;
    uint ln{};
    uint ld;
    uint n;

    MLCTree *mlc_t{nullptr};

    uint *relabel{nullptr};
    uint *k_vec;
    uint **degs;

    CoreIndex *kci;

    char *buf{nullptr};

    explicit MLCState(MultilayerGraph &mg) : ln(mg.GetLayerNumber()), n(mg.GetN()), ld(ln - 1) {
        uint n_tln = n * ln;
        size_t size, offset = 0;

        size = ln * sizeof(Graph *) + (ln + (n << 1) + n_tln) * sizeof(uint) + ln * sizeof(uint *) +
               sizeof(CoreIndex);

        buf = new char[size];

        // allocate space for g_layers
        g_layers = reinterpret_cast<Graph **> (buf);
        for (uint i = 0; i < ln; i++) g_layers[i] = &mg.GetGraph(i);
        offset += ln * sizeof(Graph *);

        // allocate space for k_vec
        k_vec = reinterpret_cast<uint *> (buf + offset);
        offset += ln * sizeof(uint);

        // allocate space for degs
        degs = reinterpret_cast<uint **> (buf + offset);
        offset += ln * sizeof(uint *);
        for (uint i = 0; i < ln; i++) {
            degs[i] = reinterpret_cast<uint *>(buf + offset);
            offset += n * sizeof(uint);
        }

        // allocate space for kci
        kci = reinterpret_cast<CoreIndex *>(buf + offset);
        offset += sizeof(CoreIndex);
        kci->Init(n, reinterpret_cast<uint *>(buf + offset));


//
//        offset += (n << 1) * sizeof(uint);
//        print(std::to_string(offset) + " " + std::to_string(size));

//        // allocate space for kvh
//        kvh = reinterpret_cast<IntLinearHeap *>(buf + offset);
//        offset += ln * sizeof(IntLinearHeap);
//        for (uint i = 0; i < ln; i++) {
//            kvh[i].Init(n, reinterpret_cast<uint *>(buf + offset));
//            offset += (3 * n) * sizeof(uint);
//
//            max_d = mg.GetGraph(i).GetMaxDeg();
//            kvh[i].SetBin(max_d, reinterpret_cast<uint *>(buf + offset));
//            offset += (max_d + 1) * sizeof(uint);
//        }

        // allocate space for aux structures
//        aux_cnt = reinterpret_cast<uint *> (buf + offset);
//        offset += ln * sizeof(uint);
//
//        aux_arr = reinterpret_cast<uint **> (buf + offset);
//        offset += ln * sizeof(uint *);
//        for (uint i = 0; i < ln; i++) {
//            aux_arr[i] = reinterpret_cast<uint *>(buf + offset);
//            offset += n * sizeof(uint);
//        }
//
//        aux_sign = reinterpret_cast<bool **> (buf + offset);
//        offset += ln * sizeof(bool *);
//        for (uint i = 0; i < ln; i++) {
//            aux_sign[i] = reinterpret_cast<bool *>(buf + offset);
//            memset(aux_sign[i], false, n * sizeof(bool));
//            offset += n * sizeof(bool);
//        }
    };

    explicit MLCState(MultilayerGraph &mg, uint *buff_, uint small_bnd) : ln(mg.GetLayerNumber()), n(small_bnd),
                                                                          ld(ln - 1) {
        uint off = 0, re_use = 0;
        size_t size, req_size[5 + ln], reuse_size = 0, avai_reuse_size = mg.GetN() * sizeof(uint), new_size = 0;

        req_size[off++] = ln * sizeof(Graph *);                              // graph_layers
        req_size[off++] = ln * sizeof(uint);                                // k_vec
        req_size[off++] = ln * sizeof(uint *);                              // degs
        for (uint i = 0; i < ln; i++) req_size[off++] = n * sizeof(uint);   // degs[i]
        req_size[off++] = sizeof(CoreIndex);                                // Coreindex
        req_size[off++] = (n << 1) * sizeof(uint);                          // Coreindex's vert, pos
//        req_size[off++] = ln * sizeof(uint);                                // aux_cnt;
//        req_size[off++] = ln * sizeof(uint *);                              // aux_arr;
//        for (uint i = 0; i < ln; i++) req_size[off++] = ln * sizeof(uint);  // aux_arr[i]
//        req_size[off++] = ln * sizeof(bool *);                              // aux_sign
//        for (uint i = 0; i < ln; i++) req_size[off++] = n * sizeof(bool);   // aux_sign[i]



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

//        // allocate space for aux structures
//        aux_cnt = reinterpret_cast<uint *> (alloc_next(off++));
//
//        aux_arr = reinterpret_cast<uint **> (alloc_next(off++));
//        for (uint i = 0; i < ln; i++) {
//            aux_arr[i] = reinterpret_cast<uint *>(alloc_next(off++));
//        }
//
//        aux_sign = reinterpret_cast<bool **> (alloc_next(off++));
//        for (uint i = 0; i < ln; i++) {
//            aux_sign[i] = reinterpret_cast<bool *>(alloc_next(off++));
//            memset(aux_sign[i], false, n * sizeof(bool));
//        }

    }

    ~MLCState() {
        delete[] buf;
//
//        if (!mlc_t->Isolate()) {
//            delete mlc_t;
//        }
    }

    // Initialize kci and degs
    void Set() const {
        uint **adj_lst, *deg;

        kci->Set();

        for (uint i = 0; i < ln; i++) {
            adj_lst = g_layers[i]->GetAdjLst();
            deg = degs[i];

            for (uint j = 0; j < n; j++) {
                deg[j] = adj_lst[j][0];
            }
        }
    }

    // Initialize kci and degs with a small ref mlc
    void Set(const uint *ref_mlc, uint ref_size) const {
        uint v, pos;

        kci->s = n - ref_size;
        kci->e = kci->s;

        memset(kci->pos, 0, n * sizeof(uint));

        for (uint j = 0; j < ref_size; j++) {
            v = ref_mlc[j], pos = kci->e + j;
            kci->vert[pos] = v;
            kci->pos[v] = pos;
        }

        uint off = ref_size;
        for (uint l = 0; l < ln; l++) {
            for (uint j = 0; j < ref_size; j++) {
                degs[l][ref_mlc[j]] = ref_mlc[off];
                off++;
            }
        }
    }

    void Set(G_small *ref_g) {
        uint **adj_lst, *deg;

        n = ref_g->GetN();
        relabel = ref_g->GetRelabel();
        kci->Set(n);

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

//    void SetKVHeap() const {
//        uint v;
//
//        for (uint i = 0; i < ln; i++) {
//            kvh[i].Clear();
//            for (uint j = kci->e; j < n; j++) {
//                v = kci->vert[j];
//                kvh[i].Insert(v, degs[i][v]);
//            }
//        }
//    }
//
//    void SetKVHeap(uint inc_k) const {
//        uint v, *deg = degs[inc_k];
//
//        kvh[inc_k].Clear();
//        for (uint j = kci->e; j < n; j++) {
//            v = kci->vert[j];
//            kvh[inc_k].Insert(v, deg[v]);
//        }
//    }

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


    bool BuildBranch(Node *r, uint i) const {
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

        if (i == ln - 1) {
            mlc_t->SetDiff(r, kci->e - old_e, kci->vert + old_e, relabel);
        }

        return kci->e < n;
    }

    bool BuildRMBranch(Node *r) const {
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

        mlc_t->SetDiff(r, kci->e - old_e, kci->vert + old_e, relabel);
        return kci->e < n;
    }

    void BatchBuildRMBranch(Node *r) const {

        bool has_child;
        Node *curr_r = r, *child;

        has_child = BuildRMBranch(curr_r);
//        print(Arr2Str(k_vec, ln));

        while (has_child) {

            child = mlc_t->GetNewTreeNode(k_vec, ld);
            mlc_t->SetRelChd(curr_r, ld, child);

            curr_r = child;

            has_child = BuildRMBranch(curr_r);
//            print(Arr2Str(k_vec, ln));
        }
    }

    void BatchBuildRMBranch(Node *r, uint level) const {

        bool has_child;
        Node *curr_r = r, *child;

        has_child = BuildRMBranch(curr_r);


        while (has_child) {

            auto sum = 0;
            for (uint i = 0; i < ln; i++) {
                sum += k_vec[i];
            }


            child = mlc_t->GetNewTreeNode(k_vec, ld);
            mlc_t->SetRelChd(curr_r, ld, child);
            if (sum >= level) break;

            curr_r = child;

            has_child = BuildRMBranch(curr_r);
            //            print(Arr2Str(k_vec, ln));
        }
    }

    void RecursiveBuild(Node *r, uint inc_k) {

        bool has_child;
        uint old_e;
        Node *child;

        old_e = kci->e;

        for (int i = (int) ld; i >= (int) inc_k; i--) {
            has_child = BuildBranch(r, i);

            if (has_child) {

                child = mlc_t->GetNewTreeNode(k_vec, i);
                mlc_t->SetRelChd(r, i, child);

                RecursiveBuild(child, i);
            }

            k_vec[i]--;
            Restore(old_e);
        }
    }
};


#endif //MLCDEC_MLCSTATE_H

////
//
//#ifdef MY_DEBUG
//if (kci->pos[arr[j]] < kci->e) {
//    cout << "false " << arr[j] << " " << Arr2Str(k_vec, ln) << " " << ld
//    << endl;
//    exit(-1);
//}
//#endif
//

//~MLCState() {
//    delete[] k_vec;
//
//    delete kci;
//    if (kvh) {
//        for (int i = 0; i < ln; i++) {
//            kvh[i].ReleaseBin();
//        }
//    }
//    delete[] kvh;
//
//    if (degs) {
//        for (int i = 0; i < ln; i++) delete[] degs[i];
//        delete[] degs;
//    }
//
//    delete[] aux_arr;
//    delete[] aux_cnt;
//    delete[] aux_sign;
//
//
//    delete[] k_vec;
//}

//
//struct statistics {
//    double init_alloc_time;
//    double new_alloc_time;
//
//    double acc_time;
//    double ini_time_1;
//    double ini_time_2;
//    double peel_time;
//
//    ll_uint path_core;
//    ll_uint rec_core;
//    ll_uint ini_core;
//    ll_uint com_core;
//
//    void init() {
//
//        init_alloc_time = 0;
//        new_alloc_time = 0;
//
//        acc_time = 0;
//        ini_time_1 = 0;
//        ini_time_2 = 0;
//        peel_time = 0;
//
//        path_core = 0;
//        rec_core = 0;
//        ini_core = 0;
//        com_core = 0;
//    }
//};
//
//struct time {
//    int mode;
//    double t;
//
//    time(int m_, double t_) {
//        mode = m_;
//        t = t_;
//    }
//};
