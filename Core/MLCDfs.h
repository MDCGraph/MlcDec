//
// Created by ldd on 2023/6/8.
//

#ifndef MLCDEC_MLCDFS_H
#define MLCDEC_MLCDFS_H


#include "../Structures/IntLinearHeap.h"
#include "../Graphs/MultilayerGraph.h"
#include "CoreIndex.h"

class MLCDfs {
public:
    explicit MLCDfs(MultilayerGraph &mg_) : mg(mg_), ln(mg.GetLayerNumber()), n(mg.GetN()) {

        k_vec = new uint[ln];

        k_core_index.Init(n);
        k_value_heaps = new IntLinearHeap[ln];
        for (uint i = 0; i < ln; i++) k_value_heaps[i].Init(n);

        degs = new uint *[ln];
        for (uint i = 0; i < ln; i++) degs[i] = new uint[n];

        aux_cnt = new uint[ln];
        aux_arr = new uint *[ln];
        for (uint i = 0; i < ln; i++) aux_arr[i] = new uint[n];

        aux_sign = new bool *[ln];
        for (uint i = 0; i < ln; i++) {
            aux_sign[i] = new bool[n];
            memset(aux_sign[i], false, n * sizeof(bool));
        }
    }

    ~MLCDfs() {
        delete[] k_vec;
        delete[] k_value_heaps;

        if (degs) {
            for (uint i = 0; i < ln; i++) delete[] degs[i];
            delete[] degs;
        }

        delete[] aux_cnt;

        if (aux_arr) {
            for (uint i = 0; i < ln; i++) delete[] aux_arr[i];
            delete[] aux_arr;
        }

        if (aux_sign) {
            for (uint i = 0; i < ln; i++) delete[] aux_sign[i];
            delete[] aux_sign;
        }
    }

    void Init() {
        uint **adj_lst, *deg;
        IntLinearHeap *k_value_heap;

        k_core_index.Set();
        for (uint i = 0; i < ln; i++) {
            auto &g = mg.GetGraph(i);

            adj_lst = g.GetAdjLst();
            deg = degs[i];
            k_value_heap = &k_value_heaps[i];

            k_value_heap->SetBin(g.GetMaxDeg());
            for (uint j = 0; j < n; j++) {
                deg[j] = adj_lst[j][0];
                k_value_heap->Insert(j, deg[j]);
            }
        }
    }

    void Peel() {
        uint k, old_s, old_e;
        uint **adj_lst, v, u;
        IntLinearHeap *k_value_heap;

        auto &s = k_core_index.s;
        auto &e = k_core_index.e;

        memset(aux_cnt, 0, ln * sizeof(uint));
        old_s = s;

        while (s < e) {
            old_e = e;
            for (uint i = 0; i < ln; i++) {
                k = k_vec[i];
                adj_lst = mg.GetGraph(i).GetAdjLst();

                for (uint j = s; j < old_e; j++) {
                    v = k_core_index.vert[j];
                    for (uint l = 1; l <= adj_lst[v][0]; l++) {
                        u = adj_lst[v][l];
                        degs[i][u]--;
                        if (k_core_index.pos[u] >= e) {
                            if (degs[i][u] < k) {
                                k_core_index.Remove(u);
                                aux_sign[i][u] = false;
                            } else if (!aux_sign[i][u]) {
                                aux_sign[i][u] = true;
                                aux_arr[i][aux_cnt[i]++] = u;
                            }
                        }
                    }
                }
            }
            s = old_e;
        }

        for (uint i = 0; i < ln; i++) {
            k_value_heap = &k_value_heaps[i];
            for (uint j = 0; j < aux_cnt[i]; j++) {
                u = aux_arr[i][j];
                if (aux_sign[i][u]) {
                    aux_sign[i][u] = false;
                    k_value_heap->Update(u, degs[i][u]);
                }
            }

            for (uint j = old_s; j < s; j++) {
                k_value_heap->Remove(k_core_index.vert[j]);
            }
        }
    }

    void Restore(uint old_e) {
        uint v, u, *deg, **adj_lst;
        IntLinearHeap *k_value_heap;

        memset(aux_cnt, 0, ln * sizeof(uint));
        for (uint i = 0; i < ln; i++) {
            deg = degs[i];
            adj_lst = mg.GetGraph(i).GetAdjLst();

            for (uint j = old_e; j < k_core_index.e; j++) {
                v = k_core_index.vert[j];
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    deg[u]++;
                    if (k_core_index.pos[u] >= k_core_index.e && !aux_sign[i][u]) {
                        aux_sign[i][u] = true;
                        aux_arr[i][aux_cnt[i]++] = u;
                    }
                }
            }

            k_value_heap = &k_value_heaps[i];
            for (uint j = 0; j < aux_cnt[i]; j++) {
                u = aux_arr[i][j];
                k_value_heap->Update(u, deg[u]);
                aux_sign[i][u] = false;
            }

            for (uint j = old_e; j < k_core_index.e; j++) {
                u = k_core_index.vert[j];
                k_value_heap->Insert(u, deg[u]);
            }
        }

        k_core_index.e = old_e;
        k_core_index.s = old_e;
    }

protected:
    MultilayerGraph &mg;
    uint *k_vec;

    CoreIndex k_core_index;
    IntLinearHeap *k_value_heaps;
    uint **degs;
    uint ln;
    uint n;

    uint *aux_cnt;
    uint **aux_arr;
    bool **aux_sign;

};


#endif //MLCDEC_MLCDFS_H
