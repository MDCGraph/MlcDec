//
// Created by ldd on 2023/10/28.
//

#include "MLCore_simp.h"


uint MLCore_simp::Extract(MultilayerGraph& mg, const vector<uint> &k_vec_, uint *mlc) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), k, k_vec[ln], length, **degs, **adj_lst;
    CoreIndex kci;

    degs = new uint *[ln];
    for (uint i = 0; i < ln; i++) {
        degs[i] = new uint[n];
    }

    // ======= Initialize =======

    memcpy(k_vec, k_vec_.data(), ln * sizeof(uint));
    for (uint i = 0; i < ln; i++) {
        adj_lst = mg.GetGraph(i).GetAdjLst();
        for (uint j = 0; j < n; j++) {
            degs[i][j] = adj_lst[j][0];
        }
    }

    kci.Init(n);
    kci.Set();

    for (uint i = 0; i < ln; i++) {
        k = k_vec[i];
        if (k) {
            for (uint j = 0; j < n; j++) {
                if (degs[i][j] < k && kci.pos[j] >= kci.e) {
                    kci.Remove(j);
                }
            }
        }
    }

    // ======= Peeling =========
    Peel(mg, kci, k_vec, degs);

    // ======= Collect results=======
    length = kci.n - kci.e;
    memcpy(mlc, kci.vert + kci.e, length * sizeof(uint));


    /* release space */
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;

    return length;
}

void MLCore_simp::Peel(MultilayerGraph &mg, CoreIndex &kci, uint *k_vec, uint **degs) {

    uint ln = mg.GetLayerNumber(), n = mg.GetN(), k, old_e, v, u, v_neighbor, **adj_lst;

    auto &s = kci.s;
    auto &e = kci.e;

    while (s + n < (e << 1)) {  // deletes more
        s = e;
        for (uint i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i).GetAdjLst();
            for (uint j = e; j < n; j++) {
                v = kci.vert[j];
                v_neighbor = 0;
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    if (kci.pos[adj_lst[v][l]] >= s) v_neighbor++;
                }
                if (v_neighbor < k) kci.Remove(v);
                else degs[i][v] = v_neighbor;
            }
        }
    }

    while (s < e) {
        old_e = e;
        for (uint i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i).GetAdjLst();
            for (uint j = s; j < old_e; j++) {
                v = kci.vert[j];
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (kci.pos[u] >= e) {
                        degs[i][u]--;
                        if (degs[i][u] < k) kci.Remove(u);
                    }
                }
            }
        }
        s = old_e;
    }
}