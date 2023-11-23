//
// Created by ldd on 2023/3/6.
//

#include "MLCore.h"


MLCore::MLCore(MultilayerGraph &mg_) : mg(mg_), ln(mg.GetLayerNumber()), n(mg.GetN()) {
    k_vec = new uint[ln];
    degs = new uint *[ln];
    for (uint i = 0; i < ln; i++) {
        degs[i] = new uint[n];
    }
}

MLCore::~MLCore() {
    delete[] k_vec;
    if (degs) {
        for (uint i = 0; i < ln; i++) {
            delete[] degs[i];
        }
        delete[] degs;
    }
}

uint MLCore::Extract(const vector<uint> &k_vec_, uint *mlc, PEELING_MODE peeling_mode) {
    uint k, length, **adj_lst;
    CoreIndex core_index;

    // ======= Initialize =======

    memcpy(k_vec, k_vec_.data(), ln * sizeof(uint));
    for (uint i = 0; i < ln; i++) {
        adj_lst = mg.GetGraph(i).GetAdjLst();
        for (uint j = 0; j < n; j++) {
            degs[i][j] = adj_lst[j][0];
        }
    }

    core_index.Init(n);
    core_index.Set();

    for (uint i = 0; i < ln; i++) {
        k = k_vec[i];
        if (k) {
            for (uint j = 0; j < n; j++) {
                if (degs[i][j] < k && core_index.pos[j] >= core_index.e) {
                    core_index.Remove(j);
                }
            }
        }
    }

    // ======= Peeling =========
    if (peeling_mode == BY_LAYER) PeelByLayer(core_index);
    else PeelByVtx(core_index);

    // ======= Collect results=======
    length = core_index.n - core_index.e;
    memcpy(mlc, core_index.vert + core_index.e, length * sizeof(uint));

    return length;
}

void MLCore::PeelByLayer(CoreIndex & ci) {
    uint k, old_e, v, u, v_neighbor, **adj_lst;

    auto &s = ci.s;
    auto &e = ci.e;

    while (s + n < (e << 1)) {  // deletes more
        s = e;
        for (uint i = 0; i < ln; i++) {
            k = k_vec[i];
            adj_lst = mg.GetGraph(i).GetAdjLst();
            for (uint j = e; j < n; j++) {
                v = ci.vert[j];
                v_neighbor = 0;
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    if (ci.pos[adj_lst[v][l]] >= s) v_neighbor++;
                }
                if (v_neighbor < k) ci.Remove(v);
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
                v = ci.vert[j];
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (ci.pos[u] >= e) {
                        degs[i][u]--;
                        if (degs[i][u] < k) ci.Remove(u);
                    }
                }
            }
        }
        s = old_e;
    }
}

void MLCore::PeelByVtx(CoreIndex & ci) {
    uint old_e, v, u, v_neighbor, **adj_lst[ln];

    for (uint i = 0; i < ln; i++) {
        adj_lst[i] = mg.GetGraph(i).GetAdjLst();
    }

    auto &s = ci.s;
    auto &e = ci.e;

    while (s + n < (e << 1)) {  // deletes more
        s = e;
        for (uint j = e; j < n; j++) {
            v = ci.vert[j];
            for (uint i = 0; i < ln; i++) {
                v_neighbor = 0;
                for (uint l = 1; l <= adj_lst[i][v][0]; l++) {
                    if (ci.pos[adj_lst[i][v][l]] >= s) v_neighbor++;
                }
                if (v_neighbor < k_vec[i]) {
                    ci.Remove(v);
                    break;
                } else degs[i][v] = v_neighbor;
            }
        }
    }

    while (s < e) {
        old_e = e;
        for (uint j = s; j < old_e; j++) {
            v = ci.vert[j];
            for (uint i = 0; i < ln; i++) {
                for (uint l = 1; l <= adj_lst[i][v][0]; l++) {
                    u = adj_lst[i][v][l];
                    if (ci.pos[u] >= e) {
                        degs[i][u]--;
                        if (degs[i][u] < k_vec[i]) ci.Remove(u);
                    }
                }
            }
        }
        s = old_e;
    }
}