//
// Created by ldd on 2023/3/29.
//

#include "CoreCube.h"

CoreCube::CoreCube(MultilayerGraph &mg_) : mg(mg_), ln(mg.GetLayerNumber()), n(mg.GetN()) {

    uint max_deg, deg;

    layer_idt = new bool[ln];
    selected_layer = new uint[ln];

    sups = new uint *[ln];
    for (uint i = 0; i < ln; i++) {
        sups[i] = new uint[n];
    }

    q1 = new uint[n];
    q2 = new uint[n];

    max_deg = mg.GetGraph(0).GetMaxDeg();
    for (uint i = 1; i < ln; i++) {
        deg = mg.GetGraph(i).GetMaxDeg();
        if (deg > max_deg) {
            max_deg = deg;
        }
    }

    counter = new uint[max_deg + 1];

    in_q1 = new bool[n];
    in_q2 = new bool[n];

    memset(in_q1, false, n * sizeof(bool));
    memset(in_q2, false, n * sizeof(bool));
}

CoreCube::~CoreCube() {
    delete[] layer_idt;
    delete[] selected_layer;

    if (sups) {
        for (uint i = 0; i < ln; i++) {
            delete[] sups[i];
        }
        delete[] sups;
    }

    delete[] q1;
    delete[] q2;
    delete[] counter;

    delete[] in_q1;
    delete[] in_q2;
}

void CoreCube::TD(Cube &cube_) {
    uint v, *coreness, *c_coreness, **adj_lst, l;

    cube = &cube_;
    cube->Set(n, GetNCombs());

    for (z = 1; z <= ln; z++) {

        InitLInd();
        do {
            SetSLayer();

            coreness = cube->GetCoreness(GetCurrCoreOffset());
            l = selected_layer[0];

            if (z == 1) {
                adj_lst = mg.GetGraph(l).GetAdjLst();

                for (v = 0; v < n; v++) {
                    coreness[v] = adj_lst[v][0];
                }

            } else {

                memcpy(coreness, cube->GetCoreness(GetAnsCoreOffset(l)), n * sizeof(uint));

                for (l = 1; l < z; l++) {
                    c_coreness = cube->GetCoreness(GetAnsCoreOffset(selected_layer[l]));
                    for (v = 0; v < n; v++) {
                        if (c_coreness[v] < coreness[v]) {
                            coreness[v] = c_coreness[v];
                        }
                    }
                }
            }

            CoreTD();

        } while (std::prev_permutation(layer_idt, layer_idt + ln));
    }
}

void CoreCube::InitLInd() {
    for (uint i = 0; i < z; i++) layer_idt[i] = true;
    for (uint i = z; i < ln; i++) layer_idt[i] = false;
}

void CoreCube::SetSLayer() {
    uint offset = 0;
    for (uint i = 0; i < ln; i++) {
        if (layer_idt[i]) selected_layer[offset++] = i;
    }
}

uint CoreCube::GetCurrCoreOffset() {
    uint sum = 0;
    for (uint i = 0; i < z; i++) {
        sum += (uint) pow(2, selected_layer[i]);
    }
    return sum - 1;
}

uint CoreCube::GetAnsCoreOffset(uint mask) const {
    uint sum = 0;
    for (uint i = 0; i < z; i++) {
        if (selected_layer[i] != mask) {
            sum += (uint) pow(2, selected_layer[i]);
        }
    }
    return sum - 1;
}

void CoreCube::CoreTD() {
    uint v, u, sl, c0, mv, n_nbr, *sup, **adj_lst, *coreness;
    uint idx1 = 0, idx2;

    coreness = cube->GetCoreness(GetCurrCoreOffset());

    // Initialize sups(v, l), line 1
    for (uint i = 0; i < z; i++) {
        sl = selected_layer[i];

        adj_lst = mg.GetGraph(sl).GetAdjLst();
        sup = sups[sl];

        for (v = 0; v < n; v++) {
            c0 = coreness[v];

            n_nbr = 0;
            for (uint j = 1; j <= adj_lst[v][0]; j++) {
                if (coreness[adj_lst[v][j]] >= c0) n_nbr++;
            }

            sup[v] = n_nbr;
            if (n_nbr < c0 && !in_q1[v]) {
                q1[idx1++] = v;
                in_q1[v] = true;
            }
        }
    }

    while (idx1) {
        idx2 = 0;

        for (uint i = 0; i < idx1; i++) {
            v = q1[i];
            in_q1[v] = false;

            c0 = coreness[v];

            for (uint j = 0; j < z; j++) {

                sl = selected_layer[j];
                adj_lst = mg.GetGraph(sl).GetAdjLst();

                // compute M and update C
                memset(counter, 0, (coreness[v] + 1) * sizeof(uint));
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (coreness[u] >= coreness[v]) counter[coreness[v]] += 1;
                    else counter[coreness[u]] += 1;
                }

                mv = 0;
                for (uint l = coreness[v]; l >= 1; l--) {
                    if (counter[l] >= l) {
                        mv = l;
                        break;
                    }
                    counter[l - 1] += counter[l];
                }

                if (coreness[v] > mv) coreness[v] = mv;
            }


            for (uint j = 0; j < z; j++) {
                sl = selected_layer[j];
                adj_lst = mg.GetGraph(sl).GetAdjLst();
                sup = sups[sl];

                // update sup[v]
                n_nbr = 0;
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    if (coreness[u] >= coreness[v]) n_nbr++;

                    if (coreness[u] <= c0 && coreness[u] > coreness[v]) {
                        sup[u]--;

                        if (sup[u] < coreness[u] && !in_q1[u] && !in_q2[u]) {
                            q2[idx2++] = u;
                            in_q2[u] = true;
                        }
                    }
                }

                sup[v] = n_nbr;
                if (sup[v] < coreness[v] && !in_q2[v]) {
                    q2[idx2++] = v;
                    in_q2[v] = true;
                }
            }
        }

        idx1 = idx2;
        std::swap(q1, q2);
        std::swap(in_q1, in_q2);
    }

}

void CoreCube::GenHybridStorage(HYB_Cube &hyb_cube) {
    uint n_pos, n_diff, min_n_diff, min_diff_ref, l, idx, ref_id;
    uint *coreness;

    hyb_cube.Set(n, GetNCombs());

    for (z = 1; z <= ln; z++) {

        InitLInd();
        do {
            SetSLayer();

            l = GetCurrCoreOffset();
            coreness = cube->GetCoreness(l);

            n_pos = GetNPos(coreness, n);

            if (z == 1) {
                StoreAbs(coreness, n_pos, hyb_cube.GetHybCoreness(l));

            } else {

                ref_id = GetAnsCoreOffset(selected_layer[0]);
                min_n_diff = GetNDiff(coreness, cube->GetCoreness(ref_id), n);
                min_diff_ref = ref_id;

                for (idx = 1; idx < z; idx++) {

                    ref_id = GetAnsCoreOffset(selected_layer[idx]);
                    n_diff = GetNDiff(coreness, cube->GetCoreness(ref_id), n);

                    if (n_diff < min_n_diff) {
                        min_n_diff = n_diff;
                        min_diff_ref = ref_id;
                    }
                }

                if (n_pos <= min_n_diff) {
                    StoreAbs(coreness, n_pos, hyb_cube.GetHybCoreness(l));

                } else {
                    StoreRef(coreness, min_diff_ref, min_n_diff, hyb_cube.GetHybCoreness(l));
                }
            }
        } while (std::prev_permutation(layer_idt, layer_idt + ln));
    }

}

void CoreCube::StoreAbs(uint *coreness, uint n_pos, HYB_Cube1d &hyb_cube1d) const {
    uint v, idx;

    hyb_cube1d.set(-1, n_pos);

    idx = 0;
    for (v = 0; v < n; v++) {
        if (coreness[v]) hyb_cube1d.vtx_cn[idx++] = {v, coreness[v]};
    }
}

void CoreCube::StoreRef(uint *coreness, uint ref_id, uint n_diff, HYB_Cube1d &hyb_cube1d) const {
    uint v, idx, *ref_coreness;

    hyb_cube1d.set(ref_id, n_diff);
    ref_coreness = cube->GetCoreness(ref_id);

    idx = 0;
    for (v = 0; v < n; v++) {
        if (coreness[v] != ref_coreness[v]) {
            hyb_cube1d.vtx_cn[idx++] = {v, coreness[v] - ref_coreness[v]};
        }
    }
}
