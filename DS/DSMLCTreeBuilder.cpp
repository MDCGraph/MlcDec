//
// Created by ldd on 2023/11/5.
//

#include "DSMLCTreeBuilder.h"

void DSMLCTreeBuilder::Execute(MultilayerGraph & mg, DSMLCTree &mlc_tree) {
    Node *root;
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), k_vec[ln];
    CoreIndex kci;
    uint **degs, **adj_lst;
    
    /* init */
    kci.Init(n);
    kci.Set();

    degs = new uint *[ln];
    for (uint i = 0; i < ln; i++) {
        auto &g = mg.GetGraph(i);
        adj_lst = g.GetAdjLst();

        degs[i] = new uint[n];
        for (uint j = 0; j < n; j++) {
            degs[i][j] = adj_lst[j][0];
        }
    }

    memset(k_vec, 0, ln * sizeof(uint));
    root = mlc_tree.GetNewTreeNode(k_vec, 0);
    mlc_tree.SetRoot(root);

    /* run*/
    BuildSubMLCTree(mg, mlc_tree, kci, degs, root, k_vec, 0);
    
    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
}


void DSMLCTreeBuilder::BuildSubMLCTree(MultilayerGraph& mg, DSMLCTree &mlc_tree, CoreIndex& kci, uint** degs, Node* r, uint* k, uint inc_k) {
    uint * deg, mlc_diff_size, v, old_e, ln = mg.GetLayerNumber(), inc_edge;
    Node *child;
    bool set_inc;

    old_e = kci.e;

    for (auto i = (int) ln - 1; i >= (int) inc_k; i--) {
        deg = degs[i];
        set_inc = false;

        for (uint j = old_e; j < kci.n; j++) {
            v = kci.vert[j];
            if (deg[v] <= k[i]) {
                kci.Remove(v);
            }
        }

        k[i] += 1;
        //        cout << Arr2Str(k, ln) << " " << kci.e << endl;
        if (kci.e != old_e) {
            Peel(mg, kci, degs, k);
            
            mlc_diff_size = kci.e - old_e;
            if (i == ln - 1 && mlc_diff_size) {
                set_inc = true;
                mlc_tree.SetDiff(r, mlc_diff_size, kci.vert + old_e);
            }
        }

        if (kci.e < kci.n) {
            child = mlc_tree.GetNewTreeNode(k, i);
            mlc_tree.SetRelChd(r, i, child);
            BuildSubMLCTree(mg, mlc_tree, kci, degs, child, k, i);

        } // else mlc_tree.SetRelChd(r, i, nullptr);

        k[i]--;

        if (set_inc) RestoreRm(mg, kci, degs, old_e, mlc_tree.GetDiffEC(r));
        else Restore(mg, kci, degs, old_e);
    }
}

void DSMLCTreeBuilder::Peel(MultilayerGraph &mg, CoreIndex &kci, uint **degs, const uint *k) {
    uint old_e;
    uint ln = mg.GetLayerNumber(), **adj_lst, v, u;

    auto &s = kci.s;
    auto &e = kci.e;

    while (s < e) {
        old_e = e;
        for (uint i = 0; i < ln; i++) {
            adj_lst = mg.GetGraph(i).GetAdjLst();

            for (uint j = s; j < old_e; j++) {
                v = kci.vert[j];
                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];
                    degs[i][u]--;
                    if (kci.pos[u] >= e) {
                        if (degs[i][u] < k[i]) {
                            kci.Remove(u);
                        }
                    }
                }
            }
        }
        s = old_e;
    }
}

void DSMLCTreeBuilder::RestoreRm(MultilayerGraph &mg, CoreIndex &kci, uint **degs, uint old_e, uint* diff_nc) {
    uint ln = mg.GetLayerNumber(), v, u, *deg, **adj_lst;
    uint inc;

    for (uint i = 0; i < ln; i++) {
        deg = degs[i];
        adj_lst = mg.GetGraph(i).GetAdjLst();
        inc = 0;

        for (uint j = old_e; j < kci.e; j++) {
            v = kci.vert[j];
            for (uint l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                deg[u]++;

                if (kci.pos[u] >= kci.e || kci.pos[u] >= old_e && u > v) inc++;
            }
        }

        diff_nc[i] = inc;
    }

    kci.e = old_e;
    kci.s = old_e;
}


void DSMLCTreeBuilder::Restore(MultilayerGraph &mg, CoreIndex &kci, uint **degs, uint old_e) {
    uint ln = mg.GetLayerNumber(), v, u, *deg, **adj_lst;

    for (uint i = 0; i < ln; i++) {
        deg = degs[i];
        adj_lst = mg.GetGraph(i).GetAdjLst();

        for (uint j = old_e; j < kci.e; j++) {
            v = kci.vert[j];
            for (uint l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                deg[u]++;
            }
        }
    }

    kci.e = old_e;
    kci.s = old_e;
}