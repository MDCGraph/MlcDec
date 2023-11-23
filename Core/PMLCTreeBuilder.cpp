//
// Created by ldd on 2023/10/23.
//

#include "PMLCTreeBuilder.h"

/***** use *****/

void PMLCTreeBuilder::InterPathPRun(MultilayerGraph &mg, MLCTree &mlc_tree) {

    uint ln = mg.GetLayerNumber(), n = mg.GetN(), n_threads = omp_get_max_threads();
    MLCState **mlc_state;

    // Set root node
    Node *root = mlc_tree.GetNewTreeNode();
    mlc_tree.SetRoot(root);

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];

#pragma omp parallel
    {
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp master
        {

            for (uint i = 0; i < ln - 1; i++) {
#pragma omp task
                {
                    BuildRMPath(mlc_state, root, i);
                }
            }

            mlc_s->Set();
            mlc_s->SetKVec(MLCTree::GetK(root));
            mlc_s->BatchBuildRMBranch(root);

        };
    };

    // =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }
}

void PMLCTreeBuilder::InterPathPRun_with_merge(MultilayerGraph &mg, MLCTree &mlc_tree, float merge_factor) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), n_threads = omp_get_max_threads(), merge_bound = uint(
            merge_factor * (float) n);
    MLCState **mlc_state;

    // Set root node
    Node *root = mlc_tree.GetNewTreeNode();
    mlc_tree.SetRoot(root);

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];

#pragma omp parallel
    {
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp master
        {

            for (uint i = 0; i < ln - 1; i++) {
#pragma omp task
                {
                    BuildRMPath_m(mlc_state, root, i, merge_bound);
                }
            }

            mlc_s->Set();
            mlc_s->SetKVec(MLCTree::GetK(root));
            mlc_s->BatchBuildRMBranch(root);

        };
    };

    // =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }
}

void PMLCTreeBuilder::PRun_async(MultilayerGraph &mg, MLCTree &mlc_tree, uint level) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1, n_threads = omp_get_max_threads();
    uint k[ln], *re_queue, **degs, ***adj_lsts;
    uint re_top = 0, iv_off = 0, n_fin = 0, re, le;
    uint k_sum = 1, max_sr_size = 0;

    MLCState **mlc_state;

    DataBuf<uint> sg_ngh_buf;
    DataBuf<uint *> sg_adj_buf;

    PCoreIndex kci;
    Node *r, *child, **node_queue;
    Diff *diff;

    vector<SUB_ROOT> sub_roots;

    bool not_fin = false, gen_gs;

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];
    sub_roots.reserve(ln);

    bool fin = false;

    /* init structs */
    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, node_queue);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    r = mlc_tree.GetNewTreeNode(k, 0);
    mlc_tree.SetRoot(r);

    node_queue[re_top] = r;
    re_queue[re_top++] = 0;
    k[top]++;

    /* parallel run */
#pragma omp parallel
    {
        uint old_e, v, u, s, e, e_be, iv_off_pt, e1;
        uint *deg, *buf = new uint[n];

        /* init */
#pragma omp for schedule(static) nowait
        for (uint i = 0; i < n; i++) {
            kci.inva[i] = false;
            kci.upd[i] = true;
            kci.vert[i] = i;
        }

#pragma omp for schedule(static) collapse(2)
        for (uint i = 0; i < ln; i++) {
            for (uint j = 0; j < n; j++) {
                degs[i][j] = adj_lsts[i][j][0];
            }
        }

        /* recursive compute */
        while (top >= 0) {

            /* compute ml-core */
            e = 0;
            deg = degs[top];
            old_e = n_fin;

#pragma omp master
            iv_off = 0;

            //#pragma omp master
            //            print_k(k, ln, k_sum, kci.e);

            /* compute first batch of removal */
#pragma omp for schedule(static) nowait
            for (uint j = old_e; j < kci.n; j++) {
                v = kci.vert[j];
                if (deg[v] < k[top] && !__atomic_test_and_set(&kci.inva[v], __ATOMIC_SEQ_CST)) {
                    buf[e++] = v;
                }
            }

            e_be = n_fin;
            /* peel */
            PeelPud(mg, kci, degs, k, buf, e);
            s = __sync_fetch_and_add(&kci.e, e);

#pragma omp barrier
#pragma omp master
            {

                /* prepare augmentation */
                if (top == ln - 1) {

                    diff = mlc_tree.GetDiff(r);
                    diff->num = kci.e - n_fin;
                    diff->vtx_ptr = mlc_tree.AllocateDiff(diff->num);

                } else {
                    diff = nullptr;
                }

                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    gen_gs = false;

                    /* continue exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    if (top != ln - 1) top = ln - 1;

                    node_queue[re_top] = child;
                    re_queue[re_top++] = kci.e;

                    k[top]++;
                    r = child;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else if (kci.e < kci.n && k_sum == level) {

                    gen_gs = true;

                    /* stop exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    uint new_n = n - kci.e;
                    uint *m_addr = sg_ngh_buf.Allocate(ln);
                    uint *relable = sg_ngh_buf.Allocate(new_n);
                    sub_roots.emplace_back(SUB_ROOT{child, new_n, m_addr, relable});

                    if (new_n > max_sr_size) max_sr_size = new_n;

                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }


                } else {

                    gen_gs = false;

                    /* backtrace*/
                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }
                }
            }


#pragma omp barrier

            auto new_s = s - e_be;
            if (diff) {
                for (uint i = 0; i < e; i++) {
                    diff->vtx_ptr[new_s + i] = buf[i];
                }
            }

            if (le < n) {

                /* reorganize kci.vert */
                e1 = n;
#pragma omp for schedule(static)
                for (uint j = e_be; j < kci.n; j++) {
                    v = kci.vert[j];
                    if (!kci.inva[v]) {
                        buf[--e1] = v;
                    }
                }

                for (uint i = 0; i < e; i++) {
                    kci.vert[s + i] = buf[i];
                }

                iv_off_pt = le + __sync_fetch_and_add(&iv_off, n - e1) - e1;

                for (uint i = e1; i < n; i++) {
                    kci.vert[iv_off_pt + i] = buf[i];
                }

#pragma omp barrier
            }


            if (gen_gs) {
                auto &sub_r = sub_roots.back();

#pragma omp master
                {
                    memcpy(sub_r.label, kci.vert + le, sub_r.n * sizeof(uint));
                };

#pragma omp for nowait
                for (uint i = 0; i < ln; i++) {
                    uint sum = 0;
                    for (uint j = le; j < n; j++) {
                        sum += degs[i][kci.vert[j]];
                    }
                    sub_r.m[i] = sum + sub_r.n;
                }
            }

            if (not_fin) {
                break;
            }

            if (re == INT32_MAX) {
                for (uint i = 0; i < e; i++) {
                    kci.upd[buf[i]] = false;
                }

            } else {

                /* restore */
                auto to_e = le;
                for (int kk = (int) re - 1; kk >= 0; kk--) {
                    auto ree = re_queue[re_top - 1 + kk];

#pragma omp for schedule(static) // nowait
                    for (uint j = ree; j < to_e; j++) {
                        v = kci.vert[j];
                        kci.inva[v] = false;
                        kci.upd[v] = true;
                    }

#pragma omp for schedule(static) collapse(2)
                    for (uint i = 0; i < ln; i++) {
                        for (uint j = ree; j < to_e; j++) {
                            v = kci.vert[j];
                            for (uint l = 1; l <= adj_lsts[i][v][0]; l++) {
                                u = adj_lsts[i][v][l];
                                if (kci.upd[u]) __sync_fetch_and_add(&degs[i][u], 1);
                            }
                        }
                    }

                    to_e = ree;
                }
            }
#pragma omp barrier
        }


#pragma omp barrier


#pragma omp for//nowait
        for (auto &sub_root : sub_roots) {

            auto sub_r = &sub_root;
            auto small_g = ComputeSmallG(mg, sub_r, sg_adj_buf, sg_ngh_buf, buf);

            sub_r->m = reinterpret_cast<uint *>(small_g);
        }


//#pragma omp barrier


        /* ========== prepare ============*/
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg, buf, max_sr_size);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp barrier
#pragma omp master
        {

#pragma omp taskloop
            for (auto &sub_r: sub_roots) {

                auto g_small = reinterpret_cast<G_small *>(sub_r.m);

                // initialize new tasks
                uint inc_k = mlc_tree.GetIncK(sub_r.r);
                for (uint i = inc_k; i < ln - 1; i++) {
#pragma omp task
                    {
                        BuildRMPath(mlc_state, g_small, sub_r.r, i);
                    }
                }

                auto mlc_ss = mlc_state[omp_get_thread_num() << 6];

                mlc_ss->Set(g_small);
                mlc_ss->SetKVec(MLCTree::GetK(sub_r.r));
                mlc_ss->BatchBuildRMBranch(sub_r.r);
            }
        };

#pragma omp barrier
        delete[] buf;
    }

    //     =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }

    for (auto &sub_r:sub_roots) {
        delete reinterpret_cast<G_small *>(sub_r.m);
    }


    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;
    delete[] node_queue;

    sg_ngh_buf.Release();
    sg_adj_buf.Release();

}

void PMLCTreeBuilder::PRun_async_with_merge(MultilayerGraph &mg, MLCTree &mlc_tree, uint level, float merge_factor) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1, n_threads = omp_get_max_threads();
    uint k[ln], *re_queue, **degs, ***adj_lsts;
    uint re_top = 0, iv_off = 0, n_fin = 0, re, le;
    uint k_sum = 1, max_sr_size = 0;

    MLCState **mlc_state;

    DataBuf<uint> sg_ngh_buf;
    DataBuf<uint *> sg_adj_buf;

    PCoreIndex kci;
    Node *r, *child, **node_queue;
    Diff *diff;

    vector<SUB_ROOT> sub_roots;

    bool not_fin = false, gen_gs;

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];
    sub_roots.reserve(ln);

    bool fin = false;

    /* init structs */
    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, node_queue);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    r = mlc_tree.GetNewTreeNode(k, 0);
    mlc_tree.SetRoot(r);

    node_queue[re_top] = r;
    re_queue[re_top++] = 0;
    k[top]++;

    /* parallel run */
#pragma omp parallel
    {
        uint old_e, v, u, s, e, e_be, iv_off_pt, e1;
        uint *deg, *buf = new uint[n];

        /* init */
#pragma omp for schedule(static) nowait
        for (uint i = 0; i < n; i++) {
            kci.inva[i] = false;
            kci.upd[i] = true;
            kci.vert[i] = i;
        }

#pragma omp for schedule(static) collapse(2)
        for (uint i = 0; i < ln; i++) {
            for (uint j = 0; j < n; j++) {
                degs[i][j] = adj_lsts[i][j][0];
            }
        }

        /* recursive compute */
        while (top >= 0) {

            /* compute ml-core */
            e = 0;
            deg = degs[top];
            old_e = n_fin;

#pragma omp master
            iv_off = 0;

            //#pragma omp master
            //            print_k(k, ln, k_sum, kci.e);

            /* compute first batch of removal */
#pragma omp for schedule(static) nowait
            for (uint j = old_e; j < kci.n; j++) {
                v = kci.vert[j];
                if (deg[v] < k[top] && !__atomic_test_and_set(&kci.inva[v], __ATOMIC_SEQ_CST)) {
                    buf[e++] = v;
                }
            }

            e_be = n_fin;
            /* peel */
            PeelPud(mg, kci, degs, k, buf, e);
            s = __sync_fetch_and_add(&kci.e, e);

#pragma omp barrier
#pragma omp master
            {

                /* prepare augmentation */
                if (top == ln - 1) {

                    diff = mlc_tree.GetDiff(r);
                    diff->num = kci.e - n_fin;
                    diff->vtx_ptr = mlc_tree.AllocateDiff(diff->num);

                } else {
                    diff = nullptr;
                }

                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    gen_gs = false;

                    /* continue exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    if (top != ln - 1) top = ln - 1;

                    node_queue[re_top] = child;
                    re_queue[re_top++] = kci.e;

                    k[top]++;
                    r = child;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else if (kci.e < kci.n && k_sum == level) {

                    gen_gs = true;

                    /* stop exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    uint new_n = n - kci.e;
                    uint *m_addr = sg_ngh_buf.Allocate(ln);
                    uint *relable = sg_ngh_buf.Allocate(new_n);
                    sub_roots.emplace_back(SUB_ROOT{child, new_n, m_addr, relable});

                    if (new_n > max_sr_size) max_sr_size = new_n;

                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }


                } else {

                    gen_gs = false;

                    /* backtrace*/
                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }
                }
            }


#pragma omp barrier

            auto new_s = s - e_be;
            if (diff) {
                for (uint i = 0; i < e; i++) {
                    diff->vtx_ptr[new_s + i] = buf[i];
                }
            }

            if (le < n) {

                /* reorganize kci.vert */
                e1 = n;
#pragma omp for schedule(static)
                for (uint j = e_be; j < kci.n; j++) {
                    v = kci.vert[j];
                    if (!kci.inva[v]) {
                        buf[--e1] = v;
                    }
                }

                for (uint i = 0; i < e; i++) {
                    kci.vert[s + i] = buf[i];
                }

                iv_off_pt = le + __sync_fetch_and_add(&iv_off, n - e1) - e1;

                for (uint i = e1; i < n; i++) {
                    kci.vert[iv_off_pt + i] = buf[i];
                }

#pragma omp barrier
            }


            if (gen_gs) {
                auto &sub_r = sub_roots.back();

#pragma omp master
                {
                    memcpy(sub_r.label, kci.vert + le, sub_r.n * sizeof(uint));
                };

#pragma omp for nowait
                for (uint i = 0; i < ln; i++) {
                    uint sum = 0;
                    for (uint j = le; j < n; j++) {
                        sum += degs[i][kci.vert[j]];
                    }
                    sub_r.m[i] = sum + sub_r.n;
                }
            }

            if (not_fin) {
                break;
            }

            if (re == INT32_MAX) {
                for (uint i = 0; i < e; i++) {
                    kci.upd[buf[i]] = false;
                }

            } else {

                /* restore */
                auto to_e = le;
                for (int kk = (int) re - 1; kk >= 0; kk--) {
                    auto ree = re_queue[re_top - 1 + kk];

#pragma omp for schedule(static) // nowait
                    for (uint j = ree; j < to_e; j++) {
                        v = kci.vert[j];
                        kci.inva[v] = false;
                        kci.upd[v] = true;
                    }

#pragma omp for schedule(static) collapse(2)
                    for (uint i = 0; i < ln; i++) {
                        for (uint j = ree; j < to_e; j++) {
                            v = kci.vert[j];
                            for (uint l = 1; l <= adj_lsts[i][v][0]; l++) {
                                u = adj_lsts[i][v][l];
                                if (kci.upd[u]) __sync_fetch_and_add(&degs[i][u], 1);
                            }
                        }
                    }

                    to_e = ree;
                }
            }
#pragma omp barrier
        }


#pragma omp barrier


#pragma omp for//nowait
        for (auto &sub_root : sub_roots) {

            auto sub_r = &sub_root;
            auto small_g = ComputeSmallG(mg, sub_r, sg_adj_buf, sg_ngh_buf, buf);

            sub_r->m = reinterpret_cast<uint *>(small_g);
        }


#pragma omp barrier


        /* ========== prepare ============*/
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg, buf, max_sr_size);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp barrier
#pragma omp master
        {

#pragma omp taskloop
            for (auto &sub_r: sub_roots) {

                auto g_small = reinterpret_cast<G_small *>(sub_r.m);

                // initialize new tasks
                uint inc_k = mlc_tree.GetIncK(sub_r.r);
                for (uint i = inc_k; i < ln - 1; i++) {
#pragma omp task
                    {
                        BuildRMPath_m(mlc_state, g_small, sub_r.r, i, merge_factor);
                    }
                }

                auto mlc_ss = mlc_state[omp_get_thread_num() << 6];

                mlc_ss->Set(g_small);
                mlc_ss->SetKVec(MLCTree::GetK(sub_r.r));
                mlc_ss->BatchBuildRMBranch(sub_r.r);
            }
        };

#pragma omp barrier
        delete[] buf;
    }


    //     =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }

    for (auto &sub_r:sub_roots) {
        delete reinterpret_cast<G_small *>(sub_r.m);
    }


    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;
    delete[] node_queue;

    sg_ngh_buf.Release();
    sg_adj_buf.Release();

}

void PMLCTreeBuilder::InitDeg(MultilayerGraph &mg, uint **&degs, uint ***&adj_lsts) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN();

    degs = new uint *[ln];
    adj_lsts = new uint **[ln];
    for (uint i = 0; i < ln; i++) {
        degs[i] = new uint[n];
        adj_lsts[i] = mg.GetGraph(i).GetAdjLst();
    }
}

void PMLCTreeBuilder::InitQueue(MultilayerGraph &mg, uint *&re_queue, Node **&node_queue) {
    uint ln = mg.GetLayerNumber(), deg_sum = ln;
    for (uint i = 0; i < ln; i++) {
        deg_sum += mg.GetGraph(i).GetMaxDeg();
    }

    re_queue = new uint[deg_sum];
    node_queue = new Node *[deg_sum];
}

void PMLCTreeBuilder::BuildRMPath(MLCState **mlc_state, Node *r, uint inc_k) {
    uint ld;

    Node *child;

    MLCState *mlc_s = mlc_state[omp_get_thread_num() << 6];
    MLCTree *mlc_t = mlc_s->mlc_t;

    mlc_s->Set();
    mlc_s->SetKVec(MLCTree::GetK(r));
    mlc_s->k_vec[inc_k] += 1;

    ld = mlc_s->ld;

    if (mlc_s->PeelCore()) {

        child = mlc_t->GetNewTreeNode(mlc_s->k_vec, inc_k);
        mlc_t->SetRelChd(r, inc_k, child);

        mlc_s->BatchBuildRMBranch(child);

        for (auto i = inc_k; i < ld; i++) {
#pragma omp task
            {
                BuildRMPath(mlc_state, child, i);
            }
        }
    }
}

void PMLCTreeBuilder::BuildRMPath(MLCState **mlc_state, G_small *ref_g, Node *r, uint inc_k) {
    uint ld;

    Node *child;

    MLCState *mlc_s = mlc_state[omp_get_thread_num() << 6];
    MLCTree *mlc_t = mlc_s->mlc_t;

    mlc_s->Set(ref_g);
    mlc_s->SetKVec(MLCTree::GetK(r));
    mlc_s->k_vec[inc_k] += 1;

    ld = mlc_s->ld;

    if (mlc_s->PeelCore()) {

        child = mlc_t->GetNewTreeNode(mlc_s->k_vec, inc_k);
        mlc_t->SetRelChd(r, inc_k, child);

        mlc_s->BatchBuildRMBranch(child);

        for (auto i = inc_k; i < ld; i++) {
#pragma omp task
            {
                BuildRMPath(mlc_state, ref_g, child, i);
            }
        }
    }
}

void PMLCTreeBuilder::BuildRMPath_m(MLCState **mlc_state, Node *r, uint inc_k, uint merge_bound) {
    uint ld, sg_size;
    Node *child;

    MLCState *mlc_s = mlc_state[omp_get_thread_num() << 6];
    MLCTree *mlc_t = mlc_s->mlc_t;

    mlc_s->Set();
    mlc_s->SetKVec(MLCTree::GetK(r));
    mlc_s->k_vec[inc_k] += 1;

    ld = mlc_s->ld;

    sg_size = mlc_s->PeelCore();

    if (sg_size) {

        child = mlc_t->GetNewTreeNode(mlc_s->k_vec, inc_k);
        mlc_t->SetRelChd(r, inc_k, child);

        if (sg_size <= merge_bound) {
            mlc_s->RecursiveBuild(child, inc_k);

        } else {

            mlc_s->BatchBuildRMBranch(child);

            for (auto i = inc_k; i < ld; i++) {
#pragma omp task
                BuildRMPath_m(mlc_state, child, i, merge_bound);
            }
        }
    }
}

void PMLCTreeBuilder::BuildRMPath_m(MLCState **mlc_state, G_small *ref_g, Node *r, uint inc_k, float merge_factor) {
    uint ld, sg_size;
    Node *child;

    MLCState *mlc_s = mlc_state[omp_get_thread_num() << 6];
    MLCTree *mlc_t = mlc_s->mlc_t;

    mlc_s->Set(ref_g);
    mlc_s->SetKVec(MLCTree::GetK(r));
    mlc_s->k_vec[inc_k] += 1;

    ld = mlc_s->ld;

    sg_size = mlc_s->PeelCore();

    if (sg_size) {

        child = mlc_t->GetNewTreeNode(mlc_s->k_vec, inc_k);
        mlc_t->SetRelChd(r, inc_k, child);

        if (sg_size <= (uint) ((float) ref_g->GetN() * merge_factor)) {
            mlc_s->RecursiveBuild(child, inc_k);

        } else {

            mlc_s->BatchBuildRMBranch(child);

            for (auto i = inc_k; i < ld; i++) {
#pragma omp task
                BuildRMPath_m(mlc_state, ref_g, child, i, merge_factor);
            }
        }
    }
}

G_small *
PMLCTreeBuilder::ComputeSmallG(MultilayerGraph &mg, SUB_ROOT *sub_r, DataBuf<uint *> &adj_buf, DataBuf<uint> &nbr_buf,
                               uint *buf) {


    uint all_n = mg.GetN(), ln = mg.GetLayerNumber(), sum = 0, off = 0;
    auto g_small = new G_small(sub_r->n, mg.GetLayerNumber());
    uint **new_adj_lst;
    uint *new_ngh_buf;

    memset(buf, -1, all_n * sizeof(uint));
    for (uint i = 0; i < sub_r->n; i++) {
        buf[sub_r->label[i]] = i;
    }

    for (uint i = 0; i < ln; i++) {
        sum += sub_r->m[i];
    }

#pragma omp critical
    {
        new_adj_lst = adj_buf.Allocate(ln * sub_r->n);

    }

#pragma omp critical
    {

        new_ngh_buf = nbr_buf.Allocate(sum);
    };

    for (uint i = 0; i < ln; i++) {
        new_adj_lst[i] = new_ngh_buf + off;
        off += sub_r->m[i];
    }


    g_small->BuildSmallG(mg, sub_r->label, buf, new_adj_lst, new_ngh_buf);

    return g_small;
}


void PMLCTreeBuilder::MergeMLCTree(MLCTree &mlc_tree, MLCState **mlc_state, uint n) {
    ll_uint n_nodes_tot = mlc_tree.GetNumOfNodes();

    for (uint i = 1; i < n; i++) {
        auto mlc_t = mlc_state[i << 6]->mlc_t;
        n_nodes_tot += mlc_t->GetNumOfNodes();
        mlc_tree.MergeBuf(mlc_t);
    }
    mlc_tree.SetNumOfNodes(n_nodes_tot);
}

void PMLCTreeBuilder::PeelPud(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf,
                              uint &e) {

    uint ln = mg.GetLayerNumber(), old_e, **adj_lst, v, u, du;
    uint s = 0;
    int t = omp_get_thread_num();

    while (s < e) {
        old_e = e;

        for (uint layer = t; layer < ln + t; layer++) {

            auto i = layer % ln;
            adj_lst = mg.GetGraph(i).GetAdjLst();

            for (uint j = s; j < old_e; j++) {
                v = buf[j];

                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];

                    if (kci.upd[u]) {
                        du = __sync_fetch_and_sub(&degs[i][u], 1);

                        if (du == k[i]) {
                            if (!__atomic_test_and_set(&kci.inva[u], __ATOMIC_SEQ_CST)) {
                                buf[e++] = u;
                            }
                        }
                    }

                }
            }
        }

        s = old_e;
    }
}


/***** not in use *****/
void PMLCTreeBuilder::InnerCorePRun(MultilayerGraph &mg, MLCTree &mlc_tree, uint level) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1;
    uint k[ln], *re_queue, *nbr, *nbr_bk, *sorted_nbr, *count, **degs, ***adj_lsts;
    uint re_top = 0, iv_off = 0, n_fin = 0, re, nbr_off = 0, compact_nbr_off = 0;
    uint num_buckets, num_blocks, block_size, le, k_sum = 1;

    PCoreIndex kci;
    Node *r, *child, **node_queue;
    Diff *diff;

    bool not_fin = false, nbr_sort;

    /* init structs */

    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, node_queue);
    InitNbrRes(mg, nbr, sorted_nbr, nbr_bk, count);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    r = mlc_tree.GetNewTreeNode(k, 0);
    mlc_tree.SetRoot(r);

    node_queue[re_top] = r;
    re_queue[re_top++] = 0;
    k[top]++;

    /* parallel run */
#pragma omp parallel
    {
        uint old_e, v, s, e, e_be, iv_off_pt, start, e1, local_nbr_off;
        uint *deg, *buf = new uint[n];

        /* init */
#pragma omp for schedule(static) nowait
        for (uint i = 0; i < n; i++) {
            kci.inva[i] = false;
            kci.upd[i] = true;
            kci.vert[i] = i;
        }

#pragma omp for schedule(static) collapse(2)
        for (uint i = 0; i < ln; i++) {
            for (uint j = 0; j < n; j++) {
                degs[i][j] = adj_lsts[i][j][0];
            }
        }

        /* recursive compute */
        while (top >= 0) {

            /* compute ml-core */
            e = 0;
            deg = degs[top];
            old_e = n_fin;

#pragma omp master
            iv_off = 0;
//
//                        #pragma omp master
//                        print_k(k, ln, k_sum, kci.e);

            /* compute first batch of removal */
#pragma omp for schedule(static) nowait
            for (uint j = old_e; j < kci.n; j++) {
                v = kci.vert[j];

                if (deg[v] < k[top]) {
                    buf[e++] = v;
                    kci.inva[v] = true;
                }
            }

            /* peel */
            start = 0;
            e_be = n_fin;
            while (true) {
                old_e = e;

                for (uint layer = 0; layer < ln; layer++) {

                    Graph &g = mg.GetGraph(layer);
                    deg = degs[layer];

                    GenOneHopeNgh(g, kci.upd, buf, start, old_e, nbr, sorted_nbr, nbr_off, compact_nbr_off);

#pragma omp barrier
                    local_nbr_off = nbr_off;
#pragma omp master
                    {
                        if (!layer) not_fin = false;
                        nbr_sort = (nbr_off > 4096);
                        if (nbr_sort) {
                            num_blocks = get_n_blk(nbr_off);
                            num_buckets = get_n_bkt(nbr_off);
                            block_size = ((nbr_off - 1) / num_blocks) + 1;

                            SegSort(nbr, nbr_bk, sorted_nbr, nbr_off, count, num_buckets, num_blocks, block_size);
                        }

                        compact_nbr_off = 0;
                    };
#pragma omp barrier

                    if (!nbr_sort) {
#pragma omp for schedule(static) nowait
                        for (size_t l = 0; l < local_nbr_off; l++) {
                            v = nbr[l];

                            auto dv = __sync_fetch_and_sub(&deg[v], 1);
                            if (dv == k[layer]) {
                                if (!__atomic_test_and_set(&kci.inva[v], __ATOMIC_SEQ_CST)) {
                                    buf[e++] = v;
                                }
                            }
                        }
                    } else if (local_nbr_off) {

#pragma omp for schedule(static) nowait
                        for (uint i = 0; i < num_buckets; i++) {
                            for (uint j = 0; j < num_blocks; j++) {
                                uint bs = (uint) (j * block_size);
                                uint snb_off = count[j * num_buckets + i];
                                uint ct;
                                if (i == num_buckets - 1) {
                                    ct = (std::min((uint) (bs + block_size), local_nbr_off) - bs) - snb_off;
                                } else {
                                    ct = count[j * num_buckets + i + 1] - snb_off;
                                }

                                snb_off += bs;

                                for (uint l = 0; l < ct; l++) {

                                    v = sorted_nbr[snb_off + l];
                                    if (deg[v] == k[layer] && !kci.inva[v]) {
                                        buf[e++] = v;
                                        kci.inva[v] = true;
                                    }
                                    deg[v]--;
                                }
                            }
                        };
                    }

                    if (old_e != e) __atomic_test_and_set(&not_fin, std::memory_order_seq_cst);

#pragma omp master
                    nbr_off = 0;

#pragma omp barrier
                }

                if (!not_fin) {
                    break;
                }
                start = old_e;
            }

            s = __sync_fetch_and_add(&kci.e, e);

#pragma omp barrier
#pragma omp master
            {

                /* prepare augmentation */
                if (top == ln - 1) {

                    diff = mlc_tree.GetDiff(r);
                    diff->num = kci.e - n_fin;
                    diff->vtx_ptr = mlc_tree.AllocateDiff(diff->num);

                } else {
                    diff = nullptr;
                }

                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    /* continue exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    if (top != ln - 1) top = ln - 1;

                    node_queue[re_top] = child;
                    re_queue[re_top++] = kci.e;

                    k[top]++;
                    r = child;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else {

                    /* backtrace*/
                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }
                }
            }


#pragma omp barrier

            auto new_s = s - e_be;
            if (diff) {
                for (uint i = 0; i < e; i++) {
                    diff->vtx_ptr[new_s + i] = buf[i];
                }
            }

            if (not_fin) {
                break;
            }

            if (le < n) {

                /* reorganize kci.vert */
                e1 = n;
#pragma omp for schedule(static)
                for (uint j = e_be; j < kci.n; j++) {
                    v = kci.vert[j];
                    if (!kci.inva[v]) {
                        buf[--e1] = v;
                    }
                }

                for (uint i = 0; i < e; i++) {
                    kci.vert[s + i] = buf[i];
                }

                iv_off_pt = le + __sync_fetch_and_add(&iv_off, n - e1) - e1;

                for (uint i = e1; i < n; i++) {
                    kci.vert[iv_off_pt + i] = buf[i];
                }

#pragma omp barrier
            }


            if (re == INT32_MAX) {
                for (uint i = 0; i < e; i++) {
                    kci.upd[buf[i]] = false;
                }

            } else {

                /* restore */
                auto to_e = le;
                for (int kk = (int) re - 1; kk >= 0; kk--) {
                    auto ree = re_queue[re_top - 1 + kk];

                    e = 0;

#pragma omp for schedule(static)// nowait
                    for (uint j = ree; j < to_e; j++) {
                        v = kci.vert[j];
                        kci.inva[v] = false;
                        kci.upd[v] = true;
                        buf[e++] = v;
                    }

                    for (uint layer = 0; layer < ln; layer++) {

                        Graph &g = mg.GetGraph(layer);
                        deg = degs[layer];
                        GenOneHopeNgh(g, kci.upd, buf, 0, e, nbr, sorted_nbr, nbr_off, compact_nbr_off);

#pragma omp barrier
                        local_nbr_off = nbr_off;
#pragma omp master
                        {
                            nbr_sort = (nbr_off > 4096);
                            if (nbr_sort) {
                                num_blocks = get_n_blk(nbr_off);
                                num_buckets = get_n_bkt(nbr_off);
                                block_size = ((nbr_off - 1) / num_blocks) + 1;

                                SegSort(nbr, nbr_bk, sorted_nbr, nbr_off, count, num_buckets, num_blocks, block_size);
                            }

                            compact_nbr_off = 0;
                        };

#pragma omp barrier
                        if (local_nbr_off) {
                            if (!nbr_sort) {
#pragma omp for schedule(static) nowait
                                for (size_t l = 0; l < local_nbr_off; l++) {
                                    v = nbr[l];
                                    __sync_fetch_and_add(&deg[v], 1);
                                }

                            } else {

#pragma omp for schedule(static) nowait
                                for (uint i = 0; i < num_buckets; i++) {
                                    for (uint j = 0; j < num_blocks; j++) {
                                        uint bs = (uint) (j * block_size);
                                        uint snb_off = count[j * num_buckets + i];
                                        uint ct;
                                        if (i == num_buckets - 1) {
                                            ct = (std::min((uint) (bs + block_size), local_nbr_off) - bs) - snb_off;
                                        } else {
                                            ct = count[j * num_buckets + i + 1] - snb_off;
                                        }

                                        snb_off += bs;

                                        for (uint l = 0; l < ct; l++) {
                                            v = sorted_nbr[snb_off + l];
                                            deg[v]++;
                                        }
                                    }
                                };
                            }

                        }

#pragma omp master
                        nbr_off = 0;
#pragma omp barrier

                    }

                    to_e = ree;
                }
            }

#pragma omp barrier
        }

        delete[] buf;
    }


    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;
    delete[] node_queue;

    delete[] nbr_bk;
    delete[] count;
}

void PMLCTreeBuilder::PRun(MultilayerGraph &mg, MLCTree &mlc_tree, uint level) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1, n_threads = omp_get_max_threads();
    uint k[ln], *re_queue, *nbr, *nbr_bk, *sorted_nbr, *count, **degs, ***adj_lsts;
    uint re_top = 0, iv_off = 0, n_fin = 0, re, nbr_off = 0, compact_nbr_off = 0;
    uint num_buckets, num_blocks, block_size, le;
    uint k_sum = 1, max_sr_size = 0;

    MLCState **mlc_state;
    MDataBuf<uint> nm;
    DataBuf<uint> sg_ngh_buf;
    DataBuf<uint *> sg_adj_buf;

    PCoreIndex kci;
    Node *r, *child, **node_queue;
    Diff *diff;

    vector<SUB_ROOT> sub_roots;

    bool not_fin = false, nbr_sort, gen_gs;

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];

    /* init structs */
    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, node_queue);
    InitNbrRes(mg, nbr, sorted_nbr, nbr_bk, count);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    r = mlc_tree.GetNewTreeNode(k, 0);
    mlc_tree.SetRoot(r);

    node_queue[re_top] = r;
    re_queue[re_top++] = 0;
    k[top]++;

    /* parallel run */
#pragma omp parallel
    {
        uint old_e, v, s, e, e_be, iv_off_pt, start, e1, local_nbr_off;
        uint *deg, *buf = new uint[n];

        /* init */
#pragma omp for schedule(static) nowait
        for (uint i = 0; i < n; i++) {
            kci.inva[i] = false;
            kci.upd[i] = true;
            kci.vert[i] = i;
        }

#pragma omp for schedule(static) collapse(2)
        for (uint i = 0; i < ln; i++) {
            for (uint j = 0; j < n; j++) {
                degs[i][j] = adj_lsts[i][j][0];
            }
        }

        /* recursive compute */
        while (top >= 0) {

            /* compute ml-core */
            e = 0;
            deg = degs[top];
            old_e = n_fin;

#pragma omp master
            iv_off = 0;

//                        #pragma omp master
//                        print_k(k, ln, k_sum, kci.e);

            /* compute first batch of removal */
#pragma omp for schedule(static) nowait
            for (uint j = old_e; j < kci.n; j++) {
                v = kci.vert[j];

                if (deg[v] < k[top]) {
                    buf[e++] = v;
                    kci.inva[v] = true;
                }
            }

            /* peel */
            start = 0;
            e_be = n_fin;
            while (true) {
                old_e = e;

                for (uint layer = 0; layer < ln; layer++) {

                    Graph &g = mg.GetGraph(layer);
                    deg = degs[layer];

                    GenOneHopeNgh(g, kci.upd, buf, start, old_e, nbr, sorted_nbr, nbr_off, compact_nbr_off);

#pragma omp barrier
                    local_nbr_off = nbr_off;
#pragma omp master
                    {
                        if (!layer) not_fin = false;
                        nbr_sort = (nbr_off > 4096);
                        if (nbr_sort) {
                            num_blocks = get_n_blk(nbr_off);
                            num_buckets = get_n_bkt(nbr_off);
                            block_size = ((nbr_off - 1) / num_blocks) + 1;

                            SegSort(nbr, nbr_bk, sorted_nbr, nbr_off, count, num_buckets, num_blocks, block_size);
                        }

                        compact_nbr_off = 0;
                    };
#pragma omp barrier

                    if (!nbr_sort) {
#pragma omp for schedule(static) nowait
                        for (size_t l = 0; l < local_nbr_off; l++) {
                            v = nbr[l];

                            auto dv = __sync_fetch_and_sub(&deg[v], 1);
                            if (dv == k[layer]) {
                                if (!__atomic_test_and_set(&kci.inva[v], __ATOMIC_SEQ_CST)) {
                                    buf[e++] = v;
                                }
                            }
                        }
                    } else if (local_nbr_off) {

#pragma omp for schedule(static) nowait
                        for (uint i = 0; i < num_buckets; i++) {
                            for (uint j = 0; j < num_blocks; j++) {
                                uint bs = (uint) (j * block_size);
                                uint snb_off = count[j * num_buckets + i];
                                uint ct;
                                if (i == num_buckets - 1) {
                                    ct = (std::min((uint) (bs + block_size), local_nbr_off) - bs) - snb_off;
                                } else {
                                    ct = count[j * num_buckets + i + 1] - snb_off;
                                }

                                snb_off += bs;

                                for (uint l = 0; l < ct; l++) {

                                    v = sorted_nbr[snb_off + l];
                                    if (deg[v] == k[layer] && !kci.inva[v]) {
                                        buf[e++] = v;
                                        kci.inva[v] = true;
                                    }
                                    deg[v]--;
                                }
                            }
                        };
                    }

                    if (old_e != e) __atomic_test_and_set(&not_fin, std::memory_order_seq_cst);

#pragma omp master
                    nbr_off = 0;

#pragma omp barrier
                }

                if (!not_fin) {
                    break;
                }
                start = old_e;
            }

            s = __sync_fetch_and_add(&kci.e, e);

#pragma omp barrier
#pragma omp master
            {

                /* prepare augmentation */
                if (top == ln - 1) {

                    diff = mlc_tree.GetDiff(r);
                    diff->num = kci.e - n_fin;
                    diff->vtx_ptr = mlc_tree.AllocateDiff(diff->num);

                } else {
                    diff = nullptr;
                }

                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    gen_gs = false;

                    /* continue exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    if (top != ln - 1) top = ln - 1;

                    node_queue[re_top] = child;
                    re_queue[re_top++] = kci.e;

                    k[top]++;
                    r = child;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else if (kci.e < kci.n && k_sum == level) {

                    gen_gs = true;

                    /* stop exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    uint new_n = n - kci.e;
                    uint *m_addr = nm.Allocate(ln);
                    uint *relable = sg_ngh_buf.Allocate(new_n);
                    sub_roots.emplace_back(SUB_ROOT{child, new_n, m_addr, relable});

                    if (new_n > max_sr_size) max_sr_size = new_n;

                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }


                } else {

                    gen_gs = false;

                    /* backtrace*/
                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];
                        r = node_queue[re_top];

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;

                        re = k[top]; // store previous k
                        k[top] = 0;
                        k[--top]++;

                        re_top++;
                    }
                }
            }


#pragma omp barrier

            auto new_s = s - e_be;
            if (diff) {
                for (uint i = 0; i < e; i++) {
                    diff->vtx_ptr[new_s + i] = buf[i];
                }
            }

            if (le < n) {

                /* reorganize kci.vert */
                e1 = n;
#pragma omp for schedule(static)
                for (uint j = e_be; j < kci.n; j++) {
                    v = kci.vert[j];
                    if (!kci.inva[v]) {
                        buf[--e1] = v;
                    }
                }

                for (uint i = 0; i < e; i++) {
                    kci.vert[s + i] = buf[i];
                }

                iv_off_pt = le + __sync_fetch_and_add(&iv_off, n - e1) - e1;

                for (uint i = e1; i < n; i++) {
                    kci.vert[iv_off_pt + i] = buf[i];
                }

#pragma omp barrier
            }


            if (gen_gs) {
                auto &sub_r = sub_roots.back();

#pragma omp master
                {
                    memcpy(sub_r.label, kci.vert + le, sub_r.n * sizeof(uint));
                    assert(sub_r.n == n - le);
                };

#pragma omp for nowait
                for (uint i = 0; i < ln; i++) {
                    uint sum = 0;
                    for (uint j = le; j < n; j++) {
                        sum += degs[i][kci.vert[j]];
                    }
                    sub_r.m[i] = sum + sub_r.n;
                }
            }

            if (not_fin) {
                break;
            }

            if (re == INT32_MAX) {
                for (uint i = 0; i < e; i++) {
                    kci.upd[buf[i]] = false;
                }

            } else {

                /* restore */
                auto to_e = le;
                for (int kk = (int) re - 1; kk >= 0; kk--) {
                    auto ree = re_queue[re_top - 1 + kk];

                    e = 0;

#pragma omp for schedule(static)// nowait
                    for (uint j = ree; j < to_e; j++) {
                        v = kci.vert[j];
                        kci.inva[v] = false;
                        kci.upd[v] = true;
                        buf[e++] = v;
                    }

                    for (uint layer = 0; layer < ln; layer++) {

                        Graph &g = mg.GetGraph(layer);
                        deg = degs[layer];
                        GenOneHopeNgh(g, kci.upd, buf, 0, e, nbr, sorted_nbr, nbr_off, compact_nbr_off);

#pragma omp barrier
                        local_nbr_off = nbr_off;
#pragma omp master
                        {
                            nbr_sort = (nbr_off > 4096);
                            if (nbr_sort) {
                                num_blocks = get_n_blk(nbr_off);
                                num_buckets = get_n_bkt(nbr_off);
                                block_size = ((nbr_off - 1) / num_blocks) + 1;

                                SegSort(nbr, nbr_bk, sorted_nbr, nbr_off, count, num_buckets, num_blocks, block_size);
                            }

                            compact_nbr_off = 0;
                        };

#pragma omp barrier
                        if (local_nbr_off) {
                            if (!nbr_sort) {
#pragma omp for schedule(static) nowait
                                for (size_t l = 0; l < local_nbr_off; l++) {
                                    v = nbr[l];
                                    __sync_fetch_and_add(&deg[v], 1);
                                }

                            } else {

#pragma omp for schedule(static) nowait
                                for (uint i = 0; i < num_buckets; i++) {
                                    for (uint j = 0; j < num_blocks; j++) {
                                        uint bs = (uint) (j * block_size);
                                        uint snb_off = count[j * num_buckets + i];
                                        uint ct;
                                        if (i == num_buckets - 1) {
                                            ct = (std::min((uint) (bs + block_size), local_nbr_off) - bs) - snb_off;
                                        } else {
                                            ct = count[j * num_buckets + i + 1] - snb_off;
                                        }

                                        snb_off += bs;

                                        for (uint l = 0; l < ct; l++) {
                                            v = sorted_nbr[snb_off + l];
                                            deg[v]++;
                                        }
                                    }
                                };
                            }

                        }

#pragma omp master
                        nbr_off = 0;
#pragma omp barrier

                    }

                    to_e = ree;
                }
            }

#pragma omp barrier
        }


#pragma omp barrier


#pragma omp for//nowait
        for (uint i = 0; i < sub_roots.size(); i++) {

            auto sub_r = &sub_roots[i];
            auto small_g = ComputeSmallG(mg, sub_r, sg_adj_buf, sg_ngh_buf, buf);

            sub_r->m = reinterpret_cast<uint *>(small_g);
        }


#pragma omp barrier


        /* ========== prepare ============*/
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg, buf, max_sr_size);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp barrier
#pragma omp master
        {

#pragma omp taskloop
            for (auto &sub_r: sub_roots) {

                auto g_small = reinterpret_cast<G_small *>(sub_r.m);

                // initialize new tasks
                uint inc_k = mlc_tree.GetIncK(sub_r.r);
                for (uint i = inc_k; i < ln - 1; i++) {
#pragma omp task
                    {
                        BuildRMPath(mlc_state, g_small, sub_r.r, i);
                    }
                }

                auto mlc_ss = mlc_state[omp_get_thread_num() << 6];

                mlc_ss->Set(g_small);
                mlc_ss->SetKVec(MLCTree::GetK(sub_r.r));
                mlc_ss->BatchBuildRMBranch(sub_r.r);
            }
        };

#pragma omp barrier
        delete[] buf;
    }




//     =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }

    for (auto &sub_r:sub_roots) {
        delete reinterpret_cast<G_small *>(sub_r.m);
    }


    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;
    delete[] node_queue;

    delete[] nbr_bk;
    delete[] count;
}

void PMLCTreeBuilder::PRun_async_all_upd(MultilayerGraph &mg, MLCTree &mlc_tree, uint level) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1, n_threads = omp_get_max_threads();
    uint k[ln], *re_queue, **degs, ***adj_lsts;
    uint re_top = 0, iv_off = 0, n_fin = 0, re, nbr_off = 0, compact_nbr_off = 0;
    uint num_buckets, num_blocks, block_size, le;
    uint k_sum = 1, max_sr_size = 0;

    MLCState **mlc_state;
    MDataBuf<uint> nm;
    DataBuf<uint> sg_ngh_buf;
    DataBuf<uint *> sg_adj_buf;

    PCoreIndex kci;
    Node *r, *child, **node_queue;
    Diff *diff;

    vector<SUB_ROOT> sub_roots;

    bool not_fin = false, nbr_sort, gen_gs;

    // Set mlc peeling state
    mlc_state = new MLCState *[n_threads << 6];

    bool fin = false;

    /* init structs */
    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, node_queue);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    r = mlc_tree.GetNewTreeNode(k, 0);
    mlc_tree.SetRoot(r);

    node_queue[re_top] = r;
    re_queue[re_top++] = 0;
    k[top]++;

    /* parallel run */
#pragma omp parallel
    {
        uint old_e, v, u, s, e, e_be, iv_off_pt, start, e1, local_nbr_off;
        uint *deg, *buf = new uint[n];

        /* init */
#pragma omp for schedule(static) nowait
        for (uint i = 0; i < n; i++) {
            kci.inva[i] = false;
            kci.vert[i] = i;
        }

#pragma omp for schedule(static) collapse(2)
        for (uint i = 0; i < ln; i++) {
            for (uint j = 0; j < n; j++) {
                degs[i][j] = adj_lsts[i][j][0];
            }
        }

        /* recursive compute */
        while (top >= 0) {

            /* compute ml-core */
            e = 0;
            deg = degs[top];
            old_e = n_fin;

#pragma omp master
            iv_off = 0;

//#pragma omp master
//            print_k(k, ln, k_sum, kci.e);

            /* compute first batch of removal */
#pragma omp for schedule(static) nowait
            for (uint j = old_e; j < kci.n; j++) {
                v = kci.vert[j];
                if (deg[v] < k[top] && !__atomic_test_and_set(&kci.inva[v], __ATOMIC_SEQ_CST)) {
                    buf[e++] = v;
                }
            }

            e_be = n_fin;
            /* peel */
            Peel(mg, kci, degs, k, buf, e);
            s = __sync_fetch_and_add(&kci.e, e);

#pragma omp barrier
#pragma omp master
            {

                /* prepare augmentation */
                if (top == ln - 1) {

                    diff = mlc_tree.GetDiff(r);
                    diff->num = kci.e - n_fin;
                    diff->vtx_ptr = mlc_tree.AllocateDiff(diff->num);

                } else {
                    diff = nullptr;
                }

                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    gen_gs = false;

                    /* continue exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);

                    if (top != ln - 1) {
                        top = ln - 1;

                        node_queue[re_top] = child;
                        re_queue[re_top++] = kci.e;
                    }


                    k[top]++;
                    r = child;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else if (kci.e < kci.n && k_sum == level) {

                    gen_gs = true;

                    /* stop exploration and generate child*/
                    child = mlc_tree.GetNewTreeNode(k, top);
                    mlc_tree.SetRelChd(r, top, child);


                    uint new_n = n - kci.e;
                    uint *m_addr = nm.Allocate(ln);
                    uint *relable = sg_ngh_buf.Allocate(new_n);
                    sub_roots.emplace_back(SUB_ROOT{child, new_n, m_addr, relable});

                    if (new_n > max_sr_size) max_sr_size = new_n;

                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        k[top] = 0;
                        k[--top]++;

                        /* restore */
                        re = re_queue[re_top - 1];
                        r = node_queue[re_top - 1];

                        if (mlc_tree.GetIncK(r) == top) {
                            re_top--;
                        }

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;
                    }


                } else {

                    gen_gs = false;

                    /* backtrace*/
                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        k[top] = 0;
                        k[--top]++;

                        /* restore */
                        re = re_queue[re_top - 1];
                        r = node_queue[re_top - 1];

                        if (mlc_tree.GetIncK(r) == top) {
                            re_top--;
                        }

                        kci.e = re;
                        kci.s = re;
                        n_fin = re;
                    }
                }
            }


#pragma omp barrier

            auto new_s = s - e_be;
            if (diff) {
                for (uint i = 0; i < e; i++) {
                    diff->vtx_ptr[new_s + i] = buf[i];
                }
            }

            if (le < n) {

                /* reorganize kci.vert */
                e1 = n;
#pragma omp for schedule(static)
                for (uint j = e_be; j < kci.n; j++) {
                    v = kci.vert[j];
                    if (!kci.inva[v]) {
                        buf[--e1] = v;
                    }
                }

                for (uint i = 0; i < e; i++) {
                    kci.vert[s + i] = buf[i];
                }

                iv_off_pt = le + __sync_fetch_and_add(&iv_off, n - e1) - e1;

                for (uint i = e1; i < n; i++) {
                    kci.vert[iv_off_pt + i] = buf[i];
                }

#pragma omp barrier
            }


            if (gen_gs) {
                auto &sub_r = sub_roots.back();

#pragma omp master
                {
                    memcpy(sub_r.label, kci.vert + le, sub_r.n * sizeof(uint));
                    assert(sub_r.n == n - le);
                };

#pragma omp for nowait
                for (uint i = 0; i < ln; i++) {
                    uint sum = 0;
                    for (uint j = le; j < n; j++) {
                        sum += degs[i][kci.vert[j]];
                    }
                    sub_r.m[i] = sum + sub_r.n;
                }
            }

            if (not_fin) {
                break;
            }

            if (re != INT32_MAX) {
                /* restore */
#pragma omp for nowait
                for (uint j = re; j < le; j++) {
                    v = kci.vert[j];
                    kci.inva[v] = false;
                }

#pragma omp for collapse(2)
                for (uint i = 0; i < ln; i++) {
                    for (uint j = re; j < le; j++) {
                        v = kci.vert[j];
                        for (uint l = 1; l <= adj_lsts[i][v][0]; l++) {
                            u = adj_lsts[i][v][l];
                            __sync_fetch_and_add(&degs[i][u], 1);
                        }
                    }
                }
            }


#pragma omp barrier
        }


#pragma omp barrier


#pragma omp for//nowait
        for (auto &sub_root : sub_roots) {

            auto sub_r = &sub_root;
            auto small_g = ComputeSmallG(mg, sub_r, sg_adj_buf, sg_ngh_buf, buf);

            sub_r->m = reinterpret_cast<uint *>(small_g);
        }


#pragma omp barrier


        /* ========== prepare ============*/
        uint off = omp_get_thread_num() << 6;
        auto mlc_s = new MLCState(mg, buf, max_sr_size);

        if (off) {
            mlc_s->mlc_t = new MLCTree(ln, n, false);
        } else {
            mlc_s->mlc_t = &mlc_tree;
        }

        mlc_state[off] = mlc_s;

#pragma omp barrier
#pragma omp master
        {

#pragma omp taskloop
            for (auto &sub_r: sub_roots) {

                auto g_small = reinterpret_cast<G_small *>(sub_r.m);

                // initialize new tasks
                uint inc_k = mlc_tree.GetIncK(sub_r.r);
                for (uint i = inc_k; i < ln - 1; i++) {
#pragma omp task
                    {
                        BuildRMPath(mlc_state, g_small, sub_r.r, i);
                    }
                }

                auto mlc_ss = mlc_state[omp_get_thread_num() << 6];

                mlc_ss->Set(g_small);
                mlc_ss->SetKVec(MLCTree::GetK(sub_r.r));
                mlc_ss->BatchBuildRMBranch(sub_r.r);
            }
        };

#pragma omp barrier
        delete[] buf;
    }




    //     =========== Merge MLCTree ===========
    if (n_threads > 1) MergeMLCTree(mlc_tree, mlc_state, n_threads);

    // =========== release space ============
    if (mlc_state) {
        for (int i = 0; i < n_threads; i++) {
            auto off = i << 6;
            if (i != 0) delete mlc_state[off]->mlc_t;
            delete mlc_state[off];
        }
        delete[] mlc_state;
    }

    for (auto &sub_r:sub_roots) {
        delete reinterpret_cast<G_small *>(sub_r.m);
    }


    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;
    delete[] node_queue;

}

void PMLCTreeBuilder::InitNbrRes(MultilayerGraph &mg, uint *&nbr, uint *&sorted_nbr, uint *&nbr_bk,
                                 uint *&count) {
    uint ln = mg.GetLayerNumber(), max_m = mg.GetGraph(0).GetM();

    for (uint i = 1; i < ln; i++) {
        max_m = std::max(max_m, mg.GetGraph(i).GetM());
    }

    nbr_bk = new uint[(max_m << 1) + max_m];
    nbr = nbr_bk + max_m;
    sorted_nbr = nbr + max_m;
    count = new uint[get_n_blk(max_m) * get_n_bkt(max_m)];
}


void PMLCTreeBuilder::GenOneHopeNgh(Graph &g, const bool *upd, const uint *buf, uint start,
                                    uint e, uint *nbr, uint *sorted_nbr, uint &nbr_off, uint &compact_nbr_off) {
    {

        uint **adj_lst, off, local_nbr_off;
        uint v, u, total_d = 0;

        adj_lst = g.GetAdjLst();   /* only modify the degs of the vertices not being removed */

        for (uint j = start; j < e; j++) {
            total_d += adj_lst[buf[j]][0];
        }

        off = __sync_fetch_and_add(&compact_nbr_off, total_d);

        local_nbr_off = 0;

        for (uint j = start; j < e; j++) {
            v = buf[j];
            for (uint l = 1; l <= adj_lst[v][0]; l++) {
                u = adj_lst[v][l];
                if (upd[u]) sorted_nbr[off + local_nbr_off++] = u;
            }
        }

        auto compact_off = __sync_fetch_and_add(&nbr_off, local_nbr_off);
        for (uint i = 0; i < local_nbr_off; i++) {
            nbr[compact_off + i] = sorted_nbr[off + i];
        }
    }
}


void PMLCTreeBuilder::SegSort(uint *nbr, uint *nbr_bk, uint *sorted_nbr, uint &nbr_off, uint *count, uint num_buckets,
                              uint num_blocks, uint block_size) {

#pragma omp taskloop
    for (uint i = 0; i < num_blocks; i++) {

        auto offset = i * block_size;
        count_sort(nbr + offset, nbr_bk + offset, sorted_nbr + offset, std::min(block_size, nbr_off - offset),
                   count + i * num_buckets, num_buckets);
    }
}


void PMLCTreeBuilder::count_sort(uint *nbr, uint *nbr_bk, uint *sorted_nbr, uint n, uint *count, uint num_buckets) {
    {
        size_t low_mask = ~((size_t) 15);
        size_t bucket_mask = num_buckets - 1;

        auto get_bkt = [&](uint i) {
            return hash32(i & low_mask) & bucket_mask;
            //            return (i & low_mask) & bucket_mask;
        };

        memset(count, 0, num_buckets * sizeof(uint));


        for (uint j = 0; j < n; j++) {
            nbr_bk[j] = get_bkt(nbr[j]);
            count[nbr_bk[j]]++;
        }


        for (uint j = 1; j < num_buckets; j++) {
            count[j] += count[j - 1];
        }

        for (uint j = 0; j < n; j++) {
            sorted_nbr[--count[nbr_bk[j]]] = nbr[j];
        }
    }
}

void PMLCTreeBuilder::Peel(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf,
                           uint &e) {

    uint ln = mg.GetLayerNumber(), old_e, s, v, u, du, **adj_lst;
    int t = omp_get_thread_num();

    s = 0;
    while (s < e) {
        old_e = e;

        for (uint layer = t; layer < ln + t; layer++) {

            auto i = layer % ln;
            adj_lst = mg.GetGraph(i).GetAdjLst();

            for (uint j = s; j < old_e; j++) {
                v = buf[j];

                for (uint l = 1; l <= adj_lst[v][0]; l++) {
                    u = adj_lst[v][l];

                    du = __sync_fetch_and_sub(&degs[i][u], 1);

                    if (du == k[i]) {
                        if (!__atomic_test_and_set(&kci.inva[u], __ATOMIC_SEQ_CST)) {
                            buf[e++] = u;
                        }
                    }
                }
            }
        }

        s = old_e;
    }
}