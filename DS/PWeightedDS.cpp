//
// Created by ldd on 2023/10/30.
//

#include "PWeightedDS.h"

float
PWeightedDS::Search(MultilayerGraph &mg, float beta, const vector<float> &w, MLCTree &mlc_t, uint *ds, uint &length,
                    vector<uint> &k) {

    uint ln = mg.GetLayerNumber(), n_threads = omp_get_max_threads(), *ds_k;
    float max_den = 0;

    WDS_search_state *ws_states[n_threads << 6];
    Node *ds_node = nullptr;


#pragma omp parallel
    {
        auto *wds_s = new WDS_search_state(&mg, w.data(), beta, &mlc_t);
        ws_states[omp_get_thread_num() << 6] = wds_s;
#pragma omp barrier
#pragma omp master
        {
#pragma omp task
            RecursiveSearch(&mg, ws_states, mlc_t.GetRoot(), &mlc_t);
        };
    };

    for (int i = 0; i < n_threads; i++) {
        auto wds_s = ws_states[i << 6];

        if (wds_s->max_den > max_den) {
            max_den = wds_s->max_den;
            ds_node = wds_s->max_den_node;
        }
    }

    if (ds_node) {
        ds_k = MLCTree::GetK(ds_node);
        k.assign(ds_k, ds_k + ln);

        length = mlc_t.Retrieve(ds_node, ds);
    }

    for (int i = 0; i < n_threads; i++) {
        delete ws_states[i << 6];
    }
    return max_den;
}

float
PWeightedDS::Search(float beta, const vector<float> &w, DSMLCTree &mlc_t, uint *ds, uint &length,
                    vector<uint> &k) {

    uint ln = mlc_t.GetLayerNumber(), n_threads = omp_get_max_threads(), *ds_k;
    float max_den = 0, bl[ln];
    Node *ds_node = nullptr;

    for (uint i = 0; i < ln; i++) {
        bl[i] = (float) pow(ln - i, beta);
    }

    WDS_state_simp *wds_states[n_threads << 6];

#pragma omp parallel
    {
        auto *wds_state = new WDS_state_simp(ln, w.data(), bl, &mlc_t);
        wds_states[omp_get_thread_num() << 6] = wds_state;

#pragma omp barrier
#pragma omp master
        {
#pragma omp task
            SearchRm(mlc_t.GetRoot(), wds_states);
        };
    };

    for (int i = 0; i < n_threads; i++) {
        auto wds_s = wds_states[i << 6];

        if (wds_s->max_den > max_den) {
            max_den = wds_s->max_den;
            ds_node = wds_s->max_den_node;
        }
    }

    if (ds_node) {
        ds_k = MLCTree::GetK(ds_node);
        k.assign(ds_k, ds_k + ln);

        length = mlc_t.Retrieve(ds_node, ds);
    }

    return max_den;
}

float
PWeightedDS::Search(MultilayerGraph &mg, float beta, const vector<float> &w, uint *ds, uint &length,
                    vector<uint> &k_vec, uint level, float merge_factor) {

    uint ln = mg.GetLayerNumber(), n = mg.GetN(), top = ln - 1, n_threads = omp_get_max_threads();
    uint re_top = 0, iv_off = 0, n_fin = 0, re, le, k_sum = 1, max_sr_size = 0, max_den_t;
    uint k[ln], n_edges[ln], *re_queue, **degs, ***adj_lsts;
//    uint n_nodes = 1;

    level = 2;
    merge_factor = 0.1;

    WDS_peeling_state *wds_state[n_threads << 6];
    float max_den = 0, wd[ln], bl[ln];
    float den, max_ww = 0;

    DataBuf<uint> sg_ngh_buf;
    DataBuf<uint *> sg_adj_buf;

    PCoreIndex kci;
    vector<SUB_ROOT_K> sub_roots;

    bool not_fin = false, gen_gs;

    // Set mlc peeling state
    sub_roots.reserve(ln);
    k_vec.resize(ln);

    /* init structs */
    kci.Init(n, true);
    InitDeg(mg, degs, adj_lsts);
    InitQueue(mg, re_queue, level);
    memset(k, 0, ln * sizeof(uint));

    /* prepare for dfs */
    re_queue[re_top++] = 0;
    k[top]++;

    for (uint i = 0; i < ln; i++) {
        bl[i] = (float) pow(ln - i, beta);
    }

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
//            print_k(k, ln);

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
                le = kci.e;

                /* gen next k */
                if (kci.e < kci.n && k_sum < level) {

                    gen_gs = false;
//                    n_nodes++;

                    /* continue exploration and generate child*/

                    if (top != ln - 1) top = ln - 1;
                    re_queue[re_top++] = kci.e;

                    k[top]++;
                    n_fin = kci.e;
                    re = INT32_MAX;

                    k_sum += 1;

                } else if (kci.e < kci.n && k_sum == level) {

                    gen_gs = true;
//                    n_nodes++;

                    /* stop exploration and generate child*/
                    uint new_n = n - kci.e;
                    uint *curr_k = sg_ngh_buf.Allocate(ln);
                    memcpy(curr_k, k, ln * sizeof(uint));

                    uint *m_addr = sg_ngh_buf.Allocate(ln);
                    uint *relable = sg_ngh_buf.Allocate(new_n);
                    sub_roots.emplace_back(SUB_ROOT_K{curr_k, top, new_n, m_addr, relable});

                    if (new_n > max_sr_size) max_sr_size = new_n;

                    if (top == 0) {
                        not_fin = true;  /* return */

                    } else {

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];

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

                        memset(n_edges, 0, ln * sizeof(uint));

                        k_sum = k_sum - k[top] + 1;

                        re_top -= k[top];
                        re = re_queue[re_top];

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
                    n_edges[i] = sum >> 1;
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
                        kci.upd[kci.vert[j]] = true;
                    }


#pragma omp for schedule(static)
                    for (uint l = 0; l < ln; l++) {
                        uint inc = 0;
                        for (uint j = ree; j < to_e; j++) {
                            v = kci.vert[j];
                            for (uint r = 1; r <= adj_lsts[l][v][0]; r++) {
                                u = adj_lsts[l][v][r];
                                if (!kci.inva[u] || (kci.upd[u] && u > v)) {
                                    inc++;
                                }
                            }
                        }
                        n_edges[l] += inc;
                        wd[l] = (float) n_edges[l] / (float) (n - ree) * w[l];
                    }
#pragma omp master
                    {
                        std::sort(wd, wd + ln);
                        for (int i = (int) ln - 1; i >= 0; i--) {
                            den = wd[i] * bl[i];

                            if (den > max_ww) {
                                max_ww = den;
                            }
                        }

                        if (max_ww > max_den) {
                            max_den = max_ww;
                            length = (n - ree);
                            memcpy(ds, kci.vert + ree, length * sizeof(uint));
                            memcpy(k_vec.data(), k, ln * sizeof(uint));
                            k_vec[ln - 1] = kk;
                            k_vec[top]--;

                        }
                    };

#pragma omp for schedule(static) nowait
                    for (uint j = ree; j < to_e; j++) {
                        v = kci.vert[j];
                        kci.inva[v] = false;
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
        MDataBuf<uint> k_buf;
        auto wds_s = new WDS_peeling_state(&mg, &k_buf, w.data(), beta, max_den, buf, max_sr_size);

        wds_state[off] = wds_s;

#pragma omp barrier
#pragma omp master
        {

#pragma omp taskloop
            for (auto &sub_r: sub_roots) {

                auto g_small = reinterpret_cast<G_small *>(sub_r.m);

                // initialize new tasks

                for (uint i = sub_r.inc_k; i < ln - 1; i++) {
#pragma omp task
                    {
                        BuildRMPath_m(sub_r.k, wds_state, g_small, i, merge_factor);
                    }
                }

                auto wds_ss = wds_state[omp_get_thread_num() << 6];

                wds_ss->Set(g_small);
                wds_ss->SetKVec(sub_r.k);
                wds_ss->Clear();
                wds_ss->SearchRMBranch();
            }
        };

#pragma omp barrier
        delete[] buf;
    }

    // =========== release space ============


//    for (int i = 0; i < n_threads; i++) {
//        n_nodes += wds_state[i << 6]->n_nodes;
//    }

    max_den_t = -1;
    for (int i = 0; i < n_threads; i++) {
        auto wds_s = wds_state[i << 6];

        if (wds_s->max_den > max_den) {
            max_den = wds_s->max_den;
            max_den_t = i;
        }
    }

    if (max_den_t != -1) {
        auto wds_s = wds_state[max_den_t << 6];
        length = wds_s->length;
        k_vec.assign(wds_s->max_den_k, wds_s->max_den_k + ln);

        for (uint i = 0; i < length; i++) {
            ds[i] = wds_s->max_den_relabel[wds_s->core[i]];
        }
    }

    /* parallel computing ml-core*/

    for (auto &sub_r:sub_roots) {
        delete reinterpret_cast<G_small *>(sub_r.m);
    }

    for (int i = 0; i < n_threads; i++) {
        delete wds_state[i << 6];
    }

    /* free*/
    for (uint i = 0; i < ln; i++) delete[] degs[i];
    delete[] degs;
    delete[] adj_lsts;
    delete[] re_queue;

    sg_ngh_buf.Release();
    sg_adj_buf.Release();

    return max_den;

}

void PWeightedDS::RecursiveSearch(MultilayerGraph *mg, WDS_search_state **wds_states, Node *r, MLCTree *mlc_t) {

    uint inc_k = mlc_t->GetIncK(r), max_d = mg->GetLayerNumber() - 1;
    Node *child;
    WDS_search_state *wds_s;

    for (uint i = inc_k; i < max_d; i++) {
        child = mlc_t->GetRelChd(r, i);

        if (child) {
#pragma omp task
            RecursiveSearch(mg, wds_states, child, mlc_t);
        }
    }

    wds_s = wds_states[omp_get_thread_num() << 6];
    wds_s->Clear();
    wds_s->SearchRMPath(r);
}

void PWeightedDS::SearchRm(Node *r, WDS_state_simp **wds_states) {
    WDS_state_simp* wds_state = wds_states[omp_get_thread_num() << 6];
    DSMLCTree* mlc_t = wds_state->mlc_tree;

    uint inc_k = mlc_t->GetIncK(r), max_d = wds_state->ln - 1;
    Node *child;

    for (uint i = inc_k; i < max_d; i++) {
        child = mlc_t->GetRelChd(r, i);

        if (child) {
#pragma omp task
            SearchRm(child, wds_states);
        }
    }

    wds_state->Clear();
    wds_state->RecSearch(r);
}

void PWeightedDS::InitDeg(MultilayerGraph &mg, uint **&degs, uint ***&adj_lsts) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN();

    degs = new uint *[ln];
    adj_lsts = new uint **[ln];
    for (uint i = 0; i < ln; i++) {
        degs[i] = new uint[n];
        adj_lsts[i] = mg.GetGraph(i).GetAdjLst();
    }
}

void PWeightedDS::InitQueue(MultilayerGraph &mg, uint *&re_queue, uint level) {
    uint ln = mg.GetLayerNumber(), deg_sum = ln;
    for (uint i = 0; i < ln; i++) {
        deg_sum += mg.GetGraph(i).GetMaxDeg();
    }

    re_queue = new uint[std::min(level, deg_sum)];
}

void PWeightedDS::PeelPud(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf,
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

G_small *
PWeightedDS::ComputeSmallG(MultilayerGraph &mg, SUB_ROOT_K *sub_r, DataBuf<uint *> &adj_buf, DataBuf<uint> &nbr_buf,
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

void PWeightedDS::BuildRMPath_m(uint *curr_k, WDS_peeling_state **wds_state, G_small *ref_g, uint inc_k,
                                float merge_factor) {
    uint ld, sg_size;

    WDS_peeling_state *wds_s = wds_state[omp_get_thread_num() << 6];

    wds_s->Set(ref_g);
    wds_s->SetKVec(curr_k);
    wds_s->k_vec[inc_k] += 1;

    ld = wds_s->ld;

    sg_size = wds_s->PeelCore();

    if (sg_size) {

//        wds_s->n_nodes++;

        if (sg_size <= (uint) ((float) ref_g->GetN() * merge_factor)) {
            wds_s->Clear();
            wds_s->RecursiveSearchChild(inc_k);

        } else {

            wds_s->Clear();
            wds_s->SearchRMBranch();

            for (auto i = inc_k; i < ld; i++) {
                auto child_k = wds_s->GetDupK();
#pragma omp task
                BuildRMPath_m(child_k, wds_state, ref_g, i, merge_factor);
            }
        }
    }
}
