//
// Created by ldd on 2023/10/23.
//

#ifndef MLCDEC_PMLCTREEBUILDER_H
#define MLCDEC_PMLCTREEBUILDER_H

#include "MLCTBuilder.h"
#include "MLCState.h"
#include "MLCTree.h"
#include "MLCDfs.h"
#include "PCoreIndex.h"
#include "../Utils/Constructs.h"

#include "MLCGSmall.h"

struct SUB_ROOT {
    Node *r;
    uint n;
    uint *m;
    uint *label;

    SUB_ROOT(Node *r_, uint n_, uint *m_, uint *label_) : r(r_), n(n_), m(m_), label(label_) {}
};

class PMLCTreeBuilder {

public:

    /* use */
    static void InterPathPRun(MultilayerGraph &mg, MLCTree &mlc_tree);
    static void InterPathPRun_with_merge(MultilayerGraph &mg, MLCTree &mlc_tree, float merge_factor);

    static void PRun_async(MultilayerGraph &mg, MLCTree &mlc_tree, uint level);
    static void PRun_async_with_merge(MultilayerGraph &mg, MLCTree &mlc_tree, uint level, float merge_factor);

    /* not in use*/
    static void InnerCorePRun(MultilayerGraph &mg, MLCTree &mlc_tree, uint level = UINT32_MAX);
    static void PRun(MultilayerGraph &mg, MLCTree &mlc_tree, uint level);
    static void PRun_async_all_upd(MultilayerGraph &mg, MLCTree &mlc_tree, uint level);

private:

    static inline void InitDeg(MultilayerGraph &mg, uint **&degs, uint ***&adj_lsts);
    static inline void InitQueue(MultilayerGraph &mg, uint *&re_queue, Node **&node_queue);

    static void PeelPud(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf, uint &e);

    static void BuildRMPath(MLCState **mlc_state, Node *r, uint inc_k);
    static void BuildRMPath(MLCState **mlc_state, G_small *ref_g, Node *r, uint inc_k);

    static void BuildRMPath_m(MLCState **mlc_state, Node *r, uint inc_k, uint merge_bound);
    static void BuildRMPath_m(MLCState **mlc_state, G_small *ref_g, Node *r, uint inc_k, float merge_factor);

    static G_small *
    ComputeSmallG(MultilayerGraph &mg, SUB_ROOT *sub_r, DataBuf<uint *> &adj_buf, DataBuf<uint> &nbr_buf, uint *buf);

    static void MergeMLCTree(MLCTree &mlc_tree, MLCState **mlc_state, uint n);

    /* for all upd */
    static void Peel(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf, uint &e);

    /* for sync */
    static void
    GenOneHopeNgh(Graph &g, const bool *upd, const uint *buf, uint start, uint e, uint *nbr, uint *sorted_nbr,
                  uint &nbr_off, uint &compact_nbr_off);
    static void
    SegSort(uint *nbr, uint *nbr_bk, uint *sorted_nbr, uint &nbr_off, uint *count, uint num_buckets, uint num_blocks,
            uint block_size);

    static void count_sort(uint *nbr, uint *nbr_bk, uint *sorted_nbr, uint n, uint *count, uint num_buckets);

    static inline void InitNbrRes(MultilayerGraph &mg, uint *&nbr, uint *&sorted_nbr, uint *&nbr_bk, uint *&count);

    static inline uint get_n_bkt(uint n) {

        uint num_buckets = (size_t) ceil(pow(n, 0.5)) / 5;
        num_buckets = std::min(std::max(1 << log2_up(num_buckets), 1), 1024);
        return num_buckets;
    }

    static inline uint get_n_blk(uint n) {

        uint num_blocks = (size_t) ceil(pow(n, 0.5)) / 10;
        num_blocks = 1 << log2_up(std::min(num_blocks, (uint) 512));
        return num_blocks;
    }
    /* for debug */

    static void print_k(uint *k, uint ln, uint top, uint e, uint re_top = -1) {
        if (re_top != -1)
            print(
                    Arr2Str(k, ln) + " top = " + std::to_string(top) + " e = " + std::to_string(e) + " re_top = " +
                    std::to_string(re_top));
        else print(Arr2Str(k, ln) + " top = " + std::to_string(top) + " e = " + std::to_string(e));
        fflush(stdout);
    }

    static void print_k(uint *k, uint ln) {
        print(Arr2Str(k, ln));
        fflush(stdout);
    }


};


#endif //MLCDEC_PMLCTREEBUILDER_H
