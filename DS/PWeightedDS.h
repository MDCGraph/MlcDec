//
// Created by ldd on 2023/10/30.
//

#ifndef MLCDEC_PWEIGHTEDDS_H
#define MLCDEC_PWEIGHTEDDS_H

#include "../Core/CoreCube.h"
#include "../Core/MLCTree.h"
#include "../Core/PCoreIndex.h"
#include "../Core/PMLCTreeBuilder.h"

#include "WDSPeelingState.h"
#include "Metric.h"
#include "DSMLCTree.h"

struct SUB_ROOT_K {
    uint* k;
    uint inc_k;
    uint n;
    uint *m;
    uint *label;

    SUB_ROOT_K(uint* k_, uint inc_k_, uint n_, uint *m_, uint *label_) : k(k_), inc_k(inc_k_), n(n_), m(m_), label(label_) {}
};


class PWeightedDS {
public:
    static float Search(MultilayerGraph &mg, float beta, const vector<float> &w, MLCTree &mlc_t, uint *ds, uint &length,
                        vector<uint> &k);

    static float Search(float beta, const vector<float> &w, DSMLCTree &mlc_t, uint *ds, uint &length,
                        vector<uint> &k);

    static float Search(MultilayerGraph &mg, float beta, const vector<float> &w, uint *ds, uint &length,
                        vector<uint> &k, uint level, float merge_factor);

private:
    static void RecursiveSearch(MultilayerGraph *mg, WDS_search_state **wds_states, Node *r, MLCTree *mlc_t);

    static void SearchRm(Node *r, WDS_state_simp **wds_state);

    static inline void InitDeg(MultilayerGraph &mg, uint **&degs, uint ***&adj_lsts);
    static inline void InitQueue(MultilayerGraph &mg, uint *&re_queue, uint level);

    static void PeelPud(MultilayerGraph &mg, PCoreIndex &kci, uint **degs, const uint *k, uint *buf,
                        uint &e);
    static G_small* ComputeSmallG(MultilayerGraph &mg, SUB_ROOT_K *sub_r, DataBuf<uint *> &adj_buf, DataBuf<uint> &nbr_buf,
                                  uint *buf);
    static void BuildRMPath_m(uint* curr_k, WDS_peeling_state **wds_state, G_small *ref_g, uint inc_k, float merge_factor);

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


#endif //MLCDEC_PWEIGHTEDDS_H
