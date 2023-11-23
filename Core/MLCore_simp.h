//
// Created by ldd on 2023/10/28.
//

#ifndef MLCDEC_MLCORE_SIMP_H
#define MLCDEC_MLCORE_SIMP_H

#include "../Graphs/MultilayerGraph.h"
#include "CoreIndex.h"

class MLCore_simp {
public:

    static uint Extract(MultilayerGraph& mg, const vector<uint> &k_vec_, uint *mlc);

private:

    static void Peel(MultilayerGraph& mg, CoreIndex &kci, uint *k_vec, uint** degs);
};

#endif //MLCDEC_MLCORE_SIMP_H
