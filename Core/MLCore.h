//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_MLCORE_H
#define MLCDEC_MLCORE_H

#include "../Graphs/MultilayerGraph.h"
#include "CoreIndex.h"

enum PEELING_MODE {
    BY_LAYER, BY_VTX
};

class MLCore {
public:
    explicit MLCore(MultilayerGraph &mg_);
    ~MLCore();

    uint Extract(const vector<uint> &k_vec_, uint *mlc, PEELING_MODE peeling_mode = BY_LAYER);

private:
    MultilayerGraph &mg;
    uint * k_vec;

    uint **degs;
    uint ln;
    uint n;

    void PeelByLayer(CoreIndex & ci);
    void PeelByVtx(CoreIndex & ci);
};

#endif //MLCDEC_MLCORE_H
