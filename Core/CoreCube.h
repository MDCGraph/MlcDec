//
// Created by ldd on 2023/3/29.
//

#ifndef MLCDEC_CORECUBE_H
#define MLCDEC_CORECUBE_H

#include "../Graphs/MultilayerGraph.h"
#include "../Utils/ArrayUtils.h"
#include "Cube.h"

class CoreCube {
public:
    explicit CoreCube(MultilayerGraph &mg_);
    ~CoreCube();

    void TD(Cube &cube_);
    void GenHybridStorage(HYB_Cube &hyb_cube);

private:

    Cube *cube{};
    MultilayerGraph &mg;

    uint ln;
    uint n;

    uint z{};
    bool *layer_idt;
    uint *selected_layer;

    uint **sups;
    uint *q1;
    uint *q2;
    uint *counter;

    bool *in_q1;
    bool *in_q2;

    void InitLInd();
    void SetSLayer();

    void CoreTD();

    uint GetCurrCoreOffset();

    [[nodiscard]] uint GetAnsCoreOffset(uint mask) const;

    [[nodiscard]] uint GetNCombs() const {
        return (uint) pow(2, ln) - 1;
    }

    void StoreAbs(uint *coreness, uint n_pos, HYB_Cube1d &hyb_cube1d) const;
    void StoreRef(uint *coreness, uint ref_id, uint n_diff, HYB_Cube1d &hyb_cube1d) const;

};


#endif //MLCDEC_CORECUBE_H
