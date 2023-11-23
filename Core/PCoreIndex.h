//
// Created by ldd on 2023/9/21.
//

#ifndef MLCDEC_PCOREINDEX_H
#define MLCDEC_PCOREINDEX_H


#include "../Graphs/Graph.h"

/*
 * ============= PCoreIndex =============
 * Integer n represents the total number of vertices.
 * Array vert contains all vertices.
 * Array inva contains indicates whether each vertex is invalid
 * Offsets s and e divide vert into 3 parts:
 * |***********|****************|**************************|
 * 0-----------sta----------------e--------------------------n
 * Part 1: vert[0:s) contains all discarded vertices.
 * Part 2: vert[s:e) contains all inactive vertices.
 * Part 3: vert[e:n) contains all active vertices.
 */

struct PCoreIndex {

    uint n{0};
    uint s{0};
    uint e{0};

    uint *vert{nullptr};

    bool *inva{nullptr};
    bool *upd{nullptr};

    PCoreIndex() = default;

    ~PCoreIndex() {
        delete[] vert;
        delete[] inva;
        delete[] upd;
    }

    void Init(uint n_, bool use_update = false) {
        n = n_;
        vert = new uint[n];
        inva = new bool[n];

        if (use_update) upd = new bool[n];
    }

    void Set() {
        s = 0;
        e = 0;
        for (uint i = 0; i < n; i++) {
            vert[i] = i;
            inva[i] = false;
        }
    }

};


#endif //MLCDEC_PCOREINDEX_H
