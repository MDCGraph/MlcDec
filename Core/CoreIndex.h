//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_COREINDEX_H
#define MLCDEC_COREINDEX_H

#include "../Graphs/Graph.h"


/*
 * ============= CoreIndex =============
 * Integer n represents the total number of vertices.
 * Array vert contains all vertices.
 * Array pos contains the positions of vertices in array vert.
 * Offsets sta and t divide vert into 3 parts:
 * |***********|****************|**************************|
 * 0-----------sta----------------e--------------------------n
 * Part 1: vert[0:s) contains all discarded vertices.
 * Part 2: vert[s:e) contains all inactive vertices.
 * Part 3: vert[e:n) contains all active vertices.
 */

struct CoreIndex {

    bool released{true};

    uint n{0};
    uint s{0};
    uint e{0};

    uint *pos{nullptr};
    uint *vert{nullptr};

    CoreIndex() = default;

    ~CoreIndex() {
        if (released) {
            delete[] pos;
            delete[] vert;
        }
    }

    void Init(uint n_) {
        n = n_;
        pos = new uint[n];
        vert = new uint[n];
    }

    void Init(uint n_, uint *buf) {
        released = false;

        n = n_;
        pos = buf;
        vert = buf + n;
    }

    // SetMLCS CoreIndex with n vertices.
    // All vertices are active initially.
    void Set() {
        s = 0;
        e = 0;
        for (uint i = 0; i < n; i++) {
            pos[i] = i;
            vert[i] = i;
        }
    }

    // SetMLCS CoreIndex with n vertices.
    // All vertices are active initially.
    void Set(uint n_) {
        s = 0;
        e = 0;
        n = n_;

        for (uint i = 0; i < n; i++) {
            pos[i] = i;
            vert[i] = i;
        }
    }

    void Copy(CoreIndex &ci, uint offset) {
        s = ci.s;
        e = ci.e;
        memcpy(&vert[offset], &ci.vert[offset], (n - offset) * sizeof(uint));
        memcpy(pos, ci.pos, n * sizeof(uint));
    }

    // Move v from Part 3 to Part 2, i.e., convert v from active state to inactive state.
    inline void Remove(uint v) {
        uint v_pos;

        v_pos = pos[v];
        vert[v_pos] = vert[e];
        vert[e] = v;
        pos[vert[v_pos]] = v_pos;
        pos[v] = e++;
    }

};


#endif //MLCDEC_COREINDEX_H
