//
// Created by ldd on 2023/3/7.
//

#ifndef MLCDEC_KC_H
#define MLCDEC_KC_H

#include "../Graphs/Graph.h"

class KC {
public:
    static uint CoreDec(Graph &g, uint *core_number);
    static uint Degeneracy(Graph & g);
};


#endif //MLCDEC_KC_H
