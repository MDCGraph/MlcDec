//
// Created by ldd on 2023/6/16.
//

#ifndef MLCDEC_MLCTBUILDER_H
#define MLCDEC_MLCTBUILDER_H

#include "../Header.h"

enum MLCTree_builder {
    NAIVE, NODE_ELIMINATION, SUBTREE_ELIMINATION, SUBTREE_MERGE, OPT
};


static string MlCTreeBuilder2Str(MLCTree_builder builder) {
    if (builder == NAIVE) return "NAIVE";
    else if (builder == NODE_ELIMINATION) return "NE";
    return "";
}

#endif //MLCDEC_MLCTBUILDER_H
