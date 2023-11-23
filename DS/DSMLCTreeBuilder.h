//
// Created by ldd on 2023/11/5.
//

#ifndef MLCDEC_DSMLCTREEBUILDER_H
#define MLCDEC_DSMLCTREEBUILDER_H

#include "DSMLCTree.h"
#include "../Graphs/MultilayerGraph.h"
#include "../Core/CoreIndex.h"

class DSMLCTreeBuilder {
public:
    static void Execute(MultilayerGraph & mg, DSMLCTree &mlc_tree);

private:
    static void BuildSubMLCTree(MultilayerGraph& mg, DSMLCTree &mlc_tre, CoreIndex& kci, uint** degs, Node* r, uint* k, uint inc_k);
    static void Peel(MultilayerGraph& mg, CoreIndex& kci, uint** degs, const uint* k);
    static void RestoreRm(MultilayerGraph& mg, CoreIndex& kci, uint** degs, uint old_e, uint* diff_nc);
    static void Restore(MultilayerGraph &mg, CoreIndex &kci, uint **degs, uint old_e);

};

#endif //MLCDEC_DSMLCTREEBUILDER_H
