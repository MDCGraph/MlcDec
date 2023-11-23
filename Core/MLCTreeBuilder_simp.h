//
// Created by ldd on 2023/9/14.
//

#ifndef MLCDEC_MLCTREEBUILDER_SIMP_H
#define MLCDEC_MLCTREEBUILDER_SIMP_H

#include "MLCTBuilder.h"
#include "MLCTree.h"
#include "MLCDfs.h"
#include "MLCHashTable.h"

class MLCTreeBuilder_simp {
public:
    static void Execute(MultilayerGraph & mg, MLCTree &mlc_tree, uint level = 0);
    static void Execute(MultilayerGraph & mg, MLCHashTable& mlc_ht);


private:
    static void BuildSubMLCTree(MultilayerGraph& mg, MLCTree &mlc_tre, CoreIndex& kci, uint** degs, Node* r, uint* k, uint inc_k);
    static void BuildSubMLCTree(MultilayerGraph& mg, MLCTree &mlc_tre, CoreIndex& kci, uint** degs, Node* r, uint* k, uint inc_k, uint balance);
    static void BuildSubMLCHt(MultilayerGraph& mg, MLCHashTable &mlc_ht, CoreIndex& kci, uint** degs, uint* k, uint inc_k);

    static void Peel(MultilayerGraph& mg, CoreIndex& kci, uint** degs, const uint* k);
    static void Restore(MultilayerGraph& mg, CoreIndex& kci, uint** degs, uint old_e);

};

#endif //MLCDEC_MLCTREEBUILDER_SIMP_H
