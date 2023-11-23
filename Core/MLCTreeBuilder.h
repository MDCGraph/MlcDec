//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_MLCTREEBUILDER_H
#define MLCDEC_MLCTREEBUILDER_H

#include "MLCTBuilder.h"
#include "MLCTree.h"
#include "MLCDfs.h"

class MLCTreeBuilder : public MLCDfs{
public:
    explicit MLCTreeBuilder(MultilayerGraph &mg_);
    ~MLCTreeBuilder() = default;
    void Execute(MLCTree &mlc_tree_, MLCTree_builder opt = OPT);

private:

    MLCTree *mlc_tree{};
    void BuildSubMLCTreeNaive(uint inc_k, Node *r);
    void BuildSubMLCTreeNE(uint inc_k, Node *r, uint * f_sign);
    void BuildSubMLCTreeSE(uint inc_k, Node *r, uint * f_sign);
    void BuildSubMLCTreeSM(uint inc_k, Node *r, Node *anc, Node *pr, uint * f_sign, unordered_map<Node*, uint *> &signs);
    void BuildSubMLCTreeSESM(uint inc_k, Node *r, Node *anc, Node *pr, uint * f_sign, unordered_map<Node*, uint *> &signs);
    bool BuildBranch(uint i, Node *r);
};

#endif //MLCDEC_MLCTREEBUILDER_H
