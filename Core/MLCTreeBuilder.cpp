//
// Created by ldd on 2023/3/6.
//

#include "MLCTreeBuilder.h"

MLCTreeBuilder::MLCTreeBuilder(MultilayerGraph &mg_) : MLCDfs(mg_) {}

void MLCTreeBuilder::Execute(MLCTree &mlc_tree_, MLCTree_builder opt) {
    Node *root;

    mlc_tree = &mlc_tree_;

    Init();
    memset(k_vec, 0, ln * sizeof(uint));
    root = mlc_tree->GetNewTreeNode(k_vec, 0);
    mlc_tree->SetRoot(root);

    if (opt == NAIVE) BuildSubMLCTreeNaive(0, root);
    else if (opt == NODE_ELIMINATION) {
        uint subtree_sign[ln];
        memset(subtree_sign, -1, ln* sizeof(uint));

        BuildSubMLCTreeNE(0, root, subtree_sign);
    } else if (opt == SUBTREE_ELIMINATION) {
        uint subtree_sign[ln];
        memset(subtree_sign, -1, ln* sizeof(uint));

        BuildSubMLCTreeSE(0, root, subtree_sign);
    } else if (opt == SUBTREE_MERGE) {
        uint subtree_sign[ln];
        unordered_map<Node *, uint *> signs;
        memset(subtree_sign, -1, ln* sizeof(uint));

        BuildSubMLCTreeSM(0, root, nullptr, nullptr, subtree_sign, signs);
        for (auto node_sign:signs) {
            delete[] node_sign.second;
        }
    } else if (opt == OPT) {
        uint subtree_sign[ln];
        unordered_map<Node *, uint *> signs;
        memset(subtree_sign, -1, ln* sizeof(uint));

        BuildSubMLCTreeSESM(0, root, nullptr, nullptr, subtree_sign, signs);
        for (auto node_sign:signs) {
            delete[] node_sign.second;
        }
    }
}


void MLCTreeBuilder::BuildSubMLCTreeNaive(uint inc_k, Node *r) {
    bool has_child;
    uint old_e;
    Node *child;

    old_e = k_core_index.e;

    for (auto i = (int) ln - 1; i >= (int) inc_k; i--) {

        has_child = BuildBranch(i, r);

        if (has_child) {
            child = mlc_tree->GetNewTreeNode(k_vec, i);
            mlc_tree->SetRelChd(r, i, child);
            BuildSubMLCTreeNaive(i, child);
        }

        k_vec[i]--;
        Restore(old_e);
    }
}

void MLCTreeBuilder::BuildSubMLCTreeNE(uint inc_k, Node *r, uint *f_sign) {
    bool has_child;
    uint old_e, old_k, *k, sign[ln];
    Node *child;

    old_e = k_core_index.e;
    k = MLCTree::GetK(r);

    for (uint i = 0; i <= inc_k; i++) {
        sign[i] = k_value_heaps[i].GetMinValue();
    }
    for (uint i = inc_k + 1; i < ln; i++) sign[i] = -1;


    for (int i = (int) ln - 1; i >= (int) inc_k; i--) {

        old_k = k_vec[i];

        // node elimination
        if (i == inc_k and inc_k == ln - 1 && sign[i] > k[i]) {
            k[i] = sign[i];
            k_vec[i] = sign[i];
        }

        has_child = BuildBranch(i, r);

        if (has_child) {
            child = mlc_tree->GetNewTreeNode(k_vec, i);
            mlc_tree->SetRelChd(r, i, child);

            BuildSubMLCTreeNE(i, child, sign);
        }

        // restoreNaive
        k_vec[i] = old_k;
        Restore(old_e);
    }

    for (uint j = 0; j <= inc_k; j++) {
        f_sign[j] = std::min(sign[j], f_sign[j]);
    }
}

void MLCTreeBuilder::BuildSubMLCTreeSE(uint inc_k, Node *r, uint *f_sign) {
    bool has_child;
    uint old_e, old_k, *k, sign[ln];
    Node *child;

    old_e = k_core_index.e;
    k = MLCTree::GetK(r);

    for (uint i = 0; i <= inc_k; i++) {
        sign[i] = k_value_heaps[i].GetMinValue();
    }
    for (uint i = inc_k + 1; i < ln; i++) sign[i] = -1;


    for (int i = (int) ln - 1; i >= (int) inc_k; i--) {

        old_k = k_vec[i];

        // subtree elimination
        if (i == inc_k && sign[i] > k[i]) {
            k[i] = sign[i];
            k_vec[i] = sign[i];
        }

        has_child = BuildBranch(i, r);

        if (has_child) {
            child = mlc_tree->GetNewTreeNode(k_vec, i);
            mlc_tree->SetRelChd(r, i, child);

            BuildSubMLCTreeSE(i, child, sign);
        }

        // restoreNaive
        k_vec[i] = old_k;
        Restore(old_e);
    }

    for (uint j = 0; j <= inc_k; j++) {
        f_sign[j] = std::min(sign[j], f_sign[j]);
    }
}

void MLCTreeBuilder::BuildSubMLCTreeSM(uint inc_k, Node *r, Node *anc, Node *pr, uint *f_sign,
                                       unordered_map<Node *, uint *> &signs) {
    bool has_child, is_sub_eli;
    uint old_e, old_k, sign[ln];
    Node *cpr, *child;

    old_e = k_core_index.e;

    for (uint i = 0; i <= inc_k; i++) {
        sign[i] = k_value_heaps[i].GetMinValue();
    }
    for (uint i = inc_k + 1; i < ln; i++) sign[i] = -1;


    for (int i = (int) ln - 1; i >= (int) inc_k; i--) {

        old_k = k_vec[i];

        has_child = BuildBranch(i, r);

        if (has_child) {

            // subtree merge
            if (i == inc_k) cpr = pr;
            else if (anc) cpr = mlc_tree->GetRelChd(anc, i);
            else cpr = nullptr;

            while (cpr && MLCTree::GetK(cpr)[i] < k_vec[i]) {
                cpr = mlc_tree->GetRelChd(cpr, i);
            }
            is_sub_eli = cpr && ArrGe(signs[cpr], k_vec, i);

            if (is_sub_eli) {

                child = cpr;
                mlc_tree->SetRelChd(r, i, child);

                for (uint j = 0; j <= inc_k; j++) {
                    sign[j] = std::min(sign[j], signs[cpr][j]);
                }

            } else {
                child = mlc_tree->GetNewTreeNode(k_vec, i);
                mlc_tree->SetRelChd(r, i, child);

                BuildSubMLCTreeSM(i, child, r, cpr, sign, signs);
            }
        }

        k_vec[i] = old_k;
        Restore(old_e);
    }

    uint *arr = new uint[inc_k + 1];
    memcpy(arr, sign, (inc_k + 1) * sizeof(uint));
    signs.emplace(r, arr);

    for (uint j = 0; j <= inc_k; j++) {
        f_sign[j] = std::min(sign[j], f_sign[j]);
    }
}


void MLCTreeBuilder::BuildSubMLCTreeSESM(uint inc_k, Node *r, Node *anc, Node *pr, uint *f_sign,
                                         unordered_map<Node *, uint *> &signs) {
    bool has_child, is_sub_eli;
    uint old_e, old_k, *k, sign[ln];
    Node *cpr, *child;

    old_e = k_core_index.e;
    k = MLCTree::GetK(r);

    for (uint i = 0; i <= inc_k; i++) {
        sign[i] = k_value_heaps[i].GetMinValue();
    }
    for (uint i = inc_k + 1; i < ln; i++) sign[i] = -1;


    for (int i = (int) ln - 1; i >= (int) inc_k; i--) {

        // record the ith component of the current p_vec
        old_k = k_vec[i];

        // subtree elimination
        if (i == inc_k && sign[i] > k[i]) {
            k[i] = sign[i];
            k_vec[i] = sign[i];
        }

        has_child = BuildBranch(i, r);

        if (has_child) {

            // subtree merge
            if (i == inc_k) cpr = pr;
            else if (anc) cpr = mlc_tree->GetRelChd(anc, i);
            else cpr = nullptr;

            while (cpr && MLCTree::GetK(cpr)[i] < k_vec[i]) {
                cpr = mlc_tree->GetRelChd(cpr, i);
            }
            is_sub_eli = cpr && ArrGe(signs[cpr], k_vec, i);

            if (is_sub_eli) {

                child = cpr;
                mlc_tree->SetRelChd(r, i, child);

                for (uint j = 0; j <= inc_k; j++) {
                    sign[j] = std::min(sign[j], signs[cpr][j]);
                }

            } else {
                child = mlc_tree->GetNewTreeNode(k_vec, i);
                mlc_tree->SetRelChd(r, i, child);

                BuildSubMLCTreeSESM(i, child, r, cpr, sign, signs);
            }
        }

        // restore
        k_vec[i] = old_k;
        Restore(old_e);
    }

    uint *arr = new uint[inc_k + 1];
    memcpy(arr, sign, (inc_k + 1) * sizeof(uint));
    signs.emplace(r, arr);

    for (uint j = 0; j <= inc_k; j++) {
        f_sign[j] = std::min(sign[j], f_sign[j]);
    }
}

bool MLCTreeBuilder::BuildBranch(uint i, Node *r) {
    uint old_e, mlc_diff_size, num_of_removed_vtx;
    IntLinearHeap *k_value_heap;

    k_value_heap = &k_value_heaps[i];
    num_of_removed_vtx = k_value_heap->GetElements(k_vec[i], aux_arr[0]);

    k_vec[i] += 1;
    old_e = k_core_index.e;

    if (num_of_removed_vtx) {
        for (uint j = 0; j < num_of_removed_vtx; j++) {
            k_core_index.Remove(aux_arr[0][j]);
        }
        Peel();
    }

    mlc_diff_size = k_core_index.e - old_e;
    if (i == ln - 1 && mlc_diff_size) {
        mlc_tree->SetDiff(r, mlc_diff_size, k_core_index.vert + old_e);
    }

    return k_core_index.e < n;
}
