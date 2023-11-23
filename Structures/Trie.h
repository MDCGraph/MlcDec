//
// Created by ldd on 2023/3/7.
//

#ifndef MLCDEC_TRIE_H
#define MLCDEC_TRIE_H


#include "DataBuf.h"

/*
 * Below implements a trie structure used for de-duplicating multi-layer cores,
 * referencing the source code of the QUICK program developed by Guimei Liu released with the paper:
 *
 * "Liu G, Wong L. Effective pruning techniques for mining quasi-cliques[J].
 *  ECML/PKDD (2), 2008, 5212: 33-49."
 *
 */


struct TrieNode {
    uint val{static_cast<uint>(-1)};
    TrieNode *l_chd{nullptr};
    TrieNode *r_sib{nullptr};
    bool is_end{false};
};


class Trie {
public:
    Trie() : n_dis(0), root(nullptr) {}

    ~Trie() {
        trie_buf.Release();
    }

    void Add(const uint *s, uint length) {
        uint i, j;
        TrieNode *node, *new_node, *parent{nullptr}, *left_sib{nullptr};

        i = 0;
        node = root;

        while (i < length) {

            while (node && node->val < s[i]) {
                left_sib = node;
                node = node->r_sib;
            }
            if (!node || node->val > s[i]) {
                new_node = trie_buf.Allocate();
                new_node->val = s[i];
                new_node->r_sib = node;

                if (left_sib) left_sib->r_sib = new_node;
                else if (parent) parent->l_chd = new_node;

                if (i == 0 && !left_sib) root = new_node;
                parent = new_node;

                for (j = i + 1; j < length; j++) {
                    new_node = trie_buf.Allocate();
                    new_node->val = s[j];
                    parent->l_chd = new_node;
                    parent = new_node;
                }

                parent->is_end = true;
                n_dis += 1;

                break;

            } else {
                parent = node;
                node = node->l_chd;
                left_sib = nullptr;
                i++;
            }
        }

        if (parent && !parent->is_end) {
            parent->is_end = true;
            n_dis += 1;
        }
    }

    [[nodiscard]] ll_uint GetCount() const {
        return n_dis;
    }

private:
    MDataBuf<TrieNode> trie_buf;
    TrieNode *root;

    ll_uint n_dis;
};

#endif //MLCDEC_TRIE_H
