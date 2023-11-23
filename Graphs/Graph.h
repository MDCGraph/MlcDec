//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_GRAPH_H
#define MLCDEC_GRAPH_H


#include "../Header.h"


struct edge {
    uint s{0};
    uint t{0};

    edge() = default;

    edge(uint s_, uint t_) : s(s_), t(t_) {};

    bool operator<(const edge &e) const {
        return s < e.s || (s == e.s && t < e.t);
    }
};

class Graph {
public:
    Graph() = default;
    ~Graph();

    void BuildFromEdgeLst(edge *edge_buf, uint num_of_vtx, uint num_of_edge);
    void BuildFromGraph(uint n_, uint m_, uint **adj_lst_, uint *adj_lst_buf_);

    [[nodiscard]] uint GetN() const {
        return n;
    }

    [[nodiscard]] uint GetM() const {
        return m;
    }

    [[nodiscard]] uint GetMaxDeg() const {
        return max_deg;
    }

    [[nodiscard]] uint **GetAdjLst() {
        return adj_lst;
    }

    void SetNull() {
        adj_lst_buf = nullptr;
        adj_lst = nullptr;
    }

protected:
    uint *adj_lst_buf{nullptr};
    uint **adj_lst{nullptr};

    uint m{0};
    uint n{0};

    uint max_deg{0};
};


#endif //MLCDEC_GRAPH_H
