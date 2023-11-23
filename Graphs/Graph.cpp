//
// Created by ldd on 2023/3/6.
//

#include "Graph.h"


Graph::~Graph() {
    delete[] adj_lst;
    delete[] adj_lst_buf;
}

// Using CSR to store graph.
void Graph::BuildFromEdgeLst(edge *edge_buf, uint num_of_vtx, uint num_of_edge) {
    uint i, j, pj;

    n = num_of_vtx;
    adj_lst = new uint *[n];
    adj_lst_buf = new uint [num_of_edge + n];  // +1

    sort(edge_buf, edge_buf + num_of_edge);

    i = 0; // index of edge_buf
    j = 0; // index of adj_lst_buf
    for (uint v = 0; v < n; v++) {
        adj_lst[v] = &adj_lst_buf[j];
        if (edge_buf[i].s > v || i >= num_of_edge) adj_lst_buf[j++] = 0;
        else {
            pj = j++; // index to store vtx degree
            adj_lst_buf[j++] = edge_buf[i++].t;
            while (i < num_of_edge && edge_buf[i].s == v) {
                if (edge_buf[i].t == edge_buf[i - 1].t) i++;  // duplicated edge
                else adj_lst_buf[j++] = edge_buf[i++].t;
            }
            adj_lst_buf[pj] = j - pj - 1;
            m += adj_lst_buf[pj];
            if (adj_lst_buf[pj] > max_deg) max_deg = adj_lst_buf[pj];
        }
    }
    //adj_lst_buf[j] = 0; // ended with 0, representing the intra-deg for all (newly added) vertices with no neighbors.
}


void Graph::BuildFromGraph(uint n_, uint m_, uint **adj_lst_, uint *adj_lst_buf_) {
        n = n_;
        m = m_;
        adj_lst = adj_lst_;
        adj_lst_buf = adj_lst_buf_;
}
