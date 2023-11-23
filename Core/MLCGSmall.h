//
// Created by ldd on 2023/10/24.
//

#ifndef MLCDEC_MLCGSMALL_H
#define MLCDEC_MLCGSMALL_H


#include "../Graphs/Graph.h"

class G_small {
public:
    G_small(uint n_, uint ln_) : n(n_), n_layers(ln_) {
        graph_layers = new Graph[n_layers];
    }

    ~G_small() {
        for (uint i = 0; i < n_layers; i++) graph_layers[i].SetNull();
        delete[] graph_layers;
    }

    void BuildSmallG(MultilayerGraph&mg, uint *relabel_, const uint *vtx2id, uint **adj_buf, uint *nbr_buf) {
        uint vtx, nbr, **adj_lst, **new_adj_lst, off = 0;
        relabel = relabel_;
        new_adj_lst = adj_buf;

        for (uint i = 0; i < n_layers; i++) {
            adj_lst = mg.GetGraph(i).GetAdjLst();

            for (uint j = 0; j < n; j++) {
                vtx = relabel[j];

                new_adj_lst[j] = nbr_buf + off;
                new_adj_lst[j][0] = 0;

                for (uint k = 1; k <= adj_lst[vtx][0]; k++) {
                    nbr = adj_lst[vtx][k];
                    if (vtx2id[nbr] != -1) {
                        new_adj_lst[j][++ new_adj_lst[j][0]] = vtx2id[nbr];
                    }
                }

                off += new_adj_lst[j][0] + 1;

            }

            graph_layers[i].BuildFromGraph(n , off - n, new_adj_lst, new_adj_lst[0]);
            new_adj_lst += n;
        }


    }

    [[nodiscard]] inline uint GetLayerNumber() const {
        return n_layers;
    }

    [[nodiscard]] inline uint GetN() const {
        return n;
    }

    [[nodiscard]] inline uint* GetRelabel() const {
        return relabel;
    }

    inline Graph &GetGraph(uint i) {
        return graph_layers[i];
    }

    inline Graph *GetGraphAddr(uint i) {
        return &graph_layers[i];
    }

private:
    Graph *graph_layers{nullptr};

    uint n_layers{0};
    uint n{0};

    uint *relabel{nullptr};
};


#endif //MLCDEC_MLCGSMALL_H
