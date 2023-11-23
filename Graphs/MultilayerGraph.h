//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_MULTILAYERGRAPH_H
#define MLCDEC_MULTILAYERGRAPH_H

#include "../Core/KC.h"
#include "GOrder.h"
#include "Graph.h"

const uint DEFAULT_EDGE_BUF_SIZE = 5000;

class MultilayerGraph {
public:
    MultilayerGraph() = default;
    ~MultilayerGraph();

    void LoadFromFile(const string &input_path);
    void SetGraphOrder(G_ORDER ordering);
    void PrintStatistics();

    void LoadId2VtxMap(ll_uint * id2vtx);
    void LoadVtx2IdMap(unordered_map<ll_uint, uint> &vtx2id);

    static MultilayerGraph* Load(const string &file);

    [[nodiscard]] inline uint GetLayerNumber() const {
        return n_layers;
    }

    [[nodiscard]] inline uint GetN() const {
        return n;
    }

    inline Graph &GetGraph(uint i) {
        return graph_layers[order[i]];
    }

    // for test only
    [[nodiscard]] inline uint* GetOrder() const {
        return order;
    }

private:
    Graph *graph_layers{nullptr};
    uint *order{nullptr};

    uint n_layers{0};
    uint n{0};

    string map_file;

    static void GetGraphFile(const string &graph_path, vector<string> &graph_files);
    static uint LoadLayer(const string &graph_file, edge *&edge_buf, unordered_map<ll_uint, uint> &vtx2id,
                         std::basic_ofstream<char> &map_file_out);
};


#endif //MLCDEC_MULTILAYERGRAPH_H
