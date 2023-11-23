//
// Created by ldd on 2023/6/10.
//

#ifndef MLCDEC_METRIC_H
#define MLCDEC_METRIC_H

#include "../Graphs/MultilayerGraph.h"
#include "../Utils/ArrayUtils.h"

class Metric {

public:
    static float GetDensity(Graph &g, const uint *vtx, uint length, const bool *sign) {
        uint n_edges = 0, v, u, **adj_lst;

        adj_lst = g.GetAdjLst();

        for (uint i = 0; i < length; i++) {
            v = vtx[i];
            for (uint j = 1; j <= adj_lst[v][0]; j++) {
                u = adj_lst[v][j];
                if (sign[u] && v > u) {
                    n_edges += 1;
                }
            }
        }

        return (float) n_edges / (float) length;
    }

    static float
    MinLayerDensity(MultilayerGraph &mg, const float *w, const uint *vtx, uint length, const uint *is_count) {

        uint n = mg.GetN(), ln = mg.GetLayerNumber();

        bool sign[n];
        float min_ww = DBL_MAX, ww;

        memset(sign, false, n * sizeof(bool));
        SetSign(vtx, length, sign);

        for (uint i = 0; i < ln; i++) {
            if (is_count[i]) {
                ww = GetDensity(mg.GetGraph(i), vtx, length, sign) * w[i];
                if (ww < min_ww) min_ww = ww;
            }
        }

        return min_ww;
    }

    static float
    MLDensity(MultilayerGraph &mg, const float *w, float beta, const uint *vtx, uint length) {
        uint n = mg.GetN(), ln = mg.GetLayerNumber();

        bool sign[n];
        float max_ww  = 0, ww[ln], den;

        memset(sign, false, n * sizeof(bool));
        SetSign(vtx, length, sign);

        for (uint i = 0; i < ln; i++) {
            ww[i] = GetDensity(mg.GetGraph(i), vtx, length, sign) * w[i];
        }
        std::sort(&ww[0], &ww[ln]);
        for (int i = (int) ln - 1; i >= 0; i--) {
            den = ww[i] * (float) pow(ln - i, beta);
            if (den > max_ww) {
                max_ww = den;
            }
        }
        return max_ww;
    }

    static void Density(MultilayerGraph &mg, uint *arr, uint length, vector<float> &res) {
        uint n = mg.GetN(), ln = mg.GetLayerNumber();
        bool sign[n];

        memset(sign, false, n * sizeof(bool));
        SetSign(arr, length, sign);

        for (uint i = 0; i < ln; i++) {
            res.emplace_back(GetDensity(mg.GetGraph(i), arr, length, sign));
        }
    }

    static void LocalCC(Graph &g, const uint *arr, uint length, vector<float> &res) {
        uint v, u, n_common_edge, d, **adj_lst = g.GetAdjLst();

//#pragma omp parallel for
        for (uint i = 0; i < length; i++) {

            v = arr[i];
            d = adj_lst[v][0];
            n_common_edge = 0;

            for (uint j = 1; j <= d; j++) {
                u = adj_lst[v][j];
                n_common_edge += CountIntersect(adj_lst[u] + 1, adj_lst[u][0], adj_lst[v] + 1, d);
            }

            res.emplace_back(((float) n_common_edge) / (float) (d * (d - 1)));
        }
    }


    static void LocalCC(MultilayerGraph &mg, const uint *arr, uint length, vector<vector<float>> &res) {
        uint n = mg.GetN(), ln = mg.GetLayerNumber(), pos[n], arr_remap[length];

        memset(pos, -1, n * sizeof(uint));
        for (uint i = 0; i < length; i++) {
            pos[arr[i]] = i;
            arr_remap[i] = i;
        }

        res.resize(ln);
        for (uint i = 0; i < ln; i++) {
            Graph sg;
            BuildSg(mg.GetGraph(i), arr, length, pos, sg);
            LocalCC(sg, arr_remap, length, res[i]);
        }
    }

private:

    static void BuildSg(Graph &g, const uint *arr, uint length, const uint *pos, Graph &sg) {
        uint v, u, num_of_edge, edge_buf_size, **adj_lst = g.GetAdjLst();
        edge *edge_buf;

        num_of_edge = 0;
        edge_buf_size = DEFAULT_EDGE_BUF_SIZE;
        edge_buf = new edge[edge_buf_size];

        for (uint i = 0; i < length; i++) {
            v = arr[i];

            for (uint j = 1; j <= adj_lst[v][0]; j++) {

                u = adj_lst[v][j];
                if (pos[u] != -1) {

                    if (num_of_edge + 1 > edge_buf_size) {
                        edge_buf_size = edge_buf_size << 1;
                        auto tmp_edge_buf = new edge[edge_buf_size];

                        memcpy(tmp_edge_buf, edge_buf, num_of_edge * sizeof(edge));
                        delete[] edge_buf;
                        edge_buf = tmp_edge_buf;
                    }

                    edge_buf[num_of_edge++] = edge(i, pos[u]);
                }
            }
        }

        sg.BuildFromEdgeLst(edge_buf, length, num_of_edge);
        delete edge_buf;
    }
};

#endif //MLCDEC_METRIC_H
