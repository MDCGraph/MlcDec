//
// Created by ldd on 2023/3/7.
//

#include "KC.h"

uint KC::CoreDec(Graph &g, uint *core_number) {
    uint n = g.GetN(), max_deg = g.GetMaxDeg(), start, num, v, pu, pw, w, u, **adj_lst;
    uint pos[n], vert[n], bin[max_deg + 1];

    adj_lst = g.GetAdjLst();
    for (uint i = 0; i < n; i++) {
        core_number[i] = adj_lst[i][0];
    }
    memset(bin, 0, (max_deg + 1) * sizeof(uint));

    // Bin sort
    for (uint i = 0; i < n; i++) {
        bin[core_number[i]]++;
    }

    start = 0;
    for (uint d = 0; d <= max_deg; d++) {
        num = bin[d];
        bin[d] = start;
        start += num;
    }

    for (uint i = 0; i < n; i++) {
        pos[i] = bin[core_number[i]];
        vert[pos[i]] = i;
        bin[core_number[i]]++;
    }

    for (uint d = max_deg; d >= 1; d--) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;

    // Compute core numbers.
    for (uint i = 0; i < n; i++) {
        v = vert[i];
        for (uint j = 1; j <= adj_lst[v][0]; j++) {
            u = adj_lst[v][j];
            if (core_number[u] > core_number[v]) {
                pu = pos[u], pw = bin[core_number[u]], w = vert[pw];

                vert[pw] = u, vert[pu] = w;
                pos[u] = pw, pos[w] = pu;

                bin[core_number[u]]++;
                core_number[u]--;
            }
        }
    }

    // Return the maximum core number.
    return core_number[v];
}

uint KC::Degeneracy(Graph &g) {
    uint cn[g.GetN()];
    return CoreDec(g, cn);
}