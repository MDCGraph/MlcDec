//
// Created by ldd on 2023/6/18.
//

#ifndef MLCDEC_GRAPHHELPER_H
#define MLCDEC_GRAPHHELPER_H

#include "../Graphs/MultilayerGraph.h"
#include "../DS/Metric.h"

static void GetStatistics(MultilayerGraph &mg, const string &output) {
    uint ln = mg.GetLayerNumber(), n = mg.GetN();
    uint m[ln], cn[ln], min_m = INT32_MAX, max_m = 0, sum_m = 0;
    float density[ln], avg_m;

    for (uint i = 0; i < ln; i++) {
        m[i] = mg.GetGraph(i).GetM() >> 1;
        density[i] = (float) m[i] / (float) n;
    }

    for (uint i = 0; i < ln; i++) {
        cn[i] = KC::Degeneracy(mg.GetGraph(i));
    }

    for (uint i = 0; i < ln; i++) {
        if (m[i] < min_m) min_m = m[i];
        if (m[i] > max_m) max_m = m[i];
        sum_m += m[i];
    }

    avg_m = (float) sum_m / (float) ln;

    auto f = ofstream(output + "_info.txt");
    f << "|L| = " << ln << endl;
    f << "|V| = " << n << endl;
    f << "|E| = " << sum_m << endl;
    f << "min_|E_i| = " << min_m << endl;
    f << "max_|E_i| = " << max_m << endl;
    f << "avg_|E_i| = " << avg_m << endl;
    f << "E_i = " << Arr2Str(m, ln) << endl;
    f << "coreness = " << Arr2Str(cn, ln) << endl;
    f << "density = " << Arr2Str(density, ln) << endl;
    f.close();
}

#endif //MLCDEC_GRAPHHELPER_H
