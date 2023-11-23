//
// Created by ldd on 2023/6/10.
//

#ifndef MLCDEC_WDSHELPER_H
#define MLCDEC_WDSHELPER_H

#include "../Core/PMLCTreeBuilder.h"
#include "../DS/PWeightedDS.h"
#include "../Utils/StringUtils.h"
#include "../Core/MLCTreeBuilder_simp.h"
#include "../Core/MLCore_simp.h"
#include "../Utils/Timer.h"

void
SearchWDS(MultilayerGraph &mg, float beta, vector<float> &w, const string &output, uint level, float merge_factor, uint nt) {

    string suffix = "_" + to_string(nt) + "_" + Vec2Str(w) + "_" + std::to_string(beta) + "_wds.txt";

    uint n = mg.GetN(), ln = mg.GetLayerNumber(), ds[n], length;
    float density;
    vector<uint> k;

    vector<float> den;
    vector<vector<float>> lcc;

    omp_set_num_threads((int) nt);
    Timer timer;
    density = PWeightedDS::Search(mg, beta, w, ds, length, k, level, merge_factor);
    timer.Stop();

    Metric::Density(mg, ds, length, den);

    cout << "n_threads = " << nt << ", runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
         << ", layer_density = " << Vec2Str(den) << endl;

    Metric::LocalCC(mg, ds, length, lcc);
    auto out = ofstream(output + suffix);
    out << Arr2Str(ds, length) << endl;
    out << "runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
        << ", layer_density = " << Vec2Str(den) << endl;
    for (uint i = 0; i < ln; i++) {
        out << "local_cc_l" << i << " = " << Vec2Str(lcc[i]) << endl;
    }
    out.close();
}

void SearchWDS(MultilayerGraph &mg, MLCTree &mlc_t, float beta, vector<float> &w, const string &output, uint nt) {

    string suffix = "_" + to_string(nt) + "_" + Vec2Str(w) + "_" + std::to_string(beta) + "_mlct_wds.txt";

    uint n = mg.GetN(), ln = mg.GetLayerNumber(), ds[n], length;
    float density;
    vector<uint> k;

    vector<float> den;
    vector<vector<float>> lcc;

    omp_set_num_threads((int) nt);
    Timer timer;
    density = PWeightedDS::Search(mg, beta, w, mlc_t, ds, length, k);
    timer.Stop();

    Metric::Density(mg, ds, length, den);

    cout << "n_threads = " << nt << ", runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
         << ", layer_density = " << Vec2Str(den) << endl;

    Metric::LocalCC(mg, ds, length, lcc);
    auto out = ofstream(output + suffix);
    out << Arr2Str(ds, length) << endl;
    out << "runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
        << ", layer_density = " << Vec2Str(den) << endl;
    for (uint i = 0; i < ln; i++) {
        out << "local_cc_l" << i << " = " << Vec2Str(lcc[i]) << endl;
    }
    out.close();
}

void SearchWDS(MultilayerGraph &mg, DSMLCTree &mlc_t, float beta, vector<float> &w, const string &output, uint nt) {

    string suffix = "_" + to_string(nt) + "_" + Vec2Str(w) + "_" + std::to_string(beta) + "_ds_mlct_wds.txt";

    uint n = mg.GetN(), ln = mg.GetLayerNumber(), ds[n], length;
    float density;
    vector<uint> k;

    vector<float> den;
    vector<vector<float>> lcc;

    omp_set_num_threads((int) nt);
    Timer timer;
    density = PWeightedDS::Search(beta, w, mlc_t, ds, length, k);
    timer.Stop();

    Metric::Density(mg, ds, length, den);

    cout << "n_threads = " << nt << ", runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
    << ", layer_density = " << Vec2Str(den) << endl;

    Metric::LocalCC(mg, ds, length, lcc);
    auto out = ofstream(output + suffix);
    out << Arr2Str(ds, length) << endl;
    out << "runtime = " << timer.GetTimeInSec() << "s, density = " << density << ", k = " << Vec2Str(k)
    << ", layer_density = " << Vec2Str(den) << endl;
    for (uint i = 0; i < ln; i++) {
        out << "local_cc_l" << i << " = " << Vec2Str(lcc[i]) << endl;
    }
    out.close();
}


void ComputeCohesiveness(MultilayerGraph &mg, const string &csg_file, const string &output) {

    string line;
    vector<ll_uint> csg_vtx;
    unordered_map<ll_uint, uint> vtx2id;
    uint length, ln = mg.GetLayerNumber();

    getline(ifstream(csg_file), line);
    Str2LLUVec(line, csg_vtx);

    mg.LoadVtx2IdMap(vtx2id);
    length = csg_vtx.size();
    uint ds[length];

    for (uint i = 0; i < length; i++) {
        ds[i] = vtx2id[csg_vtx[i]];
    }

    vector<float> den;
    vector<vector<float>> lcc;

    Metric::Density(mg, ds, length, den);
    Metric::LocalCC(mg, ds, length, lcc);

    auto out = ofstream(output + GetFileName(csg_file) + "_wds.txt");
    out << "layer_density = " << Vec2Str(den) << endl;
    for (uint i = 0; i < ln; i++) {
        out << "local_cc_l" << i << " = " << Vec2Str(lcc[i]) << endl;
    }
    out.close();
}

#endif //MLCDEC_WDSHELPER_H
