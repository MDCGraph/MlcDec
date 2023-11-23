//
// Created by ldd on 2023/6/12.
//

#ifndef MLCDEC_MLCSEARCHHELPER_H
#define MLCDEC_MLCSEARCHHELPER_H

#include "../Utils/RandVecGenerator.h"
#include "../Core/MLCHashTable.h"
#include "../Core/MLCTree.h"
#include "../Core/MLCore_simp.h"


static void MLCSearch(MultilayerGraph &mg, vector<vector<uint>> &ks, const string &output) {

    uint n = mg.GetN(), core[n], length;

    auto out = ofstream(output + "_mlcs.txt");

    Timer timer;

    for (const auto &k : ks) {
        length = MLCore_simp::Extract(mg, k, core);
    }
    timer.Stop();

    cout << "#test_cases = " << ks.size() << ", mlcs_runtime = " << timer.GetTimeInSec() << "s" << endl;

    out << "#test_cases = " << ks.size() << ", mlcs_runtime = " << timer.GetTimeInSec() << "s" << endl;
    out.close();
}

static void MLCSearch(MLCTree &mlc_t, vector<vector<uint>> &ks, const string &output, uint nt) {

    uint n = mlc_t.GetN(), core[n], length, n_visits = 0;

    auto out = ofstream(output + "_mlcs_mlct_" + to_string(nt) + ".txt");

    Timer timer;
    if (nt == -1 || nt == 0) {
        for (const auto &k : ks) {
            n_visits += mlc_t.Search(k, core, length);
        }
    } else {
        omp_set_num_threads((int) nt);
        for (const auto &k : ks) {
            n_visits += mlc_t.PSearch(k, core, length);
        }
    }

    timer.Stop();

    cout << "#test_cases = " << ks.size() << ", #threads = " << nt << ", mlcs_mlct_runtime = " << timer.GetTimeInSec() << "s, #visited_nodes = "
         << n_visits << endl;

    /* save result to file */
    out << "#test_cases = " << ks.size() << ", mlcs_mlct_runtime = " << timer.GetTimeInSec() << "s, #visited_nodes = "
        << n_visits << endl;
    out.close();
}

static void MLCSearch(MLCHashTable &mlc_ht, vector<vector<uint>> &ks, const string &output) {

    uint n = mlc_ht.GetN(), core[n], length;

    auto out = ofstream(output + "_mlcs_ht.txt");

    Timer timer;
    for (auto &k : ks) {
        length = mlc_ht.Search(k, core);
    }
    timer.Stop();

    cout << "#test_cases = " << ks.size() << ", mlcs_ht_runtime = " << timer.GetTimeInSec() << "s";

    /* save result to file */
    out << "#test_cases = " << ks.size() << ", mlcs_ht_runtime = " << timer.GetTimeInSec() << "s";
    out.close();
}


#endif //MLCDEC_MLCSEARCHHELPER_H
