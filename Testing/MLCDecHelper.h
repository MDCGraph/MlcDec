//
// Created by ldd on 2023/6/12.
//

#ifndef MLCDEC_MLCDECHELPER_H
#define MLCDEC_MLCDECHELPER_H

#include "../Core/CoreCube.h"
#include "../Core/MLCTreeBuilder.h"
#include "../Core/MLCTreeBuilder_simp.h"
#include "../Core/PMLCTreeBuilder.h"
#include "../Core/MLCTreeBuildOption.h"

#include "../Utils/MemoryUtils.h"
#include "../DS/DSMLCTreeBuilder.h"

static void
RunMLCTreeBuilder(MultilayerGraph &mg, const string &output, bool flush) {

    auto suffix = "_mlct_" + RunOpt2Str(SERIAL);

    double runtime;
    long double mem;
    ll_uint n_nodes;

    Timer timer;

    MLCTree mlc_t(mg.GetLayerNumber(), mg.GetN());
    MLCTreeBuilder_simp::Execute(mg, mlc_t);
    timer.Stop();

    n_nodes = mlc_t.GetNumOfNodes();
    runtime = timer.GetTimeInSec();
    mem = GetPeakRSSInMB();

    cout << RunOpt2Str(SERIAL) << "_runtime = " << runtime << "s, #nodes = " << n_nodes << ", mem = " << mem
         << "MB" << endl;

    /* write to file */
    auto out = ofstream(output + suffix + ".txt");
    out << "runtime = " << runtime << "s, #nodes = " << n_nodes << ", mem = " << mem
         << "MB" << endl;
    out.close();

    /* save mlc-tree */
    if (flush) mlc_t.Flush(output + suffix);
}


static void
RunDSMLCTreeBuilder(MultilayerGraph &mg, const string &output, bool flush) {

    auto suffix = "_ds_mlct_" + RunOpt2Str(SERIAL);

    double runtime;
    long double mem;
    ll_uint n_nodes;

    Timer timer;

    DSMLCTree mlc_t(mg.GetLayerNumber(), mg.GetN());
    DSMLCTreeBuilder::Execute(mg, mlc_t);
    timer.Stop();

    n_nodes = mlc_t.GetNumOfNodes();
    runtime = timer.GetTimeInSec();
    mem = GetPeakRSSInMB();

    cout << RunOpt2Str(SERIAL) << "_runtime = " << runtime << "s, #nodes = " << n_nodes << ", mem = " << mem
    << "MB" << endl;

    /* write to file */
    auto out = ofstream(output + suffix + ".txt");
    out << "runtime = " << runtime << "s, #nodes = " << n_nodes << ", mem = " << mem
    << "MB" << endl;
    out.close();

    /* save mlc-tree */
    if (flush) mlc_t.Flush(output + suffix);
}


static void
RunMLCTreeBuilder(MultilayerGraph &mg, const string &output, bool flush, RUN_OPTION opt, uint nt, uint l, float alpha) {

    auto suffix = "_mlct_" + RunOpt2Str(opt) + "_" + std::to_string(nt);

    double runtime;
    long double mem;
    ll_uint n_nodes;

    omp_set_num_threads((int) nt);

    Timer timer;
    MLCTree mlc_t(mg.GetLayerNumber(), mg.GetN());
    if (opt == P_OPT) PMLCTreeBuilder::PRun_async_with_merge(mg, mlc_t, l, alpha);
    else if (opt == IP_START) PMLCTreeBuilder::PRun_async(mg, mlc_t, l);
    else if (opt == IP_MERGE) PMLCTreeBuilder::InterPathPRun_with_merge(mg, mlc_t, alpha);
    else if (opt == IP) PMLCTreeBuilder::InterPathPRun(mg, mlc_t);
    timer.Stop();

    n_nodes = mlc_t.GetNumOfNodes();
    runtime = timer.GetTimeInSec();
    mem = GetPeakRSSInMB();

    cout << "#threads = " << nt << ", level = " << l << ", alpha = " << alpha << ", " << RunOpt2Str(opt)
         << " runtime = " << timer.GetTimeInSec() << "s, #nodes = " << mlc_t.GetNumOfNodes() << ", mem = " << mem
         << "MB" << endl;

    /* write to file */
    auto out = ofstream(output + suffix + "_" + std::to_string(l) + "_" + to_string(alpha) + ".txt");
    out << "runtime = " << runtime << "s, #nodes = " << n_nodes << ", mem = " << mem
        << "MB" << endl;
    out.close();

    /* save mlc-tree */
    if (flush) mlc_t.Flush(output + suffix);
}

static void RunMLCHashTableBuilder(MultilayerGraph &mg, uint bucket_size, const string &output, bool flush = true) {
    auto suffix = "_mlcht";
    auto out = ofstream(output + suffix + ".txt");

    double runtime;
    long double mem;
    ll_uint n_nodes, n_distinct_cores;

    Timer timer;

    MLCHashTable mlc_ht(mg.GetN(), mg.GetLayerNumber(), bucket_size);
    MLCTreeBuilder_simp::Execute(mg, mlc_ht);
    timer.Stop();

    n_nodes = mlc_ht.GetNumOfMLC();
    runtime = timer.GetTimeInSec();
    mem = GetPeakRSSInMB();

    n_distinct_cores = mlc_ht.CountDistinct();

    cout << "ht runtime = " << runtime << "s, #nodes = " << n_nodes << ", #distinct_cores = " << n_distinct_cores
         << ", mem = " << mem << "MB" << endl;

    out << "runtime = " << runtime << "s, #nodes = " << n_nodes << ", #distinct_cores = " << n_distinct_cores
        << ", mem = " << mem << "MB" << endl;
    out.close();

    if (flush) {
        mlc_ht.Flush(output + suffix);
    }
}


static void RunCorecubeGenerator(MultilayerGraph &mg, const string &output, bool flush) {
    auto suffix = "_corecube";
    auto out = ofstream(output + suffix + ".txt");
    double dec_time, compact_time;
    long double mem;

    Timer timer;

    Cube cube{};
    CoreCube cube_constructor(mg);
    cube_constructor.TD(cube);
    timer.Stop();
    dec_time = timer.GetTimeInSec();

    HYB_Cube hyb_cube{};
    cube_constructor.GenHybridStorage(hyb_cube);
    timer.Stop();
    compact_time = timer.GetTimeInSec() - dec_time;

    mem = GetPeakRSSInMB();

    cout << "dec_time = " << dec_time << "s, compact_time = " << compact_time << "s, n_comb = " << hyb_cube.GetNComb()
         << ", mem = " << mem << "MB" << endl;
    out << "dec_time = " << dec_time << "s, compact_time = " << compact_time << "s, n_comb = " << hyb_cube.GetNComb()
        << ", mem = " << mem << "MB" << endl;

    out.close();

    if (flush) {
        hyb_cube.Flush(output + suffix);
    }
}

#endif //MLCDEC_MLCDECHELPER_H
