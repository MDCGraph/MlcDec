#include "Testing/WDSHelper.h"
#include "Testing/MLCSearchHelper.h"
#include "Testing/MLCDecHelper.h"
#include "Testing/GraphHelper.h"
#include "CmdParser.h"

void set_limit(uint size) {
    rlimit limit{};
    limit.rlim_cur = size * 1024 * 1024;
    limit.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_STACK, &limit)) {
        printf("Set limit failed\n");
        exit(-1);
    }
}


int main(int argc, char **argv) {

#if USE_JEMALLOC
    cout << "use jemalloc" << endl;
#endif

    set_limit(500);  // 500MB

    Params params;
    ParseCmd(params, argc, argv);

    assert_non_empty_graph_name(params);
    MultilayerGraph mg;

    if (params.path.empty()) {
        assert_default_graph_name(params);
        mg.LoadFromFile(params.default_path + params.dataset + '/');
    } else {
        assert_valid_graph_path(params);
        mg.LoadFromFile(params.path);
    }

    mg.PrintStatistics();
    mg.SetGraphOrder(params.graph_ordering);

    warn_empty_output_path(params);
    const string output_prefix = params.output + params.dataset + "_" + GraphOrder2Str(params.graph_ordering);

    for (auto mode : params.mode) {

        if (mode == PRT_GRAPH_INFO) {
            GetStatistics(mg, output_prefix);

        } else if (mode == RUN_MLC_SEARCH) {
            assert_valid_kvec_file(params);

            vector<vector<uint>> vec;
            LoadNNegIntVec(mg.GetLayerNumber(), params.kvf, vec);

            MLCSearch(mg, vec, output_prefix);

        } else if (mode == RUN_MLCT_MLC_SEARCH) {
            assert_valid_tree_file(params);
            assert_valid_kvec_file(params);

            vector<vector<uint>> vec;
            LoadNNegIntVec(mg.GetLayerNumber(), params.kvf, vec);

            MLCTree *mlct = MLCTree::Load(params.mlct_file);
            MLCSearch(*mlct, vec, output_prefix, params.nt);

        } else if (mode == RUN_HT_MLC_SEARCH) {
            assert_valid_ht_file(params);
            assert_valid_kvec_file(params);

            vector<vector<uint>> vec;
            LoadNNegIntVec(mg.GetLayerNumber(), params.kvf, vec);

            MLCHashTable *mlc_ht = MLCHashTable::Load(params.mlcht_file);
            MLCSearch(*mlc_ht, vec, output_prefix);

        } else if (mode == BUILD_MLC_TREE) {
            if (params.prun_opt == SERIAL) {
                RunMLCTreeBuilder(mg, output_prefix, params.flush);
            } else {
                verify_prun_param(params);
                if (params.nt == -1) {
                    for (auto nt:{1, 2, 4, 8, 16, 32, 40, 80}) {
                        RunMLCTreeBuilder(mg, output_prefix, params.flush, params.prun_opt, nt, params.l, params.alpha);
                    }
                } else {
                    RunMLCTreeBuilder(mg, output_prefix, params.flush, params.prun_opt, params.nt, params.l,
                                      params.alpha);
                }
            };

        } else if (mode == BUILD_DS_MLC_TREE) {
            RunDSMLCTreeBuilder(mg, output_prefix, params.flush);

        } else if (mode == BUILD_CORE_CUBE) {
            RunCorecubeGenerator(mg, output_prefix, params.flush);

        } else if (mode == BUILD_MLC_HT) {
            cout << "Bucket size is set to " << params.bucket_size << "." << endl;
            RunMLCHashTableBuilder(mg, params.bucket_size, output_prefix, params.flush);

        } else if (mode == WDS_SEARCH) {
            assert_valid_beta(params);
            assert_valid_w(params);

            if (!params.mlct_file.empty()) {

                MLCTree *mlct = MLCTree::Load(params.mlct_file);
                SearchWDS(mg, *mlct, params.beta, params.w, output_prefix, params.nt);
            } else {
                SearchWDS(mg, params.beta, params.w, output_prefix, params.l, params.alpha, params.nt);
            }

        } else if (mode == WDS_TREE_SEARCH) {
            assert_valid_beta(params);
            assert_valid_w(params);
            assert_valid_tree_file(params);

            DSMLCTree *mlct = DSMLCTree::Load(params.mlct_file);
            SearchWDS(mg, *mlct, params.beta, params.w, output_prefix, params.nt);

        } else if (mode == DCC_WDS_SEARCH) {
            assert_valid_beta(params);
            assert_valid_w(params);
            assert_valid_cc_file(params);

            /*
             * TODO:
             *
             *  HYB_Cube *hyb_cube = HYB_Cube::Load(params.cc_file);
             *  SearchWDS(mg, *hyb_cube, params.beta, params.w, output_prefix);
             *
             */

        } else if (mode == GEN_UINT_VEC) {
            assert_valid_nkv(params);
            assert_valid_ub(params);

            GenNNegIntVec(mg.GetLayerNumber(), params.nkv, params.ub.data(), output_prefix);

        } else if (mode == CMP_COHESIVENESS) {
            assert_valid_csg_file(params);
            ComputeCohesiveness(mg, params.csg_file, output_prefix);
        }
    }

    return 0;
}
