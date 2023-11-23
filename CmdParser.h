//
// Created by ldd on 2023/3/28.
//

#ifndef MLCDEC_CMDPARSER_H
#define MLCDEC_CMDPARSER_H

#include "Utils/StringUtils.h"
#include "Params.h"
#include <getopt.h>

// running mode
#define PRT_GRAPH_INFO 0

#define RUN_MLC_SEARCH 8
#define RUN_MLCT_MLC_SEARCH 9
#define RUN_HT_MLC_SEARCH 10

#define BUILD_DS_MLC_TREE 11
#define BUILD_MLC_TREE 12
#define BUILD_MLC_TREE_PARALLEL 13
#define BUILD_CORE_CUBE 14
#define BUILD_MLC_HT 15

#define WDS_SEARCH 16
#define WDS_TREE_SEARCH 17
#define DCC_WDS_SEARCH 18

#define GEN_UINT_VEC 20
#define CMP_COHESIVENESS 21

#define PARAMETER_SENSITIVITY_TEST 30
#define SERIAL_MLCT_BUILD_TEST 31

// Parameter
#define PROVIDE_GRAPH_PATH 50
#define PROVIDE_GRAPH_NAME 51
#define PROVIDE_OUTPUT_PATH 52
#define PROVIDE_GRAPH_ORDER 53
#define PROVIDE_KVEC_FILE 54
#define PROVIDE_MLCT_FILE 55
#define PROVIDE_CC_FILE 56
#define PROVIDE_MLCHT_FILE 57

#define PROVIDE_THREAD_NUM 61
#define PROVIDE_ALPHA 62
#define PROVIDE_GAMMA 63
#define PROVIDE_CACHE_LEVEL 64
#define PROVIDE_BUCKET_SIZE 65

#define PROVIDE_BETA 70
#define PROVIDE_W 71
#define PROVIDE_NKV 72
#define PROVIDE_UB 73
#define PROVIDE_CSG_FILE 74

#define PROVIDE_MLCT_BUILDER 80
#define PROVIDE_PARALLEL_MLCT_BUILDER 81


#define PROVIDE_PRUN_MODE 85
#define PROVIDE_LEVEL 86

#define FLUSH 100


#define DIR_CREATE_MODE 0777


struct String2Mode {
    const char *name;
    int id;
};

static option long_options[] = {

        {"p",      required_argument, nullptr, PROVIDE_GRAPH_PATH},
        {"g",      required_argument, nullptr, PROVIDE_GRAPH_NAME},
        {"o",      required_argument, nullptr, PROVIDE_GRAPH_ORDER},
        {"output", required_argument, nullptr, PROVIDE_OUTPUT_PATH},
        {"kvf",    required_argument, nullptr, PROVIDE_KVEC_FILE},
        {"mlctf",  required_argument, nullptr, PROVIDE_MLCT_FILE},
        {"mlchtf", required_argument, nullptr, PROVIDE_MLCHT_FILE},
        {"ccf",    required_argument, nullptr, PROVIDE_CC_FILE},

        {"b",      required_argument, nullptr, PROVIDE_MLCT_BUILDER},
        {"pb",     required_argument, nullptr, PROVIDE_PARALLEL_MLCT_BUILDER},

        {"t",      required_argument, nullptr, PROVIDE_THREAD_NUM},
        {"alpha",  required_argument, nullptr, PROVIDE_ALPHA},
        {"gamma",  required_argument, nullptr, PROVIDE_GAMMA},
        {"cl",     required_argument, nullptr, PROVIDE_CACHE_LEVEL},
        {"bks",    required_argument, nullptr, PROVIDE_BUCKET_SIZE},

        {"beta",   required_argument, nullptr, PROVIDE_BETA},
        {"w",      required_argument, nullptr, PROVIDE_W},

        {"nkv",    required_argument, nullptr, PROVIDE_NKV},
        {"ub",     required_argument, nullptr, PROVIDE_UB},
        {"csgf",   required_argument, nullptr, PROVIDE_CSG_FILE},

        {"m",      required_argument, nullptr, PROVIDE_PRUN_MODE},
        {"l",      required_argument, nullptr, PROVIDE_LEVEL},
        {"f",      no_argument,       nullptr, FLUSH},
};

const String2Mode execution_mode_options[] = {
        {"info",   PRT_GRAPH_INFO},

        {"mlcs",   RUN_MLC_SEARCH},
        {"mlcs_t", RUN_MLCT_MLC_SEARCH},
        {"mlcs_h", RUN_HT_MLC_SEARCH},

        {"bcc",    BUILD_CORE_CUBE},
        {"bt",     BUILD_MLC_TREE},
        {"bdt",    BUILD_DS_MLC_TREE},
        {"btp",    BUILD_MLC_TREE_PARALLEL},
        {"bht",    BUILD_MLC_HT},

        {"ds",     WDS_SEARCH},
        {"dts",    WDS_TREE_SEARCH},
        {"dds",    DCC_WDS_SEARCH},

        {"gk",     GEN_UINT_VEC},
        {"cm",     CMP_COHESIVENESS},

        {"param",  PARAMETER_SENSITIVITY_TEST},
        {"smlct",  SERIAL_MLCT_BUILD_TEST},

        {nullptr, -1}
};

const String2Mode graph_ordering_options[] = {
        {"d",     DEN_INC},
        {"d_r",   DEN_DEC},
        {"c",     CORE_INC},
        {"c_r",   CORE_DEC},
        {"e",     E_INC},
        {"e_r",   E_DEC},
        {"rd",    DEFAULT},
        {nullptr, -1}
};

const String2Mode mlc_tree_builder_options[] = {
        {"naive", NAIVE},
        {"ne",    NODE_ELIMINATION},
        {"se",    SUBTREE_ELIMINATION},
        {"sm",    SUBTREE_MERGE},
        {"opt",   OPT},

        {nullptr, -1}
};

const String2Mode prun_options[] = {
        {"s",      SERIAL},
        {"serial", SERIAL},
        {"ip",     IP},
        {"ips",    IP_START},
        {"ipm",    IP_MERGE},
        {"opt",    P_OPT},

        {nullptr,  -1}
};

const char *datasets[] = {"example", "homo", "obamainisrael", "sacchcere", "dblp", "higgs", "friendfeed",
                          "friendfeedtwitter", "amazon", "wiki", "dblp-coauthor", "flickr-growth", nullptr};

static int GetMode(const String2Mode *options, const char *str_mode) {
    for (int i = 0; options[i].name; i++) {
        if (!strcasecmp(options[i].name, str_mode)) {
            return options[i].id;
        }
    }
    return -1;
}


void ParseCmd(Params &params, int argc, char *argv[]) {
    int c, option_index;

    // read options
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case PROVIDE_GRAPH_PATH:
                if (optarg) params.path = optarg;
                break;
            case PROVIDE_GRAPH_NAME:
                if (optarg) params.dataset = optarg;
                break;
            case PROVIDE_OUTPUT_PATH:
                if (optarg) params.output = optarg;
                break;
            case PROVIDE_GRAPH_ORDER:
                if (optarg) {
                    auto oid = GetMode(graph_ordering_options, optarg);
                    if (oid >= 0) params.graph_ordering = static_cast<G_ORDER>(oid);
                    else {
                        cerr << "Invalid graph ordering." << endl;
                        exit(-1);
                    }
                }
                break;
            case PROVIDE_KVEC_FILE:
                if (optarg) params.kvf = optarg;
                break;
            case PROVIDE_MLCT_FILE:
                if (optarg) params.mlct_file = optarg;
                break;
            case PROVIDE_MLCHT_FILE:
                if (optarg) params.mlcht_file = optarg;
                break;
            case PROVIDE_CC_FILE:
                if (optarg) params.cc_file = optarg;
                break;
            case PROVIDE_THREAD_NUM:
                if (optarg) params.nt = std::stoi(optarg);
                break;
            case PROVIDE_CACHE_LEVEL:
                if (optarg) params.cl = std::stoi(optarg);
                break;
            case PROVIDE_ALPHA:
                if (optarg) params.alpha = std::stof(optarg);
                break;
            case PROVIDE_GAMMA:
                if (optarg) params.gamma = std::stof(optarg);
                break;
            case PROVIDE_BETA:
                if (optarg) params.beta = std::stof(optarg);
                break;
            case PROVIDE_W:
                if (optarg) Str2FVec(optarg, params.w);
                break;
            case PROVIDE_NKV:
                if (optarg) params.nkv = stoi(optarg);
                break;
            case PROVIDE_UB:
                if (optarg) Str2UIVec(optarg, params.ub);
                break;
            case PROVIDE_CSG_FILE:
                if (optarg) params.csg_file = optarg;
                break;
            case PROVIDE_MLCT_BUILDER:
                if (optarg) {
                    auto bid = GetMode(mlc_tree_builder_options, optarg);
                    if (bid >= 0) params.builder = static_cast<MLCTree_builder>(bid);
                    else {
                        cerr << "Invalid MLC-tree builder." << endl;
                        exit(-1);
                    }
                }
                break;
            case PROVIDE_PRUN_MODE:
                if (optarg) {
                    auto bid = GetMode(prun_options, optarg);
                    if (bid >= 0) params.prun_opt = static_cast<RUN_OPTION>(bid);
                    else {
                        cerr << "Invalid parallel run option." << endl;
                        exit(-1);
                    }
                }
                break;
            case PROVIDE_LEVEL:
                if (optarg) params.l = stoi(optarg);
                break;
            case PROVIDE_BUCKET_SIZE:
                if (optarg) params.bucket_size = stoi(optarg);
                break;
            case FLUSH:
                params.flush = true;
                break;
            default:
                cerr << "Invalid option." << endl;
                exit(-1);
        }
    }

    // read operands
    if (argc > optind) {
        while (true) {
            auto mode = GetMode(execution_mode_options, argv[optind]);
            if (mode == -1) {
                cerr << "Invalid execution mode." << endl;
                exit(-1);
            }

            params.mode.emplace_back(mode);
            if (argc == ++optind) break;
        }

    } else {
        cerr << "Execution mode is missing." << endl;
        exit(-1);
    }
}

static bool is_invalid_path(const string &file) {
    struct stat buffer{};

    return stat(file.c_str(), &buffer);
}

static bool is_valid_datasets(const char *dataset) {
    for (int i = 0; datasets[i]; i++) {
        if (!strcasecmp(datasets[i], dataset)) {
            return true;
        }
    }
    return false;
}

static void assert_non_empty_graph_name(Params &params) {
    if (params.dataset.empty()) {
        cerr << "Graph name is missing." << endl;
        exit(-1);
    }
}

static void assert_default_graph_name(Params &params) {
    if (!is_valid_datasets(params.dataset.c_str())) {
        cerr << "Graph " << params.dataset << " doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_graph_path(Params &params) {

    if (is_invalid_path(params.path)) {
        cerr << "Provided graph path \"" << params.path << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void warn_empty_output_path(Params &params) {

    if (params.output.empty()) {
        params.output = default_output;
        if (is_invalid_path(params.output)) {
            mkdir(params.output.c_str(), DIR_CREATE_MODE);
        }
        cerr << "Output path is not provided, and is set to \"" << default_output << "\" by default." << endl;
    } else if (is_invalid_path(params.output)) {
        cerr << "Provided output path \"" << params.output << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_kvec_file(Params &params) {

    if (params.kvf.empty()) {
        cerr << "Coreness vector file should be provided." << endl;
        exit(-1);
    } else if (is_invalid_path(params.kvf)) {
        cerr << "Provided coreness vector file \"" << params.kvf << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_tree_file(Params &params) {

    if (params.mlct_file.empty()) {
        cerr << "MLC-tree file should be provided." << endl;
        exit(-1);
    } else if (is_invalid_path(params.mlct_file)) {
        cerr << "Provided MLC-tree file \"" << params.mlct_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}


static void assert_valid_ht_file(Params &params) {

    if (params.mlcht_file.empty()) {
        cerr << "MLC-HT file should be provided." << endl;
        exit(-1);
    } else if (is_invalid_path(params.mlcht_file)) {
        cerr << "Provided MLC-HT file \"" << params.mlcht_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}


static void assert_valid_cc_file(Params &params) {

    if (params.cc_file.empty()) {
        cerr << "Corecube file should be provided." << endl;
        exit(-1);
    } else if (is_invalid_path(params.cc_file)) {
        cerr << "Provided Corecube file \"" << params.cc_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_csg_file(Params &params) {
    if (params.csg_file.empty()) {
        cerr << "Cohesive subgraph vertex file should be provided." << endl;
        exit(-1);
    } else if (is_invalid_path(params.csg_file)) {
        cerr << "Provided cohesive subgraph vertex file \"" << params.csg_file << "\" doesn't exist." << endl;
        exit(-1);
    }
}

static void assert_valid_nkv(Params &params) {
    if (!params.nkv) {
        params.nkv = default_nkv;
        cerr << "The number of the generated coreness vectors is not provided, and is set to " << default_nkv
             << " by default." << endl;

    }
}

static void assert_valid_ub(Params &params) {
    if (params.ub.empty()) {
        cerr << "The upper bound of the generated coreness vectors should be provided." << endl;
        exit(-1);
    }
}

static void assert_valid_nt(Params &params) {
    if (!params.nt) {
        cerr << "The number of threads should be provided." << endl;
        exit(-1);
    }
}

static void assert_valid_alpha(Params &params) {
    if (params.alpha < EPSILON) {
        params.alpha = default_alpha;
        cerr << "The value of alpha is not provided, and is set to " << default_alpha << " by default." << endl;
    }
}

static void assert_valid_l(Params &params) {
    if (params.l < 0) {
        params.l = default_l;
        cerr << "The value of level is not provided, and is set to " << default_l << " by default." << endl;
    }
}


static void assert_valid_beta(Params &params) {
    if (params.beta < -EPSILON) {
        cerr << "The value of beta should be provided." << endl;
        exit(-1);
    }
}


static void assert_valid_w(Params &params) {
    if (params.w.empty()) {
        cerr << "The weight vector should be provided." << endl;
        exit(-1);
    }
}


static void verify_prun_param(Params &params) {
    if (params.prun_opt == P_OPT) {
        assert_valid_alpha(params);
        assert_valid_l(params);
    } else if (params.prun_opt == IP_START) {
        assert_valid_l(params);
        if (params.alpha > EPSILON) params.alpha = 0;
    } else if (params.prun_opt == IP_MERGE) {
        assert_valid_alpha(params);
        if (params.l > 0) params.l = -1;
    } else if (params.prun_opt == IP) {
        if (params.alpha > EPSILON) params.alpha = 0;
        if (params.l > 0) params.l = -1;
    }
}

#endif //MLCDEC_CMDPARSER_H
