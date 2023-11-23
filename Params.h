//
// Created by ldd on 2023/3/28.
//

#ifndef MLCDEC_PARAMS_H
#define MLCDEC_PARAMS_H


#include "Core/MLCTBuilder.h"

struct Params {

    // execution mode
    vector<uint> mode;

    // arguments
    string default_path = "../datasets/layer/";
    string path;
    string dataset;
    string output;

    string kvf;
    string cc_file;
    string mlct_file;
    string mlcht_file;

    string csg_file;

    uint nkv = 0;
    vector<uint> ub;

    uint nt = -1;
    uint cl = 0;

    float gamma = 0;

    float beta = 0;
    vector<float> w;

    bool flush = false;

    G_ORDER graph_ordering = DEFAULT; // default

    MLCTree_builder builder;

    // new
    RUN_OPTION prun_opt = SERIAL;
    uint l = -1;
    float alpha = 0;
    uint bucket_size = 1024;
};



const char *default_output = "../output/";
const uint default_nkv = 1000;
const uint default_l= 2;
const float default_alpha = 0.1;


#endif //MLCDEC_PARAMS_H
