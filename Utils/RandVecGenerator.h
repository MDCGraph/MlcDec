//
// Created by ldd on 2023/4/1.
//

#ifndef MLCDEC_RANDVECGENERATOR_H
#define MLCDEC_RANDVECGENERATOR_H

#include "../Header.h"

static void GenNNegIntVec(uint dim, uint size, const uint *ub, vector<vector<uint>> &vec) {

    uint k_vec[dim];

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint> dis[dim];

    for (uint i = 0; i < dim; i++) {
        dis[i] = std::uniform_int_distribution<uint>{0, ub[i]};
    }

    vec.reserve(size);
    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < dim; j++) k_vec[j] = dis[j](gen);
        vec.emplace_back(&k_vec[0], k_vec + dim);
    }
}

static void GenNNegIntVec(uint dim, uint size, const uint *ub, const string &output) {

    vector<vector<uint>> vec;
    string file;

    GenNNegIntVec(dim, size, ub, vec);

    file = output + "_kv_" + to_string(dim) + "_" + Arr2Str(ub, dim) + "_" + to_string(size) + ".txt";
    auto f = ofstream(file);

    f.write(reinterpret_cast<char *>(&size), sizeof(uint));
    for (uint i = 0; i < size; i++) {
        f.write(reinterpret_cast<char *>(vec[i].data()), long(dim * sizeof(uint)));
    }
    f.close();
}


static uint LoadNNegIntVec(uint dim, const string &file, vector<vector<uint>> &vec) {
    uint size;

    auto f = ifstream(file);
    f.read(reinterpret_cast<char *> (&size), sizeof(uint));

    vec.resize(size);
    for (uint i = 0; i < size; i++) {
        vec[i].resize(dim);
        f.read(reinterpret_cast<char *>(vec[i].data()), long(dim * sizeof(uint)));
    }
    f.close();

    return size;
}

#endif //MLCDEC_RANDVECGENERATOR_H
