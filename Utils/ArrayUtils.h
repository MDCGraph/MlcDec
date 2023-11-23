//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_ARRAYUTILS_H
#define MLCDEC_ARRAYUTILS_H

#include "../Header.h"

template<class T>
inline bool ArrGe(const T *arr1, const T *arr2, uint length) {
    for (uint i = 0; i < length; i++) {
        if (arr1[i] < arr2[i]) {
            return false;
        }
    }
    return true;
}

template<class T>
inline bool ArrEq(const T *arr1, const T *arr2, uint length) {
    for (uint i = 0; i < length; i++) {
        if (arr1[i] != arr2[i]) {
            return false;
        }
    }
    return true;
}

template<class T>
static string Arr2Str(const T *arr, uint length) {
    uint i;
    string s = "[";

    for (i = 0; i < length - 1; i++) {
        s += to_string(arr[i]) + ",";
    }

    if (length) s += to_string(arr[i]);

    s += "]";

    return s;  // format: [1,2,3,4]
}

template<class T>
static string Vec2Str(const vector<T> &arr) {
    uint i, length = (uint) arr.size();

    string s = "[";
    for (i = 0; i < length - 1; i++) {
        s += to_string(arr[i]) + ",";
    }
    if (length) s += to_string(arr[i]);

    s += "]";

    return s;
}

template<class T>
static uint FilterIndexGe(const T *arr, uint length, T *filtered_res, T threshold) {
    uint n_res = 0;
    for (uint i = 0; i < length; i++) {
        if (arr[i] >= threshold) {
            filtered_res[n_res++] = i;
        }
    }
    return n_res;
}

static void GetCombs(uint idx, uint *vec, uint len, const uint *bounds, uint step, vector<vector<uint>> &comb) {
    comb.emplace_back(vec, vec + len);
    for (uint i = idx; i < len; i++) {
        vec[i] += step;
        if (vec[i] < bounds[i]) GetCombs(i, vec, len, bounds, step, comb);
        vec[i] -= step;
    }
}

static void GetCombs(uint len, const uint *bounds, uint step, vector<vector<uint>> &comb) {
    uint vec[len];

    memset(vec, 0, len * sizeof(uint));
    GetCombs(0, vec, len, bounds, step, comb);
}

template<class T>
static uint GetNPos(const T *arr, uint length) {
    uint n = 0;
    for (uint i = 0; i < length; i++) {
        if (arr[i]) n += 1;
    }
    return n;
}

template<class T>
static uint GetNDiff(const T *arr1, const T *arr2, uint length) {
    uint n = 0;
    for (uint i = 0; i < length; i++) {
        if (arr1[i] != arr2[i]) n += 1;
    }
    return n;
}

static void SetSign(const uint *arr, uint length, bool *sign) {
    for (uint i = 0; i < length; i++) {
        sign[arr[i]] = true;
    }
}

static uint CountIntersect(const uint *arr1, uint length1, const uint *arr2, uint length2) {

    uint i = 0, j = 0, n = 0;

    while (i < length1 && j < length2) {
        if (arr1[i] > arr2[j]) j++;
        else if (arr1[i] < arr2[j]) i++;
        else {
            n += 1;
            i++, j++;
        }
    }
    return n;

}

#endif //MLCDEC_ARRAYUTILS_H
