//
// Created by ldd on 2023/3/6.
//

#ifndef MLCDEC_HEADER_H
#define MLCDEC_HEADER_H

#include <sys/stat.h>
#include <sys/resource.h>

#include <iomanip>
#include <iostream>
#include <fstream>

#include <dirent.h>

#include <limits>
#include <cfloat>

#include <functional>
#include <algorithm>
#include <random>
#include <cmath>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

#include <stdint.h>
#include <string>
#include <cstring>

#include <omp.h>
#include <thread>
#include <atomic>
#include <future>
#include <condition_variable>

#include <cassert>

using std::string;
using std::ifstream;
using std::ofstream;

using std::cerr;
using std::cout;

using std::stoi;
using std::sort;
using std::ceil;
using std::max;
using std::abs;
using std::endl;
using std::to_string;

using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using std::hash;
using std::mt19937;
using std::random_device;
using std::uniform_int_distribution;

typedef long long unsigned ll_uint;

#define EPSILON 1E-6

#ifdef USE_JEMALLOC
#include "jemalloc/jemalloc.h"
#endif

static void print(const string &info) {
    cout << info << endl;
}

#endif //MLCDEC_HEADER_H