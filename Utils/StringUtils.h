//
// Created by ldd on 2023/6/13.
//

#ifndef MLCDEC_STRINGUTILS_H
#define MLCDEC_STRINGUTILS_H


static void Split(const string &str, vector<string> &tokens, bool skip_first, char delim = ' ') {
    size_t last_pos, pos;

    tokens.clear();

    last_pos = str.find_first_not_of(delim, 0);
    if (skip_first) {
        last_pos = str.find(delim, last_pos);
        last_pos = str.find_first_not_of(delim, last_pos);
    }

    while (last_pos != string::npos) {
        pos = str.find(delim, last_pos);
        tokens.emplace_back(move(str.substr(last_pos, pos - last_pos)));
        last_pos = str.find_first_not_of(delim, pos);
    }
}

static void Str2LLUVec(string s, vector<ll_uint> &vec, char delim = ',') {
    vector<string> str_arr;
    if (*s.begin() == '[') {
        s = s.substr(1, s.size() - 2);  // remove '[' and ']'
    }

    Split(s, str_arr, false, delim);
    for (auto &ss:str_arr) vec.emplace_back(std::stoll(ss));
}

static void Str2UIVec(string s, vector<uint> &vec, char delim = ',') {
    vector<string> str_arr;
    if (*s.begin() == '[') {
        s = s.substr(1, s.size() - 2);  // remove '[' and ']'
    }

    Split(s, str_arr, false, delim);
    for (auto &ss:str_arr) vec.emplace_back(stoi(ss));
}

static void Str2FVec(string s, vector<float> &vec, char delim = ',') {
    vector<string> str_arr;
    if (*s.begin() == '[') {
        s = s.substr(1, s.size() - 2);  // remove '[' and ']'
    }
    Split(s, str_arr, false, delim);
    for (auto &ss:str_arr) vec.emplace_back(stof(ss));
}

static string GetFileName(const string &path) {
    uint pos = path.rfind('/') + 1;
    return path.substr(pos, path.rfind('.') - pos);
}

#endif //MLCDEC_STRINGUTILS_H
