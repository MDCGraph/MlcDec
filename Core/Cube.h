//
// Created by ldd on 2023/6/10.
//

#ifndef MLCDEC_CUBE_H
#define MLCDEC_CUBE_H

class Cube {

public:
    Cube() = default;

    ~Cube() {
        if (coreness) {
            for (uint i = 0; i < n_comb; i++) {
                delete[] coreness[i];
            }
            delete[] coreness;
        }
    }

    void Set(uint n_, uint n_comb_) {
        n = n_;
        n_comb = n_comb_;
        coreness = new uint *[n_comb];
        for (uint i = 0; i < n_comb; i++) {
            coreness[i] = new uint[n];
        }
    }

    uint *GetCoreness(uint id) {
        return coreness[id];
    }

private:
    uint n;
    uint n_comb;
    uint **coreness;
};

typedef pair<uint, uint> vc_pair;  // vertex coreness pair

struct HYB_Cube1d {
    uint pre{};
    uint n_vcp{};
    vc_pair *vtx_cn{nullptr};

    ~HYB_Cube1d() {
        delete[] vtx_cn;
    }

    void set(uint pre_, uint n_vcp_) {
        pre = pre_;
        n_vcp = n_vcp_;
        vtx_cn = new vc_pair[n_vcp];
    }
};

class HYB_Cube {

public:
    HYB_Cube() = default;

    ~HYB_Cube() {
        delete[] coreness;
    }

    [[nodiscard]] uint GetN() const {
        return n;
    }

    [[nodiscard]] uint GetNComb() const {
        return n_comb;
    }

    void Set(uint n_, uint n_comb_) {
        n = n_;
        n_comb = n_comb_;
        coreness = new HYB_Cube1d[n_comb];
    }

    HYB_Cube1d &GetHybCoreness(uint id) {
        return coreness[id];
    }

    void GetCoreness(vector<uint> &sorted_ls, uint *cn) const {
        uint id;
        HYB_Cube1d *hyb_cube1d;

        memset(cn, 0, n * sizeof(uint));
        id = GetCoreOffset(sorted_ls);

        hyb_cube1d = &coreness[id];
        while (true) {

            for (uint i = 0; i < hyb_cube1d->n_vcp; i++) {
                cn[hyb_cube1d->vtx_cn[i].first] += hyb_cube1d->vtx_cn[i].second;
            }

            if (hyb_cube1d->pre != -1) {
                hyb_cube1d = &coreness[hyb_cube1d->pre];
            } else break;
        }
    }

    uint Search(const vector<uint> &vec, uint *core) {
        uint cn[n], length;
        vector<uint> sorted_ls;

        for (uint i = 0; i < vec.size(); i++) {
            if (vec[i]) sorted_ls.emplace_back(i);
        }
        GetCoreness(sorted_ls, cn);

        length = 0;
        for (uint i = 0; i < n; i++) {
            if (cn[i] >= vec[sorted_ls[0]]) {
                core[length++] = i;
            }
        }
        return length;
    }

    void Flush(const string &file) {
        auto f = std::ofstream(file);

        f.write(reinterpret_cast<char *> (&n), sizeof(uint));
        f.write(reinterpret_cast<char *> (&n_comb), sizeof(uint));

        for (uint i = 0; i < n_comb; i++) {
            auto &hyb_cube1d = coreness[i];
            f.write(reinterpret_cast<char *> (&hyb_cube1d.pre), sizeof(uint));
            f.write(reinterpret_cast<char *> (&hyb_cube1d.n_vcp), sizeof(uint));
            f.write(reinterpret_cast<char *> (hyb_cube1d.vtx_cn), long(hyb_cube1d.n_vcp * sizeof(uint)));
        }

        f.close();
    }

    static HYB_Cube *Load(const string &file) {
        uint n_, n_comb_, pre, n_vcp;

        auto hyb_cube = new HYB_Cube;

        auto f = std::ifstream(file);
        f.read(reinterpret_cast<char *> (&n_), sizeof(uint));
        f.read(reinterpret_cast<char *> (&n_comb_), sizeof(uint));

        hyb_cube->Set(n_, n_comb_);

        for (uint i = 0; i < n_comb_; i++) {
            auto &hyb_cube1d = hyb_cube->coreness[i];

            f.read(reinterpret_cast<char *> (&pre), sizeof(uint));
            f.read(reinterpret_cast<char *> (&n_vcp), sizeof(uint));
            hyb_cube1d.set(pre, n_vcp);

            f.read(reinterpret_cast<char *> (hyb_cube1d.vtx_cn), long(n_vcp * sizeof(uint)));
        }

        f.close();

        return hyb_cube;
    }

    bool Equals(HYB_Cube &hyb_cube) {
        if (n != hyb_cube.n || n_comb != hyb_cube.n_comb) {
            return false;
        }

        for (uint i = 0; i < n_comb; i++) {
            auto &hyb_1d = coreness[i];
            auto &hyb_cube_hyb_1d = coreness[i];

            if (hyb_1d.pre != hyb_cube_hyb_1d.pre || hyb_1d.vtx_cn != hyb_cube_hyb_1d.vtx_cn) {
                return false;
            }

            sort(hyb_1d.vtx_cn, hyb_1d.vtx_cn + hyb_1d.n_vcp);
            sort(hyb_cube_hyb_1d.vtx_cn, hyb_cube_hyb_1d.vtx_cn + hyb_1d.n_vcp);

            if (memcmp(hyb_1d.vtx_cn, hyb_cube_hyb_1d.vtx_cn, hyb_1d.n_vcp * sizeof(vc_pair)) != 0) {
                return false;
            }
        }
        return true;
    }

private:
    uint n;
    uint n_comb;
    HYB_Cube1d *coreness;

    static uint GetCoreOffset(const vector<uint> &sorted_layers) {
        uint sum = 0;

        for (auto i:sorted_layers) sum += (uint) pow(2, i);
        return sum - 1;
    }
};


#endif //MLCDEC_CUBE_H
