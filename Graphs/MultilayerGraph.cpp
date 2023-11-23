//
// Created by ldd on 2023/3/6.
//


#include "MultilayerGraph.h"


MultilayerGraph::~MultilayerGraph() {
    delete[] graph_layers;
    delete[] order;
}

void MultilayerGraph::LoadFromFile(const string &input_path) {
    vector<string> graph_files;
    vector<edge *> edge_buf;
    vector<uint> edge_size;

    map_file = input_path + "__vtx_map.txt";
    unordered_map<ll_uint, uint> vtx2id;
    std::basic_ofstream<char> map_file_out;

    GetGraphFile(input_path, graph_files);

    n_layers = graph_files.size();
    graph_layers = new Graph[n_layers];

    edge_buf.resize(n_layers);
    edge_size.resize(n_layers);

    map_file_out = ofstream(map_file);
    for (uint i = 0; i < n_layers; i++) {
        edge_size[i] = LoadLayer(input_path + graph_files[i], edge_buf[i], vtx2id, map_file_out);
    }
    map_file_out.close();

    n = vtx2id.size();
    for (uint i = 0; i < n_layers; i++) {
        graph_layers[i].BuildFromEdgeLst(edge_buf[i], n, edge_size[i]);
        delete[] edge_buf[i];
    }

    order = new uint[n_layers];
    for (uint i = 0; i < n_layers; i++) {
        order[i] = i;
    }
}

void MultilayerGraph::SetGraphOrder(G_ORDER ordering) {
    
    if (ordering == DEFAULT) return;

    if (ordering == DEN_INC || ordering == DEN_DEC) {

        float density[n_layers];
        for (uint i = 0; i < n_layers; i++) {
            density[i] = (float) (graph_layers[i].GetM() >> 1) / (float) graph_layers[i].GetN();
        }

        if (ordering == DEN_INC) {
            sort(order, order + n_layers, [&density](uint id1, uint id2) { return density[id1] < density[id2]; });
        } else {
            sort(order, order + n_layers, [&density](uint id1, uint id2) { return density[id1] > density[id2]; });
        }

        cout << "Graph order set to: [ ";
        for (uint i = 0; i < n_layers; i++) {
            cout << "l" << order[i] << "(" << density[order[i]] << ") ";
        }
        cout << "]." << endl;

    } else if (ordering == CORE_INC || ordering == CORE_DEC) {

        uint degeneracy[n_layers];
        for (uint i = 0; i < n_layers; i++) {
            degeneracy[i] = KC::Degeneracy(graph_layers[i]);
        }

        if (ordering == CORE_INC) {
            sort(order, order + n_layers,
                 [&degeneracy](uint id1, uint id2) { return degeneracy[id1] < degeneracy[id2]; });
        } else {
            sort(order, order + n_layers,
                 [&degeneracy](uint id1, uint id2) { return degeneracy[id1] > degeneracy[id2]; });
        }

        cout << "Graph order set to: [ ";
        for (uint i = 0; i < n_layers; i++) {
            cout << "l" << order[i] << "(" << degeneracy[order[i]] << ") ";
        }
        cout << "]." << endl;

    } else if (ordering == E_INC || ordering == E_DEC) {

        uint n_edges[n_layers];
        for (uint i = 0; i < n_layers; i++) {
            n_edges[i] = graph_layers[i].GetM();
        }

        if (ordering == E_INC) {
            sort(order, order + n_layers,
                 [&n_edges](uint id1, uint id2) { return n_edges[id1] < n_edges[id2]; });
        } else {
            sort(order, order + n_layers,
                 [&n_edges](uint id1, uint id2) { return n_edges[id1] > n_edges[id2]; });
        }

//        cout << "Graph order set to: [ ";
//        for (uint i = 0; i < n_layers; i++) {
//            cout << "l" << order[i] << "(" << n_edges[order[i]] << ") ";
//        }
//        cout << "]." << endl;
    }
}


void MultilayerGraph::PrintStatistics() {
    uint max_m = 0, sum_m = 0, m;

    for (uint i = 0; i < n_layers; i++) {
        m = graph_layers[i].GetM() >> 1;
        max_m = std::max(max_m, m);
        sum_m += m;
    }

    std::cout << "|L| = " << n_layers << ", |V| = " << n << ", |E| = " << sum_m << ", max_|E_i| = " << max_m << std::endl;
}


void MultilayerGraph::GetGraphFile(const string &graph_path, vector<string> &graph_files) {
    struct stat buffer{};
    string &&conf_file = graph_path + "mlg.conf";

    if (stat((conf_file).c_str(), &buffer)) {
        // no configure file provided

        DIR *dir;
        dirent *env;
        if ((dir = opendir(graph_path.c_str()))) {
            while ((env = readdir(dir))) {
                graph_files.emplace_back(env->d_name);
            }
            closedir(dir);
        }

    } else {
        // configure file provided
        string line;

        auto fin = ifstream(conf_file);
        while (fin.peek() != EOF) {

            getline(fin, line);

            if (!line.empty()) {
                graph_files.emplace_back(line);
            }
        }
        fin.close();
    }
}

uint MultilayerGraph::LoadLayer(const string &graph_file, edge *&edge_buf, unordered_map<ll_uint, uint> &vtx2id,
                                std::basic_ofstream<char> &map_file_out) {
    uint uid, vid;
    uint edge_buf_size, num_of_vtx, num_of_edge;
    ll_uint u, v;

    edge *tmp_edge_buf;
    std::basic_ifstream<char> graph_in;

    graph_in = ifstream(graph_file);

    num_of_vtx = (uint) vtx2id.size();

    // load graph
    num_of_edge = 0;
    edge_buf_size = DEFAULT_EDGE_BUF_SIZE;
    edge_buf = new edge[edge_buf_size];

    while (graph_in.peek() != EOF) {

        graph_in >> u >> v;

        // remove self-loop
        if (u != v) {
            auto iter1 = vtx2id.find(u);
            if (iter1 != vtx2id.end()) {
                uid = iter1->second;
            } else {
                uid = num_of_vtx++;
                vtx2id.emplace(u, uid);
                map_file_out << u << " " << uid << endl;
            }

            auto iter2 = vtx2id.find(v);
            if (iter2 != vtx2id.end()) {
                vid = iter2->second;
            } else {
                vid = num_of_vtx++;
                vtx2id.emplace(v, vid);
                map_file_out << v << " " << vid << endl;
            }

            if (num_of_edge + 2 > edge_buf_size) {
                edge_buf_size = edge_buf_size << 1;
                tmp_edge_buf = new edge[edge_buf_size];

                memcpy(tmp_edge_buf, edge_buf, num_of_edge * sizeof(edge));
                delete[] edge_buf;
                edge_buf = tmp_edge_buf;
            }

            edge_buf[num_of_edge++] = edge(uid, vid);
            edge_buf[num_of_edge++] = edge(vid, uid);
        }
    }

    graph_in.close();
    return num_of_edge;
}

void MultilayerGraph::LoadId2VtxMap(ll_uint *id2vtx) {
    ll_uint u;
    uint uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        id2vtx[uid] = u;
    }
    fin.close();
}

void MultilayerGraph::LoadVtx2IdMap(unordered_map<ll_uint, uint> &vtx2id) {
    ll_uint u;
    uint uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        vtx2id[u] = uid;
    }
    fin.close();
}


MultilayerGraph* MultilayerGraph::Load(const string &file) {
    auto mg = new MultilayerGraph();
    mg->LoadFromFile(file);
    return mg;
}
