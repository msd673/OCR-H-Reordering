#ifndef ONLINE_GRAPH_COMMON_H
#define ONLINE_GRAPH_COMMON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include <string>
#include <vector>
#include <dirent.h>
#include <sys/types.h>
#include <climits>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <omp.h>
#include <random>
#include <stack>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <sys/syscall.h>
#include <linux/perf_event.h>
#include <cstring>
#include <sys/ioctl.h>
#include <unistd.h>
#include <tuple>
#include <iomanip>
#include <errno.h>
#include <functional>

using namespace std;
using namespace chrono;

const double DAMPING_FACTOR = 0.85;
const double PR_EPS = 1e-6;
const int MAX_PR_ITER = 20;


struct CSRGraph {
    vector<int> row_ptr;
    vector<int> col_idx;
    int num_nodes;

    CSRGraph(int n = 0);
    CSRGraph(vector<int>& offsets, vector<int>& edges);
};

struct UnionFind {
    vector<int> parent;
    vector<int> sz;

    UnionFind(int n) : parent(n), sz(n, 1) {
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    bool unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return false;

        if (sz[ra] < sz[rb]) {
            swap(ra, rb);
        }

        parent[rb] = ra;
        sz[ra] += sz[rb];
        return true;
    }
};
                    
bool read_csr_direct(const string& filename,
                     vector<int>& csr_offsets,
                     vector<int>& csr_edges);

vector<int> bfs_order_csr(const CSRGraph& graph, int start);
vector<int> dfs_order_csr(const CSRGraph& graph, int start);
vector<int> sssp_order_csr(const vector<int>& csr_offsets, const vector<int>& csr_edges, int start_node);

vector<int> compute_wcc_unionfind(const CSRGraph& graph);
CSRGraph build_reverse_graph(const CSRGraph& graph);
vector<int> compute_wcc_bfs_only(const CSRGraph& graph, const CSRGraph& rev_graph);
vector<int> compute_wcc(const CSRGraph& graph);

vector<double> compute_pagerank(const CSRGraph& graph);
vector<int> compute_pagerank_order(const CSRGraph& graph);

int bfs_eccentricity(const CSRGraph& graph, int src);
int compute_diameter(const CSRGraph& graph, const vector<int>& sources);
int compute_diameter_with_sources_parallel(const CSRGraph& graph, const vector<int>& sources);


vector<int> sort_nodes_by_degree(const CSRGraph& graph);
vector<int> get_random_vertices(const CSRGraph& graph, int count);
vector<int> get_top_degree_vertices(const CSRGraph& graph, int count);
vector<int> read_vertices_from_file(const string& filename);

pair<CSRGraph, unordered_map<int, int>> reorder_by_vector(const CSRGraph& graph, const vector<int>& order);


struct CachePerfResult {
    long long l1_ref;       
    long long l1_miss;      
    long long llc_ref;      
    long long llc_miss;     
    long long cache_miss_total; 
    long long execution_time;
};

int perf_event_open_wrapper(struct perf_event_attr* hw_event, pid_t pid,
                            int cpu, int group_fd, unsigned long flags, const string& event_name);
CachePerfResult measure_cache_perf(function<void()> workload);
double measure_avg_time(int repeat, function<void()> func);

vector<string> get_all_files(const string& dir_path);

#endif