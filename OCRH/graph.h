#ifndef GRAPH_COMMON_H
#define GRAPH_COMMON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <string>
#include <vector>
#include <climits>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <omp.h>
#include <random>
#include <iomanip>
#include <atomic>
#include <cstdint>
#include <dirent.h>
#include <sys/types.h>


using namespace std;
using namespace chrono;


const int k = 5;
const int BIG_PRIME = 1000000007;
const int SIGNATURE_SIZE = 3;
const int TOP_K_HASHES = 3;



bool read_csr_direct(const string& filename, vector<int>& csr_offsets, vector<int>& csr_edges, unordered_map<int, int>& out_degree_map);

vector<int> bfs_order_csr(const vector<int>& csr_offsets, const vector<int>& csr_edges, int start_node);
vector<int> dfs_order_csr(const vector<int>& csr_offsets, const vector<int>& csr_edges, int start_node);
vector<int> sssp_order_csr(const vector<int>& csr_offsets, const vector<int>& csr_edges, int start_node);

vector<int> bfs_in_cluster_global(int start, int cid, const vector<int>& clusters, const vector<int>& off, const vector<int>& edges);
vector<int> dfs_in_cluster_global(int start, int cid, const vector<int>& clusters, const vector<int>& off, const vector<int>& edges);
vector<int> sssp_in_cluster_global(int start, int cid, const vector<int>& clusters, const vector<int>& off, const vector<int>& edges);

vector<pair<int, int>> generate_hash_params(int num_hashes);
uint64_t compute_kmer_hash(const vector<int>& kmer_vec, int a, int b);
unordered_map<int, vector<uint64_t>> compute_topk_minhash_kmers(const unordered_map<int, vector<uint64_t>>& node_to_kmers, const vector<pair<int, int>>& hash_params);
double jaccard_similarity_hash(const vector<uint64_t>& A, const vector<uint64_t>& B);

pair<vector<int>, vector<int>> kmeans_clustering(vector<int>& csr_offsets, vector<int>& csr_edges, const vector<int>& vertex_list, const unordered_map<int, vector<uint64_t>>& node_to_minhash, int k_clusters, int x_candidates, const unordered_map<int, int>& out_degree);
vector<pair<double, int>> sort_centers_by_similarity_with_score(const vector<int>& centers, const unordered_map<int, vector<uint64_t>>& node_to_minhash);


#endif