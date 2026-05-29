#include "graph.h"

int main(int argc, char* argv[]) {
    auto total_start = high_resolution_clock::now();

    
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <graph_file> <topk> <x> <k_clusters> <traversal_type> [-large]" << endl;
        return 1;
    }

    string fileName = argv[1];
    int topk = stoi(argv[2]);
    int x_candidates = stoi(argv[3]);
    int k_clusters = stoi(argv[4]);
    int traversal_type = stoi(argv[5]);  
    string scale_type = "-small";
    if (argc >= 7) {
        scale_type = argv[6];
        if (scale_type != "-large" && scale_type != "-small") {
            cerr << "Invalid scale flag: " << scale_type << ". Use -large for large graphs; small graph mode is the default." << endl;
            return 1;
        }
    }
   

    vector<int> csr_offsets, csr_edges;
    unordered_map<int, int> out_degree_map;
    int V;

    auto stage1_start = high_resolution_clock::now();
    cout << "[Stage 1] Starting to read graph file..." << endl;
    read_csr_direct(fileName, csr_offsets, csr_edges, out_degree_map);
    V = csr_offsets.size() - 1;
    auto stage1_end = high_resolution_clock::now();
    cout << "[Stage 1] Graph reading completed, time elapsed: "
        << duration_cast<milliseconds>(stage1_end - stage1_start).count()
        << "ms" << endl;

    auto stage2_start = high_resolution_clock::now();
    cout << "[Stage 2] Selecting high-degree vertices as traversal start points..." << endl;
    vector<int> top_vertices;
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
    for (const auto& [id, deg] : out_degree_map) {
        if (pq.size() < (size_t)topk) {
            pq.emplace(deg, id);
        } else if (deg > pq.top().first) {
            pq.pop();
            pq.emplace(deg, id);
        }
    }
    while (!pq.empty()) {
        top_vertices.push_back(pq.top().second);
        pq.pop();
    }
    reverse(top_vertices.begin(), top_vertices.end());
    auto stage2_end = high_resolution_clock::now();
    cout << "[Stage 2] Completed, time elapsed: "
        << duration_cast<milliseconds>(stage2_end - stage2_start).count()
        << "ms" << endl;

    auto stage3_start = high_resolution_clock::now();
    int m = top_vertices.size();
    vector<vector<int>> node_to_bfs_positions(V, vector<int>(m, -1));
    cout << "[Stage 3] Generating traversal order for " << m << " start points (type=" << traversal_type << ")..." << endl;

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < m; ++i) {
        int start = top_vertices[i];
        vector<int> order;

        if (traversal_type == 0) {
            order = bfs_order_csr(csr_offsets, csr_edges, start);
        } else if (traversal_type == 1) {
            order = dfs_order_csr(csr_offsets, csr_edges, start);
        } else if (traversal_type == 2) {
            order = sssp_order_csr(csr_offsets, csr_edges, start);
        }

        for (int j = 0; j < (int)order.size(); ++j) {
            int v = order[j];
            node_to_bfs_positions[v][i] = j;
        }
    }
    auto stage3_end = high_resolution_clock::now();
    cout << "[Stage 3] Parallel traversal completed, time elapsed: "
        << duration_cast<milliseconds>(stage3_end - stage3_start).count()
        << "ms" << endl;

    vector<int> global_sorted_nodes;
    unordered_map<int, vector<uint64_t>> node_to_minhash_kmers;

    if (scale_type ==  "-large")
    {
        
        auto stage456_start = high_resolution_clock::now();
    cout << "\nCombining stages 4+5+6, computing MinHash in one pass..." << endl;
    node_to_minhash_kmers.reserve(V);

        const int kmer_size = k;
        auto hash_params = generate_hash_params(SIGNATURE_SIZE);
        auto kmer_hash_param = hash_params.empty() ? make_pair(1, 0) : hash_params[0];
        int a_kmer = kmer_hash_param.first;
        int b_kmer = kmer_hash_param.second;

        #pragma omp parallel
        {
            unordered_map<int, vector<uint64_t>> local_map;
            local_map.reserve(1024);

            #pragma omp for schedule(dynamic)
            for (int i = 0; i < V; ++i) {
                const auto& positions = node_to_bfs_positions[i];
                int size = positions.size();

                vector<int> indices(size);
                iota(indices.begin(), indices.end(), 0);
                sort(indices.begin(), indices.end(), [&](int a, int b) {
                    int pos_a = positions[a];
                    int pos_b = positions[b];
                    pos_a = (pos_a == -1) ? INT_MAX : pos_a;
                    pos_b = (pos_b == -1) ? INT_MAX : pos_b;
                    return pos_a < pos_b;
                });

                vector<uint64_t> sketch;
                sketch.reserve(hash_params.size() * TOP_K_HASHES);
                int n = indices.size();

                for (const auto& param : hash_params) {
                    vector<pair<uint64_t, uint64_t>> hashed_kmers;
                    vector<int> kmer_vec(kmer_size);

                    for (int j = 0; j < n; ++j) {
                        for (int l = 0; l < kmer_size; ++l) {
                            int idx = (j + l) % n;
                            kmer_vec[l] = indices[idx];
                        }
                        sort(kmer_vec.begin(), kmer_vec.end());
                        uint64_t kmer_hash = compute_kmer_hash(kmer_vec, a_kmer, b_kmer);
                        uint64_t minhash_val = (uint64_t(param.first) * kmer_hash + param.second) % BIG_PRIME;
                        hashed_kmers.emplace_back(minhash_val, kmer_hash);
                    }

                    if ((int)hashed_kmers.size() > TOP_K_HASHES) {
                        nth_element(hashed_kmers.begin(), hashed_kmers.begin() + TOP_K_HASHES, hashed_kmers.end());
                        for (int t = 0; t < TOP_K_HASHES; ++t) {
                            sketch.push_back(hashed_kmers[t].second);
                        }
                    } else {
                        for (const auto& [_, kmer_hash] : hashed_kmers) {
                            sketch.push_back(kmer_hash);
                        }
                    }
                }

                sort(sketch.begin(), sketch.end());
                sketch.erase(unique(sketch.begin(), sketch.end()), sketch.end());
                local_map[i] = move(sketch);
            }

            #pragma omp critical
            {
                for (auto& p : local_map) {
                    node_to_minhash_kmers[p.first] = move(p.second);
                }
            }
        }

        auto stage456_end = high_resolution_clock::now();
        cout << "[Stages 4+5+6] Completed, time elapsed: "
            << duration_cast<milliseconds>(stage456_end - stage456_start).count()
            << "ms" << endl;
    }
    else
    {
        auto stage4_start = high_resolution_clock::now();
        cout << "\n[Stage 4] Sorting the results" << endl;
        unordered_map<int, vector<int>> node_to_sorted_bfs_index;
        node_to_sorted_bfs_index.reserve(V);

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < V; ++i) {
            const auto& positions = node_to_bfs_positions[i];
            int size = positions.size();
            vector<int> indices(size);
            iota(indices.begin(), indices.end(), 0);
            sort(indices.begin(), indices.end(), [&](int a, int b) {
                int pos_a = positions[a];
                int pos_b = positions[b];
                pos_a = (pos_a == -1) ? INT_MAX : pos_a;
                pos_b = (pos_b == -1) ? INT_MAX : pos_b;
                return pos_a < pos_b;
            });
        #pragma omp critical
            node_to_sorted_bfs_index[i] = move(indices);
        }

        auto stage4_end = high_resolution_clock::now();
        cout << "[Stage 4] Completed, time elapsed: "
            << duration_cast<milliseconds>(stage4_end - stage4_start).count()
            << "ms" << endl;

        auto stage5_start = high_resolution_clock::now();
        cout << "[Stage 5] Generating k-mer hashes (parallel)..." << endl;

        unordered_map<int, vector<uint64_t>> node_to_kmers;
        node_to_kmers.reserve(node_to_sorted_bfs_index.size());

        const int kmer_size = k;
        auto hash_params = generate_hash_params(SIGNATURE_SIZE);
        auto kmer_hash_param = hash_params.empty() ? make_pair(1, 0) : hash_params[0];
        int a_kmer = kmer_hash_param.first;
        int b_kmer = kmer_hash_param.second;

        vector<int> node_ids;
        node_ids.reserve(node_to_sorted_bfs_index.size());
        for (const auto& [id, _] : node_to_sorted_bfs_index)
            node_ids.push_back(id);

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < (int)node_ids.size(); ++i) {
            int node_id = node_ids[i];
            const auto& sorted_bfs_index = node_to_sorted_bfs_index.at(node_id);
            int n = sorted_bfs_index.size();

            vector<uint64_t> kmers;
            kmers.reserve(n);
            vector<int> kmer_vec(kmer_size);

            for (int j = 0; j < n; ++j) {
                for (int l = 0; l < kmer_size; ++l) {
                    int idx = (j + l) % n;
                    kmer_vec[l] = sorted_bfs_index[idx];
                }
                sort(kmer_vec.begin(), kmer_vec.end());
                uint64_t kmer_hash = compute_kmer_hash(kmer_vec, a_kmer, b_kmer);
                kmers.push_back(kmer_hash);
            }
            
            #pragma omp critical
            node_to_kmers[node_id] = move(kmers);
        }

        auto stage5_end = high_resolution_clock::now();
        cout << "[Stage 5] k-mer hash generation completed, time elapsed: "
            << duration_cast<milliseconds>(stage5_end - stage5_start).count()
            << "ms" << endl;

        auto stage6_start = high_resolution_clock::now();
        cout << "[Stage 6] Computing MinHash signatures..." << endl;
        node_to_minhash_kmers = compute_topk_minhash_kmers(node_to_kmers, hash_params);
        auto stage6_end = high_resolution_clock::now();
        cout << "[Stage 6] Completed, time elapsed: "
            << duration_cast<milliseconds>(stage6_end - stage6_start).count()
            << "ms" << endl;
    }

    auto stage7_start = high_resolution_clock::now();
    cout << "[Stage 7] Starting k-means clustering, number of clusters:  " << k_clusters << "..." << endl;
    vector<int> vertex_list(V);
    iota(vertex_list.begin(), vertex_list.end(), 0);
    auto [clusters, centers] = kmeans_clustering(csr_offsets, csr_edges, vertex_list, node_to_minhash_kmers, k_clusters, x_candidates, out_degree_map);
    auto stage7_end = high_resolution_clock::now();
    cout << "[Stage 7] Clustering completed, time elapsed: " 
        << duration_cast<milliseconds>(stage7_end - stage7_start).count() 
        << "ms" << endl;

    auto stage8_start = high_resolution_clock::now();
    cout << "[Stage 8] Sorting clustering results (traversal type=" << traversal_type << ")" << endl;
    vector<pair<double, int>> center_order = sort_centers_by_similarity_with_score(centers, node_to_minhash_kmers);
    vector<vector<int>> cluster_orders(center_order.size());
    global_sorted_nodes.reserve(V);
    int center_order_count = static_cast<int>(center_order.size());

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < center_order_count; ++i) {
        int cluster_id = center_order[i].second;
        int center_node = centers[cluster_id];

        if (traversal_type == 0) {
            cluster_orders[i] = bfs_in_cluster_global(center_node, cluster_id, clusters, csr_offsets, csr_edges);
        } else if (traversal_type == 1) {
            cluster_orders[i] = dfs_in_cluster_global(center_node, cluster_id, clusters, csr_offsets, csr_edges);
        } else if (traversal_type == 2) {
            cluster_orders[i] = sssp_in_cluster_global(center_node, cluster_id, clusters, csr_offsets, csr_edges);
        }
    }

    for (int i = 0; i < center_order_count; ++i) {
        global_sorted_nodes.insert(global_sorted_nodes.end(),
                          cluster_orders[i].begin(),
                          cluster_orders[i].end());
    }
    auto stage8_end = high_resolution_clock::now();
    cout << "[Stage 8] Completed, time elapsed: " 
        << duration_cast<milliseconds>(stage8_end - stage8_start).count() 
        << "ms" << endl;

    
    auto stage9_start = high_resolution_clock::now();
    cout << "[Stage 9] Outputting results..." << endl;
    size_t pos = fileName.find_last_of("/\\");
    string baseName = (pos == string::npos) ? fileName : fileName.substr(pos + 1);
    size_t dot_pos = baseName.find_last_of(".");
    string name_no_ext = (dot_pos == string::npos) ? baseName : baseName.substr(0, dot_pos);

    string type_str = traversal_type == 0 ? "BFS" : (traversal_type == 1 ? "DFS" : "SSSP");
    string outputFileName = name_no_ext + "_" + type_str + "_kmeans";
    if (scale_type == "-large") {
        outputFileName += "_large";
    }
    outputFileName += ".txt";

    ofstream outfile(outputFileName);
    for (int node : global_sorted_nodes) {
        outfile << node << "\n";
    }
    outfile.close();
    auto stage9_end = high_resolution_clock::now();
    cout << "[Stage 9] Completed, time elapsed: " 
         << duration_cast<milliseconds>(stage9_end - stage9_start).count() 
         << "ms" << endl;

    
    auto total_end = high_resolution_clock::now();
    cout << "\n[Summary] All processing completed, total time elapsed: " 
         << duration_cast<milliseconds>(total_end - total_start).count() 
         << "ms" << endl;

    return 0;
}