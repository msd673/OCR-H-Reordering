#include "graph.h"


pair<vector<int>, vector<int>> kmeans_clustering(
    vector<int>& csr_offsets,
    vector<int>& csr_edges,
    const vector<int>& vertex_list,
    const unordered_map<int, vector<uint64_t>>& node_to_minhash_kmers,
    int k_clusters,
    int x_candidates,
    const unordered_map<int, int>& out_degree_map
) {
    using namespace std::chrono;
    auto total_start = high_resolution_clock::now();

    const int n = vertex_list.size();
    vector<int> clusters(n, -1);
    vector<int> centers;
    mt19937 rng(random_device{}());
    
    
    auto init_start = high_resolution_clock::now();
    centers.reserve(k_clusters);
    vector<int> candidates = vertex_list;
    shuffle(candidates.begin(), candidates.end(), rng);
    if (candidates.size() > static_cast<size_t>(k_clusters)) {
        centers.assign(candidates.begin(), candidates.begin() + k_clusters);
    } else {
        centers = candidates;
        k_clusters = centers.size();
        cout << "[K-means]  Insufficient number of nodes, number of clusters adjusted to:  " << k_clusters << endl;
    }
    auto init_end = high_resolution_clock::now();
    cout << "[K-means] Initial center selection time:  " 
         << duration_cast<milliseconds>(init_end - init_start).count() << "ms" << endl;
   
    const int max_iters = 200;
    const double stop_change_ratio = 0.05;

    
    vector<vector<int>> cluster_nodes(k_clusters);
    for (auto& group : cluster_nodes) {
        group.reserve(n / k_clusters + 1);
    }

    bool changed = true;
    for (int iter = 0; iter < max_iters && changed; ++iter) {
        auto iter_start = high_resolution_clock::now();
        cout << "[K-means] Iteration " << iter + 1 << "/" << max_iters << "..." << endl;

        
        auto assign_start = high_resolution_clock::now();
        atomic<int> change_cnt(0);
        vector<int> local_clusters(n, -1);

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            const int node = vertex_list[i];
            double best_sim = -1.0;
            int best_cluster = -1;

            for (int j = 0; j < k_clusters; ++j) {
                if (node == centers[j]) {
                    best_sim = 1.0;
                    best_cluster = j;
                    break;
                }
            }
            if (best_cluster != -1) {
                local_clusters[i] = best_cluster;
                if (clusters[i] != best_cluster) change_cnt++;
                continue;
            }

            
            const auto& node_kmer = node_to_minhash_kmers.at(node);
            for (int j = 0; j < k_clusters; ++j) {
                double sim = jaccard_similarity_hash(
                    node_kmer,
                    node_to_minhash_kmers.at(centers[j])
                );
                if (sim > best_sim) {
                    best_sim = sim;
                    best_cluster = j;
                }
            }

            local_clusters[i] = best_cluster;
            if (clusters[i] != best_cluster) change_cnt++;
        }

        auto assign_end = high_resolution_clock::now();
        double change_ratio = (double)change_cnt / n;
        cout << "[K-means] Node assignment time: "
            << duration_cast<milliseconds>(assign_end - assign_start).count()
            << "ms, changed nodes: " << change_cnt
            << " (" << fixed << setprecision(2) << 100.0 * change_ratio << "%)" << endl;

        if (change_ratio <= stop_change_ratio) {
            cout << "[K-means] Change rate is below " << 100 * stop_change_ratio
                << "%, considered converged, stopping iteration early." << endl;
            break;
        }

        
        clusters = std::move(local_clusters);
        for (auto& group : cluster_nodes) group.clear();
        for (int i = 0; i < n; ++i) {
            if (clusters[i] >= 0 && clusters[i] < k_clusters) {
                cluster_nodes[clusters[i]].push_back(vertex_list[i]);
            }
        }

        
        auto center_start = high_resolution_clock::now();
        #pragma omp parallel for schedule(dynamic)
        for (int c = 0; c < k_clusters; ++c) {
            auto& group = cluster_nodes[c];
            if (group.empty()) continue;

            vector<pair<int, int>> node_degrees; 
            for (int node : group) {
                auto deg_iter = out_degree_map.find(node);
                int degree = (deg_iter != out_degree_map.end()) ? deg_iter->second : 0;
                node_degrees.emplace_back(-degree, node);
            }

            
            sort(node_degrees.begin(), node_degrees.end());
            vector<int> candidates;
            int take = min(x_candidates, (int)node_degrees.size());
            for (int i = 0; i < take; ++i) {
                candidates.push_back(node_degrees[i].second);
            }
            if (candidates.empty()) candidates.push_back(group[0]);

            
            vector<const vector<uint64_t>*> candidate_kmers;
            candidate_kmers.reserve(candidates.size());
            for (int cand_node : candidates) {
                candidate_kmers.push_back(&node_to_minhash_kmers.at(cand_node));
            }

            
            int best_node = candidates[0];
            double best_avg_sim = -1.0;
            for (size_t i = 0; i < candidates.size(); ++i) {
                int cand_node = candidates[i];
                const auto& cand_kmer = *candidate_kmers[i];
                
                double total_sim = 0.0;
                for (size_t j = 0; j < candidates.size(); ++j) {
                    if (i == j) continue; 
                    total_sim += jaccard_similarity_hash(cand_kmer, *candidate_kmers[j]);
                }
                
                
                double avg_sim = total_sim / (candidates.size() - 1);
                if (avg_sim > best_avg_sim) {
                    best_avg_sim = avg_sim;
                    best_node = cand_node;
                }
            }

            
            centers[c] = best_node;
        }

        auto center_end = high_resolution_clock::now();
        cout << "[K-means] Center update time: "
            << duration_cast<milliseconds>(center_end - center_start).count() << "ms" << endl;

        auto iter_end = high_resolution_clock::now();
        cout << "[K-means] Iteration " << iter + 1 << " total time: "
            << duration_cast<milliseconds>(iter_end - iter_start).count() << "ms" << endl;
    }

    auto total_end = high_resolution_clock::now();
    cout << "[K-means] K-means total time: " 
         << duration_cast<milliseconds>(total_end - total_start).count() << "ms" << endl;

    return {clusters, centers};
}

vector<pair<double, int>> sort_centers_by_similarity_with_score(
    const vector<int>& centers,
    const unordered_map<int, vector<uint64_t>>& node_to_minhash_kmers) {

    int n = centers.size();
    vector<bool> visited(n, false);

    vector<double> max_sim_to_MST(n, -1.0); 

    vector<int> order;
    order.reserve(n);

    int current = 0;
    visited[current] = true;
    order.push_back(current);

    for (int i = 0; i < n; ++i) {
        if (i == current) continue;
        double sim = jaccard_similarity_hash(
            node_to_minhash_kmers.at(centers[current]),
            node_to_minhash_kmers.at(centers[i])
        );
        max_sim_to_MST[i] = sim;
    }

    for (int step = 1; step < n; ++step) {
        int next = -1;
        double max_sim = -1.0;

        for (int i = 0; i < n; ++i) {
            if (visited[i]) continue;
            if (max_sim_to_MST[i] > max_sim) {
                max_sim = max_sim_to_MST[i];
                next = i;
            }
        }
        if (next == -1) break;

        visited[next] = true;
        order.push_back(next);
        for (int i = 0; i < n; ++i) {
            if (visited[i]) continue;
            double new_sim = jaccard_similarity_hash(
                node_to_minhash_kmers.at(centers[next]),
                node_to_minhash_kmers.at(centers[i])
            );
            max_sim_to_MST[i] = max(max_sim_to_MST[i], new_sim);
        }
    }

    vector<pair<double, int>> result;
    for (int i = 0; i < (int)order.size(); ++i) {
        result.emplace_back(i, order[i]);
    }
    return result;
}