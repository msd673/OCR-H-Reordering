#include "eval.h"

int main(int argc, char* argv[]) {
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <graph_file> <order_files_dir> <is_directed(0/1)> <traversal_id> [repeat_count]\n";
        cerr << "traversal_id: 1-BFS 2-DFS 3-SSSP 4-WCC 5-PageRank 6-Diameter\n";
        return 1;
    }

    const string graph_file = argv[1];
    const string order_files_dir = argv[2];
    bool is_directed = (stoi(argv[3]) != 0);
    int traversal_id = stoi(argv[4]);
    int repeat_count = (argc >= 6) ? stoi(argv[5]) : 5;
    if (repeat_count <= 0) {
        cerr << "repeat_count must be positive." << endl;
        return 1;
    }

    cout << "is_directed: " << boolalpha << is_directed << endl;
    cout << "traversal_type: " << traversal_id << " (1-BFS 2-DFS 3-SSSP 4-WCC 5-PageRank 6-Diameter)" << endl;

    vector<int> csr_offsets, csr_edges;
    read_csr_direct(graph_file, csr_offsets, csr_edges);
    CSRGraph original_graph(csr_offsets, csr_edges);
    cout << "-----0-----" << endl;

    vector<int> test_nodes;
    if (traversal_id >= 1 && traversal_id <= 3) {
        test_nodes = get_random_vertices(original_graph, 40);
    }

    vector<pair<string, vector<int>>> all_orders = {
        {"original", {}},
        {"degree_sort", sort_nodes_by_degree(original_graph)}
    };

    for (const auto& file : get_all_files(order_files_dir)) {
        size_t pos = file.find_last_of('/');
        string name = (pos != string::npos) ? file.substr(pos + 1) : file;
        all_orders.emplace_back(name, read_vertices_from_file(file));
    }

    vector<pair<string, CSRGraph>> prepared_graphs;
    prepared_graphs.reserve(all_orders.size());
    
    for (const auto& [name, order] : all_orders) {
        if (name == "original") {
            prepared_graphs.emplace_back(name, original_graph);
        } else {
            auto reordered_pair = reorder_by_vector(original_graph, order);
            CSRGraph& reordered_graph = reordered_pair.first;
            prepared_graphs.emplace_back(name, move(reordered_graph));
        }
    }

    vector<int> optimal_order;
    CSRGraph optimal_reordered_graph;

    
    if (traversal_id >= 4 && traversal_id <= 6) {
        if (traversal_id == 4) {
            
            CSRGraph re_graph = build_reverse_graph(original_graph);
            optimal_order = compute_wcc_bfs_only(original_graph, re_graph);
        }
        else if (traversal_id == 5) {
            
            optimal_order = compute_pagerank_order(original_graph);
        }
        else if (traversal_id == 6) {
            
            optimal_order = bfs_order_csr(original_graph, 0);
        }

        
        auto opt_pair = reorder_by_vector(original_graph, optimal_order);
        optimal_reordered_graph = move(opt_pair.first);
    }

    cout << "-----1-----" << endl;
    cout << "-----2-----" << endl;

    cout << "╔═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗\n";
    cout << "║          Cache & traversal Performance Analysis (Hardware Measured Possible)                                ║\n";
    cout << "╚═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝\n\n";


    if (traversal_id == 1 || traversal_id == 2 || traversal_id == 3) {
        for (int node : test_nodes) {
            cout << "┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
            cout << "│ Testing Node (Original ID): " << setw(24) << left << node << "                                                                   │\n";
            cout << "└───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘\n";

            for (auto& [name, cur] : prepared_graphs) {
                vector<CachePerfResult> cache_results;
                for (int i = 0; i < repeat_count; ++i) {
                    if (traversal_id == 1)
                        cache_results.push_back(measure_cache_perf([&]() { bfs_order_csr(cur, node); }));
                    else if (traversal_id == 2)
                        cache_results.push_back(measure_cache_perf([&]() { dfs_order_csr(cur, node); }));
                    else if (traversal_id == 3)
                        cache_results.push_back(measure_cache_perf([&]() { sssp_order_csr(cur.row_ptr, cur.col_idx, node); }));
                }

                double avg_l1_ref = 0, avg_l1_miss = 0;
                double avg_llc_ref = 0, avg_llc_miss = 0;
                double avg_cache_miss_total = 0;
                for (const auto& r : cache_results) {
                    avg_l1_ref += r.l1_ref;
                    avg_l1_miss += r.l1_miss;
                    avg_llc_ref += r.llc_ref;
                    avg_llc_miss += r.llc_miss;
                    avg_cache_miss_total += r.llc_miss;
                }
                avg_l1_ref /= repeat_count;
                avg_l1_miss /= repeat_count;
                avg_llc_ref /= repeat_count;
                avg_llc_miss /= repeat_count;
                avg_cache_miss_total /= repeat_count;

                double l1_mr = (avg_l1_ref > 0) ? 100.0 * avg_l1_miss / avg_l1_ref : 0.0;
                double l3_r = (avg_l1_ref > 0) ? 100.0 * avg_llc_ref / avg_l1_ref : 0.0;
                double cache_mr = (avg_l1_ref > 0) ? 100.0 * avg_cache_miss_total / avg_l1_ref : 0.0;

                double avg_time = 0;
                if (traversal_id == 1)
                    avg_time = measure_avg_time(repeat_count, [&](){ bfs_order_csr(cur, node); });
                else if (traversal_id == 2)
                    avg_time = measure_avg_time(repeat_count, [&](){ dfs_order_csr(cur, node); });
                else if (traversal_id == 3)
                    avg_time = measure_avg_time(repeat_count, [&](){ sssp_order_csr(cur.row_ptr, cur.col_idx, node); });

                cout << "   ├─ " << setw(15) << left << name << " "
                     << "L1-ref: " << setw(12) << right << static_cast<long long>(avg_l1_ref) << " │ "
                     << "L1-mr: " << fixed << setprecision(1) << setw(6) << right << l1_mr << "% │ "
                     << "L3-ref: " << setw(12) << right << static_cast<long long>(avg_llc_ref) << " │ "
                     << "L3-r: " << fixed << setprecision(1) << setw(6) << right << l3_r << "% │ "
                     << "Cache-mr: " << fixed << setprecision(1) << setw(6) << right << cache_mr << "% │ "
                     << "Time: " << fixed << setprecision(2) << setw(9) << right << avg_time << " ms\n";
            }

            vector<int> optimal_order;
            string optimal_name;
            CSRGraph optimal_graph(0);
            int real_opt_node = node;

            if (traversal_id == 1) {
                optimal_order = bfs_order_csr(original_graph, node);
                optimal_name = "Optimal-BFS";
            }
            // 2 = DFS
            else if (traversal_id == 2) {
                optimal_order = dfs_order_csr(original_graph, node);
                optimal_name = "Optimal-DFS";
            }
            else if (traversal_id == 3) {
                optimal_order = sssp_order_csr(original_graph.row_ptr, original_graph.col_idx, node);
                optimal_name = "Optimal-SSSP";
            }

            auto [opt_g, opt_map] = reorder_by_vector(original_graph, optimal_order);
            optimal_graph = move(opt_g);
            real_opt_node = opt_map[node];

            vector<CachePerfResult> cache_opt;
            for (int i = 0; i < repeat_count; ++i) {
                if (traversal_id == 1)
                    cache_opt.push_back(measure_cache_perf([&]() { bfs_order_csr(optimal_graph, real_opt_node); }));
                else if (traversal_id == 2)
                    cache_opt.push_back(measure_cache_perf([&]() { dfs_order_csr(optimal_graph, real_opt_node); }));
                else if (traversal_id == 3)
                    cache_opt.push_back(measure_cache_perf([&]() { 
                        sssp_order_csr(optimal_graph.row_ptr, optimal_graph.col_idx, real_opt_node); 
                    }));
            }

            double l1r = 0, l1m = 0, llcr = 0, llcm = 0;
            for (auto& r : cache_opt) {
                l1r += r.l1_ref;
                l1m += r.l1_miss;
                llcr += r.llc_ref;
                llcm += r.llc_miss;
            }
            l1r /= repeat_count;
            l1m /= repeat_count;
            llcr /= repeat_count;
            llcm /= repeat_count;

            double l1mr = (l1r > 0) ? 100.0 * l1m / l1r : 0.0;
            double l3r  = (l1r > 0) ? 100.0 * llcr / l1r : 0.0;
            double cmr  = (l1r > 0) ? 100.0 * llcm / l1r : 0.0;

            double time_opt = 0;
            if (traversal_id == 1)
                time_opt = measure_avg_time(repeat_count, [&]() { bfs_order_csr(optimal_graph, real_opt_node); });
            else if (traversal_id == 2)
                time_opt = measure_avg_time(repeat_count, [&]() { dfs_order_csr(optimal_graph, real_opt_node); });
            else if (traversal_id == 3)
                time_opt = measure_avg_time(repeat_count, [&]() { 
                    sssp_order_csr(optimal_graph.row_ptr, optimal_graph.col_idx, real_opt_node); 
                });

            cout << "   ├─ " << setw(15) << left << optimal_name << " "
                << "L1-ref: " << setw(12) << right << static_cast<long long>(l1r) << " │ "
                << "L1-mr: " << fixed << setprecision(1) << setw(6) << right << l1mr << "% │ "
                << "L3-ref: " << setw(12) << right << static_cast<long long>(llcr) << " │ "
                << "L3-r: " << fixed << setprecision(1) << setw(6) << right << l3r << "% │ "
                << "Cache-mr: " << fixed << setprecision(1) << setw(6) << right << cmr << "% │ "
                << "Time: " << fixed << setprecision(2) << setw(9) << right << time_opt << " ms\n";
            cout <<endl;    
        }
    }

    else if (traversal_id >= 4 && traversal_id <= 6) {
        vector<int> component_counts;
        if (traversal_id == 4) {
            for (auto& [name, cur] : prepared_graphs) {
                vector<int> res = compute_wcc_unionfind(cur);
                int maxc = *max_element(res.begin(), res.end());
                component_counts.push_back(maxc + 1);
            }
            vector<int> res_opt = compute_wcc_unionfind(optimal_reordered_graph);
            int maxc_opt = *max_element(res_opt.begin(), res_opt.end());
            component_counts.push_back(maxc_opt + 1);
        }

        for (int iter = 0; iter < 10; ++iter) {
            cout << "┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
            cout << "│ Testing Iter: " << setw(24) << left << iter << "                                                                   │\n";
            cout << "└───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘\n";

            for (size_t idx = 0; idx < prepared_graphs.size(); ++idx) {
                auto& [name, cur] = prepared_graphs[idx];
                vector<CachePerfResult> cache_results;
                for (int i = 0; i < repeat_count; ++i) {
                    if (traversal_id == 4) {
                        cache_results.push_back(measure_cache_perf([&]() {
                            for (int k = 0; k < 10; k++) compute_wcc_unionfind(cur);
                        }));
                    } else if (traversal_id == 5) {
                        cache_results.push_back(measure_cache_perf([&]() { compute_pagerank(cur); }));
                    } else if (traversal_id == 6) {
                        vector<int> s = get_top_degree_vertices(cur, 50);
                        cache_results.push_back(measure_cache_perf([&]() { compute_diameter_with_sources_parallel(cur, s); }));
                    }
                }

                double avg_l1_ref = 0, avg_l1_miss = 0;
                double avg_llc_ref = 0, avg_llc_miss = 0;
                double avg_cache_miss_total = 0;
                for (const auto& r : cache_results) {
                    avg_l1_ref += r.l1_ref;
                    avg_l1_miss += r.l1_miss;
                    avg_llc_ref += r.llc_ref;
                    avg_llc_miss += r.llc_miss;
                    avg_cache_miss_total += r.llc_miss;
                }
                avg_l1_ref /= repeat_count;
                avg_l1_miss /= repeat_count;
                avg_llc_ref /= repeat_count;
                avg_llc_miss /= repeat_count;
                avg_cache_miss_total /= repeat_count;

                double l1_mr = (avg_l1_ref > 0) ? 100.0 * avg_l1_miss / avg_l1_ref : 0.0;
                double l3_r = (avg_l1_ref > 0) ? 100.0 * avg_llc_ref / avg_l1_ref : 0.0;
                double cache_mr = (avg_l1_ref > 0) ? 100.0 * avg_cache_miss_total / avg_l1_ref : 0.0;

                double avg_time = 0;
                if (traversal_id == 4) {
                    avg_time = measure_avg_time(repeat_count, [&](){
                        for (int k = 0; k < 10; k++) compute_wcc_unionfind(cur);
                    }) / 10.0;
                } else if (traversal_id == 5) {
                    avg_time = measure_avg_time(repeat_count, [&](){ compute_pagerank(cur); });
                } else {
                    vector<int> s = get_top_degree_vertices(cur, 50);
                    avg_time = measure_avg_time(repeat_count, [&](){ compute_diameter_with_sources_parallel(cur, s); });
                }

                if (traversal_id == 4) {
                    cout << "[WCC-UnionFind] Total weakly connected components:  " << component_counts[idx] << endl;
                }

                cout << "   ├─ " << setw(15) << left << name << " "
                     << "L1-ref: " << setw(12) << right << static_cast<long long>(avg_l1_ref) << " │ "
                     << "L1-mr: " << fixed << setprecision(1) << setw(6) << right << l1_mr << "% │ "
                     << "L3-ref: " << setw(12) << right << static_cast<long long>(avg_llc_ref) << " │ "
                     << "L3-r: " << fixed << setprecision(1) << setw(6) << right << l3_r << "% │ "
                     << "Cache-mr: " << fixed << setprecision(1) << setw(6) << right << cache_mr << "% │ "
                     << "Time: " << fixed << setprecision(2) << setw(9) << right << avg_time << " ms\n";
            }

            vector<CachePerfResult> cache_results_optimal;
                for (int i = 0; i < repeat_count; ++i) {
                    if (traversal_id == 4) {
                        cache_results_optimal.push_back(measure_cache_perf([&]() {
                            for (int k = 0; k < 10; k++) compute_wcc_unionfind(optimal_reordered_graph);
                        }));
                    } else if (traversal_id == 5) {
                        cache_results_optimal.push_back(measure_cache_perf([&]() { 
                            compute_pagerank(optimal_reordered_graph); 
                        }));
                    } else if (traversal_id == 6) {
                        vector<int> s = get_top_degree_vertices(optimal_reordered_graph, 50);
                        cache_results_optimal.push_back(measure_cache_perf([&]() { 
                            compute_diameter_with_sources_parallel(optimal_reordered_graph, s); 
                        }));
                    }
                }

                double avg_l1_ref_opt = 0, avg_l1_miss_opt = 0;
                double avg_llc_ref_opt = 0, avg_llc_miss_opt = 0;
                double avg_cache_miss_total_opt = 0;
                for (const auto& r : cache_results_optimal) {
                    avg_l1_ref_opt += r.l1_ref;
                    avg_l1_miss_opt += r.l1_miss;
                    avg_llc_ref_opt += r.llc_ref;
                    avg_llc_miss_opt += r.llc_miss;
                    avg_cache_miss_total_opt += r.llc_miss;
                }
                avg_l1_ref_opt /= repeat_count;
                avg_l1_miss_opt /= repeat_count;
                avg_llc_ref_opt /= repeat_count;
                avg_llc_miss_opt /= repeat_count;
                avg_cache_miss_total_opt /= repeat_count;

                double l1_mr_opt = (avg_l1_ref_opt > 0) ? 100.0 * avg_l1_miss_opt / avg_l1_ref_opt : 0.0;
                double l3_r_opt = (avg_l1_ref_opt > 0) ? 100.0 * avg_llc_ref_opt / avg_l1_ref_opt : 0.0;
                double cache_mr_opt = (avg_l1_ref_opt > 0) ? 100.0 * avg_cache_miss_total_opt / avg_l1_ref_opt : 0.0;

                double t_opt = 0;
                if (traversal_id == 4) {
                    t_opt = measure_avg_time(repeat_count, [&](){
                        for (int k = 0; k < 10; k++) compute_wcc_unionfind(optimal_reordered_graph);
                    }) / 10.0;
                } else if (traversal_id == 5) {
                    t_opt = measure_avg_time(repeat_count, [&](){ 
                        compute_pagerank(optimal_reordered_graph); 
                    });
                } else {
                    vector<int> s = get_top_degree_vertices(optimal_reordered_graph, 50);
                    t_opt = measure_avg_time(repeat_count, [&](){ 
                        compute_diameter_with_sources_parallel(optimal_reordered_graph, s); 
                    });
                }

                if (traversal_id == 4) {
                    cout << "[WCC-UnionFind] Total weakly connected components:  " << component_counts.back() << endl;
                }

                cout << "   ├─ " << setw(15) << left << "Self Order" << " "
                    << "L1-ref: " << setw(12) << right << static_cast<long long>(avg_l1_ref_opt) << " │ "
                    << "L1-mr: " << fixed << setprecision(1) << setw(6) << right << l1_mr_opt << "% │ "
                    << "L3-ref: " << setw(12) << right << static_cast<long long>(avg_llc_ref_opt) << " │ "
                    << "L3-r: " << fixed << setprecision(1) << setw(6) << right << l3_r_opt << "% │ "
                    << "Cache-mr: " << fixed << setprecision(1) << setw(6) << right << cache_mr_opt << "% │ "
                    << "Time: " << fixed << setprecision(2) << setw(9) << right << t_opt << " ms\n";

                cout << endl;
            }
    } else {
        cerr << "Invalid traversal number!\n";
        return 1;
    }

    return 0;
}