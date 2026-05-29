#include "eval.h"

vector<int> bfs_order_csr(const CSRGraph& graph, int start) {
    const int n = graph.num_nodes;
    vector<int> order;
    order.reserve(n);
    vector<char> visited(n, 0);
    queue<int> q;

    visited[start] = 1;
    q.push(start);

    while (!q.empty()) {
        int current = q.front();
        q.pop();
        order.push_back(current);

        int start_idx = graph.row_ptr[current];
        int end_idx = graph.row_ptr[current + 1];
        for (int i = start_idx; i < end_idx; ++i) {
            int neighbor = graph.col_idx[i];
            if (!visited[neighbor]) {
                visited[neighbor] = 1;
                q.push(neighbor);
            }
        }
    }
    return order;
}
vector<int> dfs_order_csr(const CSRGraph& graph, int start) {
    const int n = graph.num_nodes;

    vector<int> order;
    order.reserve(n);

    vector<char> visited(n, 0);
    stack<int> st;  
    
    visited[start] = 1;
    st.push(start);

    while (!st.empty()) {
        int current = st.top();
        st.pop();

        order.push_back(current);
        
        int start_idx = graph.row_ptr[current];
        int end_idx = graph.row_ptr[current + 1];
        
        
        for (int i = end_idx - 1; i >= start_idx; --i) {
            int neighbor = graph.col_idx[i];
            
            if (!visited[neighbor]) {
                visited[neighbor] = 1;
                st.push(neighbor);
            }
        }
    }
    return order;
}

vector<int> sssp_order_csr(const vector<int>& csr_offsets,
                          const vector<int>& csr_edges,
                          int start_node) {
    int n = csr_offsets.size() - 1;
    vector<int> order;
    order.reserve(n);
    vector<char> visited(n, 0);
    vector<int> dist(n, INT_MAX);
    queue<int> q;

    dist[start_node] = 0;
    q.push(start_node);
    visited[start_node] = 1;
    order.push_back(start_node);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int i = csr_offsets[u]; i < csr_offsets[u + 1]; ++i) {
            int v = csr_edges[i];
            if (!visited[v]) {
                visited[v] = 1;
                dist[v] = dist[u] + 1;
                order.push_back(v);
                q.push(v);
            }
        }
    }
    return order;
}

vector<int> compute_wcc_unionfind(const CSRGraph& graph)
{
    const int V = graph.num_nodes;
    UnionFind uf(V);

    for (int u = 0; u < V; ++u) {
        for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
            int v = graph.col_idx[i];
            uf.unite(u, v);
        }
    }

    
    vector<int> comp_id(V, -1);
    vector<int> root_to_comp(V, -1);

    int comp_cnt = 0;
    for (int v = 0; v < V; ++v) {
        int root = uf.find(v);
        if (root_to_comp[root] == -1) {
            root_to_comp[root] = comp_cnt++;
        }
        comp_id[v] = root_to_comp[root];
    }
  
    return comp_id;
}

CSRGraph build_reverse_graph(const CSRGraph& graph)
{
    const vector<int>& row_ptr = graph.row_ptr;
    const vector<int>& col_idx = graph.col_idx;
    const int V = graph.num_nodes;
    const int E = static_cast<int>(col_idx.size());

    vector<int> in_deg(V, 0);

    for (int u = 0; u < V; ++u) {
        for (int i = row_ptr[u]; i < row_ptr[u + 1]; ++i) {
            int v = col_idx[i];
            in_deg[v]++;
        }
    }

    vector<int> rev_row_ptr(V + 1, 0);
    for (int v = 0; v < V; ++v) {
        rev_row_ptr[v + 1] = rev_row_ptr[v] + in_deg[v];
    }

    vector<int> rev_col_idx(E);
    vector<int> cur = rev_row_ptr;

    
    for (int u = 0; u < V; ++u) {
        for (int i = row_ptr[u]; i < row_ptr[u + 1]; ++i) {
            int v = col_idx[i];
            rev_col_idx[cur[v]++] = u;
        }
    }

    CSRGraph rev_graph(V);
    rev_graph.row_ptr = std::move(rev_row_ptr);
    rev_graph.col_idx = std::move(rev_col_idx);

    return rev_graph;
}

vector<int> compute_wcc_bfs_only(const CSRGraph& graph, const CSRGraph& rev_graph)
{
    const int V = graph.num_nodes;

    vector<int> comp_id(V, -1);
    vector<int> order;
    order.reserve(V);

    int comp_cnt = 0;
    queue<int> q;

    for (int s = 0; s < V; ++s) {
        if (comp_id[s] != -1) continue;

        comp_id[s] = comp_cnt;
        q.push(s);
        order.push_back(s);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
                int v = graph.col_idx[i];
                if (comp_id[v] == -1) {
                    comp_id[v] = comp_cnt;
                    q.push(v);
                    order.push_back(v);
                }
            }

            for (int i = rev_graph.row_ptr[u]; i < rev_graph.row_ptr[u + 1]; ++i) {
                int v = rev_graph.col_idx[i];
                if (comp_id[v] == -1) {
                    comp_id[v] = comp_cnt;
                    q.push(v);
                    order.push_back(v);
                }
            }
        }

        comp_cnt++;
    }

    return order;
}

vector<int> compute_wcc(const CSRGraph& graph)
{
    const vector<int>& row_ptr = graph.row_ptr;
    const vector<int>& col_idx = graph.col_idx;
    const int V = row_ptr.size() - 1;

    vector<int> deg(V, 0);
    for (int u = 0; u < V; ++u) {
        for (int i = row_ptr[u]; i < row_ptr[u + 1]; ++i) {
            int v = col_idx[i];
            deg[u]++;
            deg[v]++;
        }
    }

    vector<int> und_row_ptr(V + 1, 0);
    for (int i = 0; i < V; ++i) {
        und_row_ptr[i + 1] = und_row_ptr[i] + deg[i];
    }

    vector<int> und_col_idx(und_row_ptr[V]);
    vector<int> cur = und_row_ptr;

    for (int u = 0; u < V; ++u) {
        for (int i = row_ptr[u]; i < row_ptr[u + 1]; ++i) {
            int v = col_idx[i];
            und_col_idx[cur[u]++] = v; 
            und_col_idx[cur[v]++] = u; 
        }
    }

    vector<int> comp_id(V, -1);
    vector<int> order;
    order.reserve(V);

    int comp_cnt = 0;
    queue<int> q;

    for (int s = 0; s < V; ++s) {
        if (comp_id[s] != -1) continue;

        comp_id[s] = comp_cnt;
        q.push(s);
        order.push_back(s);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int i = und_row_ptr[u]; i < und_row_ptr[u + 1]; ++i) {
                int v = und_col_idx[i];
                if (comp_id[v] == -1) {
                    comp_id[v] = comp_cnt;
                    q.push(v);
                    order.push_back(v);
                }
            }
        }

        comp_cnt++;
    }

    return order;
}


vector<double> compute_pagerank(const CSRGraph& graph)
{
    const vector<int>& csr_offsets = graph.row_ptr;
    const vector<int>& csr_edges   = graph.col_idx;

    const int V = csr_offsets.size() - 1;
    
    vector<double> pr(V, 1.0 / V);
    vector<double> new_pr(V, 0.0);
    double base = (1.0 - DAMPING_FACTOR) / V;

    for (int iter = 0; iter < MAX_PR_ITER; ++iter) {
        fill(new_pr.begin(), new_pr.end(), base);

        for (int u = 0; u < V; ++u) {
            int out_degree = csr_offsets[u + 1] - csr_offsets[u];
            if(out_degree==0) continue;
            double val = DAMPING_FACTOR * pr[u] / out_degree;

            for (int i = csr_offsets[u]; i < csr_offsets[u+1]; ++i) {
                int v = csr_edges[i];
                new_pr[v] += val;
            }
        }

        double delta = 0.0;
        for (int u = 0; u < V; ++u) {
            delta += fabs(new_pr[u] - pr[u]);
        }
        pr.swap(new_pr);

        if (delta < PR_EPS) {
            break;
        }
    }
    return pr;
}

vector<int> compute_pagerank_order(const CSRGraph& graph)
{
    const vector<int>& csr_offsets = graph.row_ptr;
    const vector<int>& csr_edges   = graph.col_idx;

    const int V = csr_offsets.size() - 1;
    
    vector<double> pr(V, 1.0 / V);
    vector<double> new_pr(V, 0.0);
    double base = (1.0 - DAMPING_FACTOR) / V;
    
    vector<int> order;
    order.reserve(V);
    vector<bool> visited(V, false);

    fill(new_pr.begin(), new_pr.end(), base);
    for (int u = 0; u < V; ++u) {
        if (!visited[u]) {
            visited[u] = true;
            order.push_back(u);
        }
        int out_degree = csr_offsets[u + 1] - csr_offsets[u];
        if(out_degree==0) continue;
        double val = pr[u] / out_degree;

        for (int i = csr_offsets[u]; i < csr_offsets[u+1]; ++i) {
            int v = csr_edges[i];
            if (!visited[v]) {
                visited[v] = true;
                order.push_back(v);
            }
            new_pr[v] += DAMPING_FACTOR * val;
        }
    }
    cout<<"Vertices: "<<V<<" order size: "<<order.size()<<endl;
    return order;
}


int bfs_eccentricity(const CSRGraph& graph, int src) {
    const vector<int>& row_ptr = graph.row_ptr;
    const vector<int>& col_idx = graph.col_idx;
    int n = graph.num_nodes;

    vector<int> dist(n, -1);
    queue<int> q;

    dist[src] = 0;
    q.push(src);

    int max_dist = 0;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int eid = row_ptr[u]; eid < row_ptr[u + 1]; ++eid) {
            int v = col_idx[eid];
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                max_dist = max(max_dist, dist[v]);
                q.push(v);
            }
        }
    }

    return max_dist;
}

int compute_diameter(const CSRGraph& graph,
    const vector<int>& sources) {
    int n = graph.num_nodes;
    if (n == 0) return 0;
    if (n == 1) return 0;

    int approx_diameter = 0;

    for (int i = 0; i < (int)sources.size(); ++i) {
        int src = sources[i];
        int ecc = bfs_eccentricity(graph, src);
        approx_diameter = max(approx_diameter, ecc);
    }

    return approx_diameter;
}

static int bfs_eccentricity_reuse(
    const CSRGraph& graph,
    int src,
    vector<int>& dist,
    vector<int>& visited,
    vector<int>& q)
{
    int head = 0, tail = 0;
    q[tail++] = src;
    dist[src] = 0;
    visited.push_back(src);

    int max_dist = 0;

    while (head < tail) {
        int u = q[head++];

        for (int eid = graph.row_ptr[u]; eid < graph.row_ptr[u + 1]; ++eid) {
            int v = graph.col_idx[eid];
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                if (dist[v] > max_dist) max_dist = dist[v];
                q[tail++] = v;
                visited.push_back(v);
            }
        }
    }

    for (int v : visited) dist[v] = -1;
    visited.clear();

    return max_dist;
}

int compute_diameter_with_sources_parallel(
    const CSRGraph& graph,
    const vector<int>& sources)
{
    int n = graph.num_nodes;
    if (n == 0 || n == 1 || sources.empty()) return 0;

    int approx_diameter = 0;

    #pragma omp parallel
    {
        vector<int> dist(n, -1);
        vector<int> visited;
        visited.reserve(min(n, 1 << 20));
        vector<int> q(n);

        int local_max = 0;

        #pragma omp for schedule(dynamic, 1) nowait
        for (int i = 0; i < (int)sources.size(); ++i) {
            int src = sources[i];
            int ecc = bfs_eccentricity_reuse(graph, src, dist, visited, q);
            if (ecc > local_max) local_max = ecc;
        }

        #pragma omp critical
        {
            if (local_max > approx_diameter) approx_diameter = local_max;
        }
    }

    return approx_diameter;
}

vector<int> sort_nodes_by_degree(const CSRGraph& graph) {
    vector<pair<int, int>> degree_list;
    for (int u = 0; u < graph.num_nodes; ++u) {
        int degree = graph.row_ptr[u + 1] - graph.row_ptr[u];
        degree_list.emplace_back(u, degree);
    }
    sort(degree_list.begin(), degree_list.end(), [](auto& a, auto& b) {
        return a.second < b.second;
    });

    vector<int> sorted;
    for (auto& p : degree_list) sorted.push_back(p.first);
    return sorted;
}


vector<int> get_top_degree_vertices(const CSRGraph& graph, int count) {
    vector<pair<int, int>> deg;
    for (int u = 0; u < graph.num_nodes; ++u) {
        int degree = graph.row_ptr[u + 1] - graph.row_ptr[u];
        deg.emplace_back(u, degree);
    }
    sort(deg.begin(), deg.end(), [](auto& a, auto& b) {
        return a.second > b.second;
    });

    vector<int> top;
    for (int i = 0; i < count && i < (int)deg.size(); ++i)
        top.push_back(deg[i].first);
    return top;
}

vector<int> read_vertices_from_file(const string& filename) {
    ifstream infile(filename);
    vector<int> vertices;
    string line;
    int node;
    while (getline(infile, line)) {
        istringstream iss(line);
        if (iss >> node) vertices.push_back(node);
    }
    return vertices;
}

pair<CSRGraph, unordered_map<int, int>> reorder_by_vector(const CSRGraph& graph, const vector<int>& order) {
    const int num_nodes = graph.num_nodes;
    unordered_map<int, int> old_to_new;

    for (int i = 0; i < (int)order.size(); ++i) {
        old_to_new[order[i]] = i;
    }
    int next_id = order.size();
    for (int v = 0; v < num_nodes; ++v) {
        if (old_to_new.find(v) == old_to_new.end()) {
            old_to_new[v] = next_id++;
        }
    }

    vector<int> new_degree(num_nodes, 0);
    for (int u = 0; u < num_nodes; ++u) {
        new_degree[old_to_new[u]] = graph.row_ptr[u + 1] - graph.row_ptr[u];
    }

    CSRGraph new_graph(num_nodes);
    new_graph.row_ptr[0] = 0;
    for (int u = 0; u < num_nodes; ++u) {
        new_graph.row_ptr[u + 1] = new_graph.row_ptr[u] + new_degree[u];
    }
    new_graph.col_idx.resize(new_graph.row_ptr.back());

    vector<int> temp_pos(num_nodes);
    for (int u = 0; u < num_nodes; ++u) {
        temp_pos[u] = new_graph.row_ptr[u];
    }

    vector<vector<int>> temp_adj_lists(num_nodes);
    for (int u_old = 0; u_old < num_nodes; ++u_old) {
        int u_new = old_to_new[u_old];
        for (int i = graph.row_ptr[u_old]; i < graph.row_ptr[u_old + 1]; ++i) {
            int v_old = graph.col_idx[i];
            int v_new = old_to_new[v_old];
            temp_adj_lists[u_new].push_back(v_new);
        }
    }

    for (int u_new = 0; u_new < num_nodes; ++u_new) {
        sort(temp_adj_lists[u_new].begin(), temp_adj_lists[u_new].end());
        for (int v_new : temp_adj_lists[u_new]) {
            new_graph.col_idx[temp_pos[u_new]++] = v_new;
        }
    }

    return {new_graph, old_to_new};
}

