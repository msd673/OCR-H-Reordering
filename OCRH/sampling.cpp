#include "graph.h"

vector<int> bfs_order_csr(const vector<int>& csr_offsets,
                          const vector<int>& csr_edges,
                          int start_node) {
    int n = static_cast<int>(csr_offsets.size()) - 1;
    vector<int> order;
    order.reserve(n);

    vector<char> visited(n, 0);
    queue<int> q;

    q.push(start_node);
    visited[start_node] = 1;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        order.push_back(u);
        for (int i = csr_offsets[u]; i < csr_offsets[u + 1]; ++i) {
            int v = csr_edges[i];
            if (!visited[v]) {
                visited[v] = 1;
                q.push(v);
            }
        }
    }
    return order;
}

vector<int> dfs_order_csr(const vector<int>& csr_offsets,
                          const vector<int>& csr_edges,
                          int start_node) {
    int n = csr_offsets.size() - 1;
    vector<int> order;
    order.reserve(n);

    vector<char> visited(n, 0);
    stack<int> st;

    st.push(start_node);
    visited[start_node] = 1;

    while (!st.empty()) {
        int u = st.top(); st.pop();
        order.push_back(u);
        for (int i = csr_offsets[u]; i < csr_offsets[u + 1]; ++i) {
            int v = csr_edges[i];
            if (!visited[v]) {
                visited[v] = 1;
                st.push(v);
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

vector<int> bfs_in_cluster_global(
    int start_node,
    int cluster_id,
    const vector<int>& clusters,
    const vector<int>& csr_offsets,
    const vector<int>& csr_edges) {
    
    int n = csr_offsets.size() - 1;
    vector<int> order;
    order.reserve(10000000);
    vector<char> visited(n, 0);
    for (int c = 0; c < n; ++c) {
        int u_start = (c + start_node) % n;
        queue<int> q;
        if (clusters[u_start] != cluster_id || visited[u_start]) {
            continue;
        }
        
        q.push(u_start);
        visited[u_start] = 1;
        
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (u < 0 || u >= n) continue;
            order.push_back(u);

            int start_idx = csr_offsets[u];
            int end_idx = csr_offsets[u + 1];

            if (end_idx > static_cast<int>(csr_edges.size())) {
                end_idx = static_cast<int>(csr_edges.size());
            }
            
            for (int i = start_idx; i < end_idx; ++i) {
                int v = csr_edges[i];
                
                if (v < 0 || v >= n || v >= static_cast<int>(clusters.size())) {
                    cout << "CRITICAL: v is out of bounds! v=" << v 
                        << ", visited.size=" << visited.size() << endl;
                    continue; 
                }

                if (!visited[v] && clusters[v] == cluster_id) {
                    visited[v] = 1;
                    q.push(v);
                }
            }
        }
    }
    
    return order;
}

vector<int> dfs_in_cluster_global(
    int start_node,
    int cluster_id,
    const vector<int>& clusters,
    const vector<int>& csr_offsets,
    const vector<int>& csr_edges) {
    
    int n = static_cast<int>(csr_offsets.size()) - 1;
    vector<int> order;
    order.reserve(10000000);
    vector<char> visited(n, 0);
    for (int c = 0; c < n; ++c) {
        int u_start = (c + start_node) % n;
        stack<int> st;
        if (clusters[u_start] != cluster_id || visited[u_start]) {
            continue;
        }
        
        st.push(u_start);
        visited[u_start] = 1;
        
        while (!st.empty()) {
            int u = st.top(); st.pop();
            if (u < 0 || u >= n) continue;
            order.push_back(u);

            int start_idx = csr_offsets[u];
            int end_idx = csr_offsets[u + 1];
            
            if (end_idx > static_cast<int>(csr_edges.size())) {
                end_idx = static_cast<int>(csr_edges.size());
            }
            
            for (int i = start_idx; i < end_idx; ++i) {
                int v = csr_edges[i];
                
                if (v < 0 || v >= n || v >= static_cast<int>(clusters.size())) {
                    cout << "CRITICAL: v is out of bounds! v=" << v 
                        << ", visited.size=" << visited.size() << endl;
                    continue; 
                }

                if (!visited[v] && clusters[v] == cluster_id) {
                    visited[v] = 1;
                    st.push(v);
                }
            }
        }
    }
    
    return order;
}

vector<int> sssp_in_cluster_global(
    int start_node,
    int cluster_id,
    const vector<int>& clusters,
    const vector<int>& csr_offsets,
    const vector<int>& csr_edges) {
    
    int n = csr_offsets.size() - 1;
    vector<int> order;
    vector<char> visited(n, 0);

    for (int c = 0; c < n; ++c) {
        int u_start = (c + start_node) % n;
        if (clusters[u_start] != cluster_id || visited[u_start]) continue;

        queue<int> q;
        q.push(u_start);
        visited[u_start] = 1;
        order.push_back(u_start);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int i = csr_offsets[u]; i < csr_offsets[u+1]; ++i) {
                int v = csr_edges[i];
                if (v < 0 || v >= n || clusters[v] != cluster_id) continue;
                if (!visited[v]) {
                    visited[v] = 1;
                    order.push_back(v);
                    q.push(v);
                }
            }
        }
    }
    return order;
}