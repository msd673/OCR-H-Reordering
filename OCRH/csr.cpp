#include "graph.h"

bool read_csr_direct(const string& filename,
                     vector<int>& csr_offsets,
                     vector<int>& csr_edges,
                     unordered_map<int, int>& out_degree_map) {
    
    ifstream file(filename, ios::binary);
    if (!file) {
        cerr <<"Failed to open CSR file: " << filename << endl;
        return false;
    }
    file.seekg(0, ios::end);
    streampos fileSize = file.tellg();
    file.seekg(0, ios::beg);
    
    unsigned int num_nodes, num_edges;
    file.read((char*)&num_nodes, sizeof(unsigned int));
    file.read((char*)&num_edges, sizeof(unsigned int));
    if (fileSize !=
        sizeof(uint32_t)*2 +
        sizeof(uint32_t)*num_nodes +
        sizeof(uint32_t)*num_edges) {
        cerr << "File size does not match header information!" << endl;
    }
    
    if (!file) {
        cerr << "Failed to read CSR header" << endl;
        return false;
    }

    cout << "Reading CSR: num_nodes=" << num_nodes << ", num_edges=" << num_edges << endl;

    try {
        csr_offsets.resize(num_nodes + 1);
        csr_edges.resize(num_edges);
    } catch (const bad_alloc& e) {
        cerr << "Memory allocation failed: " << e.what() << endl;
        cerr << "Required: " << ((num_nodes + 1 + num_edges) * sizeof(int) / 1024 / 1024) << " MB" << endl;
        return false;
    }
    
    file.read((char*)csr_offsets.data(), num_nodes * sizeof(int));

    csr_offsets[num_nodes] = num_edges;

    file.read((char*)csr_edges.data(), num_edges * sizeof(int));

    
    if (csr_offsets[0] != 0) {
        cerr << "warning: offsets[0] != 0" << endl;
    }
    
    if (csr_offsets[num_nodes] != num_edges) {
        cerr << "warning:  Last element of offsets != number of edges" << endl;
        csr_offsets.back() = (int)csr_edges.size();
    }
    
    cout<<"csr_offsets size: "<< csr_offsets.size()<<endl;
    cout<<"csr_edges size: "<< csr_edges.size()<<endl;

    out_degree_map.clear();
    out_degree_map.reserve(num_nodes);
    
    for (int i = 0; i < num_nodes; i++) {
        int degree = csr_offsets[i+1] - csr_offsets[i];
        if (degree > 0) {
            out_degree_map[i] = degree;
        }
    }
    
    cout << "Successfully read CSR file" << endl;
    cout << "Built out_degree_map with " << out_degree_map.size() << " nodes" << endl;
    
    return true;
}
