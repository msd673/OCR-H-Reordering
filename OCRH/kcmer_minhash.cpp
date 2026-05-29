#include "graph.h"

vector<pair<int, int>> generate_hash_params(int num_hashes) {
    vector<pair<int, int>> params;
    srand(time(0));
    for (int i = 0; i < num_hashes; ++i) {
        int a = (rand() % (BIG_PRIME-1)) + 1;
        int b = rand() % BIG_PRIME;
        params.emplace_back(a, b);
    }
    return params;
}


uint64_t compute_kmer_hash(const vector<int>& kmer_vec, int a, int b) {
    uint64_t hash_val = 0;
    for (int x : kmer_vec) {
        uint64_t h = (uint64_t(a) * x + b) % BIG_PRIME;
        hash_val = (hash_val * 31 + h) % BIG_PRIME; 
    }
    return hash_val;
}

unordered_map<int, vector<uint64_t>> compute_topk_minhash_kmers(
    const unordered_map<int, vector<uint64_t>>& node_to_kmers,
    const vector<pair<int, int>>& hash_params) {

    unordered_map<int, vector<uint64_t>> node_to_minhash_map;
    node_to_minhash_map.reserve(node_to_kmers.size());

    vector<int> node_ids;
    node_ids.reserve(node_to_kmers.size());
    for (const auto& pair : node_to_kmers) {
        node_ids.push_back(pair.first);
    }

    #pragma omp parallel
    {
        unordered_map<int, vector<uint64_t>> local_map;

        #pragma omp for nowait
        for (size_t idx = 0; idx < node_ids.size(); ++idx) {
            int node_id = node_ids[idx];
            const auto& kmers = node_to_kmers.at(node_id);
            vector<uint64_t> sketch;

            for (const auto& param : hash_params) {
                vector<pair<uint64_t, uint64_t>> hashed_kmers; 
                hashed_kmers.reserve(kmers.size());

                for (uint64_t kmer_hash : kmers) {
                    uint64_t minhash_val = (uint64_t(param.first) * kmer_hash + param.second) % BIG_PRIME;
                    hashed_kmers.emplace_back(minhash_val, kmer_hash);
                }

                if ((int)hashed_kmers.size() > TOP_K_HASHES) {
                    nth_element(
                        hashed_kmers.begin(),
                        hashed_kmers.begin() + TOP_K_HASHES,
                        hashed_kmers.end()
                    );
                    for (int i = 0; i < TOP_K_HASHES; ++i) {
                        sketch.push_back(hashed_kmers[i].second);
                    }
                } else {
                    for (const auto& [_, kmer_hash] : hashed_kmers) {
                        sketch.push_back(kmer_hash);
                    }
                }
            }
            sort(sketch.begin(), sketch.end());
            sketch.erase(unique(sketch.begin(), sketch.end()), sketch.end());
            local_map[node_id] = std::move(sketch);
        }

        #pragma omp critical
        {
            for (auto& [node_id, sketch] : local_map) {
                node_to_minhash_map[node_id] = std::move(sketch);
            }
        }
    }

    return node_to_minhash_map;
}
double jaccard_similarity_hash(const vector<uint64_t>& A,  
                         const vector<uint64_t>& B) {
    size_t i = 0, j = 0;
    size_t intersection = 0;
    
    while (i < A.size() && j < B.size()) {
        if (A[i] == B[j]) {
            intersection++;
            i++; j++;
        } 
        else if (A[i] < B[j]) i++;
        else j++;
    }
    
    size_t union_size = A.size() + B.size() - intersection;
    return union_size > 0 ? (double)intersection / union_size : 0.0;
}