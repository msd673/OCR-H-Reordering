#include "eval.h"

vector<string> get_all_files(const string& dir_path) {
    vector<string> files;
    DIR* dir;
    struct dirent* ent;
    
    if ((dir = opendir(dir_path.c_str())) != nullptr) {
        while ((ent = readdir(dir)) != nullptr) {
            string filename = ent->d_name;
            if (filename != "." && filename != "..") {
                files.push_back(dir_path + "/" + filename);
            }
        }
        closedir(dir);
    } else {
        cerr << "Failed to open directory: " << dir_path << endl;
    }
    
    return files;
}
vector<int> get_random_vertices(const CSRGraph& graph, int count) {
    vector<int> vertices;
    int n = graph.num_nodes;
    if (n <= 0 || count <= 0) return vertices;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, n - 1);

    unordered_set<int> s;
    while (s.size() < (size_t)min(count, n)) {
        s.insert(dist(gen));
    }
    return vector<int>(s.begin(), s.end());
}