#include <k2Tree.h>
#include <math.h>

k2Tree::k2Tree(uint32_t n)   {
    setNoLabels(1);
    setNoVertices(n);
}

uint8_t k2Tree::getK() {
    return k;
}

void k2Tree::setK(uint8_t k_) {
    k = k_;
}

uint8_t k2Tree::getH() {
    return h;
}

void k2Tree::setH(uint8_t h_) {
    h = h_;
}

uint32_t k2Tree::getNoLabels() const {
    return L;
}

void k2Tree::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    adj.resize(L);
    trees.resize(L);
    leaves.resize(L);
}

uint32_t k2Tree::getNoVertices() const {
    return V;
}

void k2Tree::setNoVertices(uint32_t n) {
    V = n;
    for (int i = 0; i < L; i++)
        adj[i].resize(V);
}

uint32_t k2Tree::getNoEdges() const {
    uint32_t sum = 0;
    for (const auto & l : adj)
        sum += l.size();
    return sum;
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortList(const uint32_t &a, const uint32_t &b) {
    if (a < b) return true;
    return false;
}

uint32_t k2Tree::getNoDistinctEdges() const {

    uint32_t sum = 0;

    for (auto label : adj) {
        for (auto sourceVec : label) {

            std::sort(sourceVec.begin(), sourceVec.end(), sortList);

            uint32_t prevTarget = 0;
            bool first = true;

            for (const auto &Tgt : sourceVec) {
                if (first || prevTarget != Tgt) {
                    first = false;
                    sum++;
                    prevTarget = Tgt;
                }
            }
        }
    }

    return sum;
}

void k2Tree::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                 "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                 std::to_string(edgeLabel) + ")");
    adj[edgeLabel][from].emplace_back(to);
}

bool k2Tree::edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    auto it = std::find(adj[edgeLabel][from].begin(), adj[edgeLabel][from].end(), to);
    return (it != adj[edgeLabel][from].end());
}

void k2Tree::readFromContiguousFile(const std::string &fileName) {

    std::string line;
    std::ifstream graphFile { fileName };

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)"); // subject predicate object .
    std::regex headerPat (R"((\d+),(\d+),(\d+))"); // noNodes,noEdges,noLabels

    // parse the header (1st line)
    std::getline(graphFile, line);
    std::smatch matches;
    if(std::regex_search(line, matches, headerPat)) {
        uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
        uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

        setNoLabels(noLabels);
        setNoVertices(noNodes);
    } else {
        throw std::runtime_error(std::string("Invalid graph header!"));
    }

    // parse edge data
    while(std::getline(graphFile, line)) {

        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            addEdge(subject, object, predicate);
        }
    }
    graphFile.close();
//    setK((uint8_t)4);
    createTrees();
}

void k2Tree::createTrees() {
    uint8_t to = 0;
    for (int val = 1; val < V; val *= k, to++) {}
    if (to > 1)
        setH(to);
    int n = pow(k, h);
    int g = 0;
    for (int label = 0; label < L; label++) {
        std::vector<std::vector<uint8_t>> T(h-1);
        std::vector<typename std::vector<uint32_t>::iterator> cursors;
        for (auto iter = adj[0].begin(); iter != adj[0].end(); iter++)
            cursors.emplace_back(iter->begin());
        build(adj[0], cursors, T, n, 1, 0, 0);


        adj.erase(adj.begin());
        std::cout << "finished " << label << '\n';
    }

}

uint8_t k2Tree::build (std::vector<std::vector<uint32_t>> adj, std::vector<typename std::vector<uint32_t>::iterator> &cursors, std::vector<std::vector<uint8_t>> &T, uint32_t n, uint32_t l, uint32_t p, uint32_t q) {
    std::vector<uint8_t> C;
    if (l == h) {
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (((p + i) < adj.size()) && (cursors[p + i] != adj[p + i].end()) &&
                    ((q + j) == *(cursors[p + i])))
                    C.push_back(1);
                else
                    C.push_back(0);
                if (C.back())
                    cursors[p + i]++;
            }
        }
    }
    else {
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                C.push_back(build(adj, cursors, T, n / k, l + 1, p + i * (n / k), q + j * (n / k)));
            }
        }
    }
    for (int i = 0; i < C.size(); i++){
        if (C[i] == 0)
            return 0;
    }
    T[l-1].insert(T[l-1].end(), C.begin(), C.end());
    return 1;
}