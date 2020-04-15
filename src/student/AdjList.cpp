#include "AdjList.h"

AdjList::AdjList(uint32_t n)   {
    setNoLabels(1);
    setNoVertices(n);
}

uint32_t AdjList::getNoLabels() const {
    return L;
}

void AdjList::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    adj.resize(L);
}

uint32_t AdjList::getNoVertices() const {
    return V;
}

void AdjList::setNoVertices(uint32_t n) {
    V = n;
    for (int i = 0; i < L; i++)
        adj[i].resize(V);
}

uint32_t AdjList::getNoEdges() const {
    uint32_t sum = 0;
    for (const auto & l : adj)
        sum += l.size();
    return sum;
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortLists(const uint32_t &a, const uint32_t &b) {
    if (a < b) return true;
    return false;
}

uint32_t AdjList::getNoDistinctEdges() const {

    uint32_t sum = 0;

    for (auto label : adj) {
        for (auto sourceVec : label) {

            std::sort(sourceVec.begin(), sourceVec.end(), sortLists);

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

void AdjList::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                 "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                 std::to_string(edgeLabel) + ")");
    adj[edgeLabel][from].emplace_back(to);
}

bool AdjList::edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    auto it = std::find(adj[edgeLabel][from].begin(), adj[edgeLabel][from].end(), to);
    return (it != adj[edgeLabel][from].end());
}

void AdjList::readFromContiguousFile(const std::string &fileName) {

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
}
