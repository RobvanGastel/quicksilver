#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    return IA.size();
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortPairs(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {
    uint32_t sum = 0;
    for (int i = 1; i < getNoEdges(); i++) {
        if (IA[i] == IA[i-1])
            sum++;
    }
    return sum;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

std::pair<uint32_t, uint32_t> SimpleGraph::SelectLabel(uint32_t label, bool reverse) {
    if (reverse) {
        return std::pair<uint32_t, uint32_t>(
            positions_adj_reverse[label][0],
            positions_adj_reverse[label][positions_adj_reverse[label].size()-1]
        );

    } else {
        return std::pair<uint32_t, uint32_t>(
            positions_adj[label][0],
            positions_adj[label][positions_adj[label].size()-1]
        );
    }
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    LabelCount.resize(L);
    LabelSource.resize(L);
    LabelTarget.resize(L);
    for(int i = 0; i < L; i++) {
        LabelCount[i] = 0;
        std::vector<uint32_t> zeroes(V, 0);
        LabelSource[i] = zeroes;
        LabelTarget[i] = zeroes;
    }
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    readInitialInfo(fileName);
    initialize_positions_adj();

    std::string line;
    std::ifstream graphFile { fileName };

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)"); // subject predicate object .
    std::getline(graphFile, line);
    std::smatch matches;

    // create positions_adj
    std::vector<std::vector<uint32_t>> adj = positions_adj;
    std::vector<std::vector<uint32_t>> adj_reverse = positions_adj_reverse;
    while(std::getline(graphFile, line)) {
        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            uint32_t i = adj[predicate][subject]++;
            uint32_t i_reverse = adj_reverse[predicate][object]++;

            IA[i] = object;
            IA_reverse[i] = subject;
        }
    }

    graphFile.close();

    // std::pair<uint32_t, uint32_t> x = SelectLabel(1, false);
    // std::cout << x.first << " " << x.second << std::endl;
}

void SimpleGraph::initialize_positions_adj() {
    positions_adj.resize(L);
    positions_adj_reverse.resize(L);

    for (uint32_t label = 0; label < L; label++) {
        positions_adj[label].resize(V+1);
        positions_adj_reverse[label].resize(V+1);
    }

    // add label posisitons to adj
    positions_adj[0][0] = 0;
    positions_adj_reverse[0][0] = 0;

    for (uint32_t label = 0; label < L; label++){
        uint32_t count = positions_adj[label][0] + LabelCount[label];

        if (label < L-1) {
            positions_adj[label+1][0] = count;
            positions_adj_reverse[label+1][0] = count;
        }

        positions_adj[label][V] = count;
        positions_adj_reverse[label][V] = count;
    }

    // add target positions to adj
    for (uint32_t label = 0; label < L; label++){
        uint32_t sourceIndex = positions_adj[label][0];
        uint32_t targetIndex = positions_adj[label][0];
        
        for (uint32_t i = 1; i < L; i++) {
            sourceIndex += LabelSource[label][i-1];
            targetIndex += LabelTarget[label][i-1];
            positions_adj[label][i] = sourceIndex;
            positions_adj_reverse[label][i] = targetIndex;
        }
    }
}

void SimpleGraph::readInitialInfo(const std::string &fileName) {
    std::string line;
    std::ifstream graphFile { fileName };

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)"); // subject predicate object .
    std::regex headerPat (R"((\d+),(\d+),(\d+))"); // noNodes,noEdges,noLabels

    // parse the header (1st line)
    std::getline(graphFile, line);
    std::smatch matches;
    if(std::regex_search(line, matches, headerPat)) {
        uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
        uint32_t noEdges = (uint32_t) std::stoul(matches[2]);
        uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

        setNoVertices(noNodes);
        setNoLabels(noLabels);
        IA.resize(noEdges);
        IA_reverse.resize(noEdges);
    } else {
        throw std::runtime_error(std::string("Invalid graph header!"));
    }

    // parse edge data
    while(std::getline(graphFile, line)) {
        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            LabelCount[predicate]++;
            LabelSource[predicate][subject]++;
            LabelTarget[predicate][object]++;
        }
    }

    graphFile.close();
}