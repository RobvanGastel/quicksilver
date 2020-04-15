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

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    LabelCount.resize(L);
    LabelSource.resize(L);
    LabelTarget.resize(L);
    subjects.resize(L);
    objects.resize(L);
    for(int i = 0; i < L; i++) {
        LabelCount[i] = 0;
        std::vector<uint32_t> zeroes(V, 0);
        LabelSource[i] = zeroes;
        LabelTarget[i] = zeroes;
    }
}


bool SimpleGraph::edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel) {

    return false;
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    readInitialInfo(fileName);
    initialize_positions_adj();
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

        setNoVertices(noNodes);
        setNoLabels(noLabels);
    } else {
        throw std::runtime_error(std::string("Invalid graph header!"));
    }

    std::vector<std::vector<uint32_t>> offset;
    std::vector<std::vector<uint32_t>> rev_offset;
    offset.resize(L);
    rev_offset.resize(L);
    for (int i = 0; i < L; i++){
        std::vector<uint32_t> zeroes(subjects[i].size(), 0);
        offset[i] = zeroes;
    }
    for (int i = 0; i < L; i++){
        std::vector<uint32_t> zeroes(objects[i].size(), 0);
        rev_offset[i] = zeroes;
    }
    // parse edge data
    while(std::getline(graphFile, line)) {

        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            int index;
            for (index = 0; index < subjects[predicate].size(); index++){
                if (subjects[predicate][index] == subject) {
                    uint32_t i = offset[predicate][index];
                    offset[predicate][index]++;
                    IA[positions_adj[predicate][index] + i] = object;
                    break;
                }
            }
            for (index = 0; index < objects[predicate].size(); index++){
                if (objects[predicate][index] == object) {
                    uint32_t i_reverse = rev_offset[predicate][index];
                    rev_offset[predicate][index]++;
                    IA_reverse[positions_adj_reverse[predicate][index] + i_reverse] = subject;
                    break;
                }
            }
        }
    }
    std::vector<uint32_t > N = findNeighbours(8125, 1, true);
    for (int i = 0; i < N.size(); i++)
        std::cout << N[i] << "\n";
    graphFile.close();
}

void SimpleGraph::initialize_positions_adj() {
    positions_adj.resize(L);
    positions_adj_reverse.resize(L);

    for (uint32_t label = 0; label < L; label++) {
        positions_adj[label].resize(subjects[label].size());
        positions_adj_reverse[label].resize(objects[label].size());
    }

    // add label positions to adj
    positions_adj[0][0] = 0;
    positions_adj_reverse[0][0] = 0;

    for (uint32_t label = 0; label < L; label++){
        uint32_t count = positions_adj[label][0] + LabelCount[label];
        if (label < L-1) {
            positions_adj[label+1][0] = count;
            positions_adj_reverse[label+1][0] = count;
        }
    }
    std::vector<uint32_t >().swap(LabelCount);

    // add target positions to adj
    for (uint32_t label = 0; label < L; label++){
        uint32_t sourceIndex = positions_adj[label][0];
        uint32_t targetIndex = positions_adj[label][0];

        for (uint32_t i = 1; i < subjects[label].size(); i++) {
            auto prev = subjects[label][i-1];
            sourceIndex += LabelSource[label][prev];
            positions_adj[label][i] = sourceIndex;
        }
        for (uint32_t i = 1; i < objects[label].size(); i++) {
            auto prev = objects[label][i-1];
            targetIndex += LabelTarget[label][prev];
            positions_adj_reverse[label][i] = targetIndex;
        }
        std::vector<uint32_t >().swap(LabelSource[label]);
        std::vector<uint32_t >().swap(LabelTarget[label]);
    }
    std::cout << "Positions done\n";
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
            subjects[predicate].emplace_back(subject);
            objects[predicate].emplace_back(object);
        }
    }
    for (int i=0; i<L; i++) {
        sort(subjects[i].begin(), subjects[i].end());
        subjects[i].erase(unique(subjects[i].begin(), subjects[i].end()), subjects[i].end());
        sort(objects[i].begin(), objects[i].end());
        objects[i].erase(unique(objects[i].begin(), objects[i].end()), objects[i].end());
    }
    std::cout << "Initial info read\n";
    graphFile.close();
}

std::vector<uint32_t> SimpleGraph::findNeighbours(uint32_t id, uint32_t label, bool reverse) {
    std::vector<uint32_t> N = {};
    bool found = false;
    bool lastIndex = false;
    uint32_t index;
    if (!reverse) {
        for (index = 0; index < subjects[label].size(); index++){
            if (subjects[label][index] == id) {
                found = true;
                if (index == subjects[label].size()-1)
                    lastIndex = true;
                break;
            }
        }
        if (!found)
            return N;
        if (!lastIndex) {
            uint32_t first = positions_adj[label][index];
            uint32_t last = positions_adj[label][index + 1];
            while (first < last) {
                N.emplace_back(IA[first]);
                first++;
            }
        }
        else {
            uint32_t first = positions_adj[label][index];
            uint32_t last;
            if (label < L - 1)
                last = positions_adj[label + 1][0];
            else
                last = IA.size();
            while (first < last) {
                N.emplace_back(IA[first]);
                first++;
            }
        }
    }
    else {
        for (index = 0; index < objects[label].size(); index++){
            if (objects[label][index] == id){
                found = true;
                if (index == objects[label].size()-1)
                    lastIndex = true;
                break;
            }
        }
        if (!found)
            return N;
        if (!lastIndex) {
            uint32_t first = positions_adj_reverse[label][index];
            uint32_t last = positions_adj_reverse[label][index + 1];
            while (first < last) {
                N.emplace_back(IA_reverse[first]);
                first++;
            }
        }
        else {
            uint32_t first = positions_adj_reverse[label][index];
            uint32_t last;
            if (label < L - 1)
                last = positions_adj_reverse[label + 1][0];
            else
                last = IA_reverse.size();
            while (first < last) {
                N.emplace_back(IA_reverse[first]);
                first++;
            }
        }
    }
    return N;
}
