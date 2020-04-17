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

// void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel, std::vector<std::vector<uint32_t>> &offset, std::vector<std::vector<uint32_t>> &rev_offset) {
//     int index;
//     for (index = 0; index < subjects[edgeLabel].size(); index++){
//         if (subjects[edgeLabel][index] == from) {
//             uint32_t i = offset[edgeLabel][index];
//             offset[edgeLabel][index]++;
//             IA[positions_adj[edgeLabel][index] + i] = to;
//             break;
//         }
//     }
//     for (index = 0; index < objects[edgeLabel].size(); index++){
//         if (objects[edgeLabel][index] == to) {
//             uint32_t i_reverse = rev_offset[edgeLabel][index];
//             rev_offset[edgeLabel][index]++;
//             IA_reverse[positions_adj_reverse[edgeLabel][index] + i_reverse] = from;
//             break;
//         }
//     }
// }

// bool SimpleGraph::edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel) {
//     std::vector<uint32_t> N = findNeighbours(from, edgeLabel, false);
//     for (int i = 0; i < N.size(); i++) {
//         if (N[i] == to)
//             return true;
//     }
//     return false;
// }

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
            // if (i > positions_adj[predicate][subject+1])
            //     std::cout << "index too big:" << i << " source:" << subject << " label:" << predicate << " target:" << object << std::endl;
            // if (i > IA.size())
            //     std::cout << predicate << "    " << subject << "   " << i << std::endl;
            IA[i] = object;
            IA_reverse[i] = subject;
        }
    }

    // // test positio_adj & IA
    // std::cout << IA[0] << " " << IA[1] << " " << IA[2] << " " << IA[3] << " " << IA[4] << " " << IA[5] << std::endl;
    // std::cout << positions_adj[0][0] << " " << positions_adj[0][1] << " " << positions_adj[0][2] << " " << positions_adj[0][3] << " " << positions_adj[0][4] << std::endl;

    graphFile.close();

    std::pair<uint32_t, uint32_t> x = SelectLabel(1, false);
    std::cout << x.first << " " << x.second << std::endl;
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

    // std::cout << LabelSource[3][9725] << "\n";

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

    // for (uint32_t label = 0; label < LabelTarget.size(); label++){
    //     uint32_t sourceIndex = positions_adj[label][0];
    //     for (uint32_t source = 1; source < LabelTarget[label].size(); source++) {
    //         sourceIndex += LabelTarget[label][source-1];
    //         positions_adj[label][source] = sourceIndex;
    //     }
    // }

    // test
    // uint32_t max_index = 0;
    // uint32_t last_index = 0;
    // uint32_t index;
    // for (uint32_t label = 0; label < positions_adj.size(); label++) {
    //     for (uint32_t source = 0; source < positions_adj[label].size(); source++) {
    //         index = positions_adj[label][source];
    //         max_index = std::max(max_index,index);
    //         if (index < last_index)
    //             std::cout << "index too small" << std::endl;
    //     }
    // }
    // std::cout << "max " << max_index << " " << IA.size() <<std::endl;

    // std::cout << positions_adj[0][0] << " " << positions_adj[0][2] <<" " << positions_adj[0][3] <<" " << positions_adj[0][4] << std::endl;
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

// std::vector<uint32_t> SimpleGraph::findNeighbours(uint32_t id, uint32_t label, bool reverse) {
//     std::vector<uint32_t> N = {};
//     bool found = false;
//     bool lastIndex = false;
//     uint32_t index;
//     if (!reverse) {
//         for (index = 0; index < subjects[label].size(); index++){
//             if (subjects[label][index] == id) {
//                 found = true;
//                 if (index == subjects[label].size()-1)
//                     lastIndex = true;
//                 break;
//             }
//         }
//         if (!found)
//             return N;
//         if (!lastIndex) {
//             uint32_t first = positions_adj[label][index];
//             uint32_t last = positions_adj[label][index + 1];
//             while (first < last) {
//                 N.emplace_back(IA[first]);
//                 first++;
//             }
//         }
//         else {
//             uint32_t first = positions_adj[label][index];
//             uint32_t last;
//             if (label < L - 1)
//                 last = positions_adj[label + 1][0];
//             else
//                 last = IA.size();
//             while (first < last) {
//                 N.emplace_back(IA[first]);
//                 first++;
//             }
//         }
//     }
//     else {
//         for (index = 0; index < objects[label].size(); index++){
//             if (objects[label][index] == id){
//                 found = true;
//                 if (index == objects[label].size()-1)
//                     lastIndex = true;
//                 break;
//             }
//         }
//         if (!found)
//             return N;
//         if (!lastIndex) {
//             uint32_t first = positions_adj_reverse[label][index];
//             uint32_t last = positions_adj_reverse[label][index + 1];
//             while (first < last) {
//                 N.emplace_back(IA_reverse[first]);
//                 first++;
//             }
//         }
//         else {
//             uint32_t first = positions_adj_reverse[label][index];
//             uint32_t last;
//             if (label < L - 1)
//                 last = positions_adj_reverse[label + 1][0];
//             else
//                 last = IA_reverse.size();
//             while (first < last) {
//                 N.emplace_back(IA_reverse[first]);
//                 first++;
//             }
//         }
//     }
//     return N;
// }

// uint32_t SimpleGraph::getIn() {
//     return subjects.size();
// }

// uint32_t SimpleGraph::getOut() {
//     return objects.size();
// }

// uint32_t SimpleGraph::getPaths() {
//     return getNoEdges();
// }

// void SimpleGraph::sortNodes() {
//     for (int i=0; i<L; i++) {
//         if (subjects[i].size() > 1) {
//             sort(subjects[i].begin(), subjects[i].end());
//             subjects[i].erase(unique(subjects[i].begin(), subjects[i].end()), subjects[i].end());
//             sort(objects[i].begin(), objects[i].end());
//             objects[i].erase(unique(objects[i].begin(), objects[i].end()), objects[i].end());
//         }
//     }
// }

// uint32_t SimpleGraph::getLabelEdgeCount(uint32_t label, bool reverse) {
//     if (!reverse) {
//         if (label < L - 1)
//             return positions_adj[label+1][0] - positions_adj[label][0];
//         else
//             return getNoEdges() - positions_adj[label][0];
//     }
//     else {
//         if (label < L - 1)
//             return positions_adj_reverse[label+1][0] - positions_adj_reverse[label][0];
//         else
//             return getNoEdges() - positions_adj_reverse[label][0];
//     }
// }

// std::vector<uint32_t> SimpleGraph::getLabelSources(uint32_t label, bool reverse) {
//     if (!reverse)
//         return subjects[label];
//     else
//         return objects[label];
// }

// std::vector<uint32_t> SimpleGraph::getLabelTargets(uint32_t label, bool reverse) {
//     if (!reverse)
//         return objects[label];
//     else
//         return subjects[label];
// }

// void SimpleGraph::setLabelSources(uint32_t label, std::vector<uint32_t> sources) {
//     subjects[label] = sources;
// }

// void SimpleGraph::setLabelTargets(uint32_t label, std::vector<uint32_t> targets) {
//     objects[label] = targets;
// }

// void SimpleGraph::setLabelCount(uint32_t label, uint32_t count) {
//     LabelCount[label] = count;
// }

// void SimpleGraph::addLabelSource(uint32_t label, uint32_t source) {
//     LabelSource[label][source]++;
// }

// void SimpleGraph::addLabelTarget(uint32_t label, uint32_t target) {
//     LabelTarget[label][target]++;
// }
