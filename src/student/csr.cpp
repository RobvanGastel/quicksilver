#include "csr.h"

csr::csr(uint32_t n)   {
    setNoVertices(n);
}


std::pair<uint32_t, uint32_t> csr::SelectLabel(uint32_t label, bool reverse) {
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

uint32_t csr::getNoVertices() const {
    return V;
}

void csr::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t csr::getNoEdges() const {
    uint32_t sum = 0;
    // for (const auto & l : adj)
    //     sum += l.size();
    return sum;
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortPairs2(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

uint32_t csr::getNoDistinctEdges() const {

    uint32_t sum = 0;

    // for (auto sourceVec : adj) {

    //     std::sort(sourceVec.begin(), sourceVec.end(), sortPairs2);

    //     uint32_t prevTarget = 0;
    //     uint32_t prevLabel = 0;
    //     bool first = true;

    //     for (const auto &labelTgtPair : sourceVec) {
    //         if (first || !(prevTarget == labelTgtPair.second && prevLabel == labelTgtPair.first)) {
    //             first = false;
    //             sum++;
    //             prevTarget = labelTgtPair.second;
    //             prevLabel = labelTgtPair.first;
    //         }
    //     }
    // }

    return sum;
}

uint32_t csr::getNoLabels() const {
    return L;
}

void csr::setNoLabels(uint32_t noLabels) {
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

void csr::initialize_positions_adj() {
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

void csr::readFromContiguousFile(const std::string &fileName) {
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

void csr::readInitialInfoFromContiguousFile(const std::string &fileName) {
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
