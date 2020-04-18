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


std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::SelectLabel(uint32_t label, bool reverse) {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (reverse) {
        for (uint32_t i = 0; i < V; i++) {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj_reverse[label][i],
                positions_adj_reverse[label][i+1]
            );
            for (uint32_t j = indexes.first; j < indexes.second; j++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(i, IA_reverse[j]));
        }
    } else {
        for (uint32_t i = 0; i < V; i++) {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj[label][i],
                positions_adj[label][i+1]
            );
            for (uint32_t j = indexes.first; j < indexes.second; j++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(i, IA[j]));
        }
    }

    return pairs;
}

std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::SelectIdLabel(uint32_t id, uint32_t label, bool reverse) {
    std::string asd = "normal";
    if (reverse)
        asd = "reverse";
    std::cout << "\n id:" << id << "   label: " << label << "   reverse: " << asd << "\n";
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    std::cout << id << " " << label << " " << positions_adj[label][id] <<std::endl;
    if (reverse) {
        auto indexes = std::pair<uint32_t, uint32_t>(
            positions_adj_reverse[label][id],
            positions_adj_reverse[label][id+1]
        );
        for (uint32_t i = indexes.first; i < indexes.second; i++)
            pairs.emplace_back(std::pair<uint32_t, uint32_t>(IA_reverse[i], id));
    } else {
        auto indexes = std::pair<uint32_t, uint32_t>(
            positions_adj[label][id],
            positions_adj[label][id+1]
        );
        for (uint32_t i = indexes.first; i < indexes.second; i++)
            pairs.emplace_back(std::pair<uint32_t, uint32_t>(id, IA[i]));
    }
    return pairs;
}

std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::SelectSTL(uint32_t source, uint32_t target, uint32_t label, bool reverse) {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (reverse) {
        auto indexes = std::pair<uint32_t, uint32_t>(
            positions_adj_reverse[label][source],
            positions_adj_reverse[label][source+1]
        );
        for (uint32_t i = indexes.first; i < indexes.second; i++) {
            if (IA_reverse[i] == target)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(target, source));
        }
    } else {
        auto indexes = std::pair<uint32_t, uint32_t>(
            positions_adj[label][source],
            positions_adj[label][source+1]
        );
        for (uint32_t i = indexes.first; i < indexes.second; i++){
            if (IA[i] == target)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(source, target));
        }
    }
    return pairs;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    // LabelCount.resize(L);
    // LabelSource.resize(L);
    // LabelTarget.resize(L);
    // for(int i = 0; i < L; i++) {
    //     LabelCount[i] = 0;
    //     std::vector<uint32_t> zeroes(V, 0);
    //     LabelSource[i] = zeroes;
    //     LabelTarget[i] = zeroes;
    // }
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {    std::string line;
    using_csr = true;

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)");  // subject predicate object .
    std::regex headerPat (R"((\d+),(\d+),(\d+))");      // noNodes,noEdges,noLabels

    std::vector<uint32_t> LabelCount;               // number of edges for each label
    std::vector<std::vector<uint32_t>> LabelSource; // number of edges for each label of each source
    std::vector<std::vector<uint32_t>> LabelTarget; // number of edges for each label of each target


    // READ INITIAL INFO
    // readInitialInfo(fileName);
    std::ifstream graphFile { fileName };

    // parse the header (1st line)
    std::getline(graphFile, line);
    std::smatch matches;
    if(std::regex_search(line, matches, headerPat)) {
        uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
        uint32_t noEdges = (uint32_t) std::stoul(matches[2]);
        uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

        setNoVertices(noNodes);
        setNoLabels(noLabels);

        LabelCount.resize(L);
        LabelSource.resize(L);
        LabelTarget.resize(L);
        for(int i = 0; i < L; i++) {
            LabelCount[i] = 0;
            std::vector<uint32_t> zeroes(V, 0);
            LabelSource[i] = zeroes;
            LabelTarget[i] = zeroes;
        }

        IA.resize(noEdges);
        IA_reverse.resize(noEdges);
    } 
    else {
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


    // INITIALIZE POSITIONS
    // initialize_positions_adj();
    positions_adj.resize(L);
    positions_adj_reverse.resize(L);

    for (uint32_t label = 0; label < L; label++) {
        positions_adj[label].resize(V+1);
        positions_adj_reverse[label].resize(V+1);
    }

    // add label posisitons to adj
    // positions_adj[0][0] = 0;
    // positions_adj_reverse[0][0] = 0;

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
        
        for (uint32_t i = 1; i < V; i++) {
            sourceIndex += LabelSource[label][i-1];
            targetIndex += LabelTarget[label][i-1];
            positions_adj[label][i] = sourceIndex;
            positions_adj_reverse[label][i] = targetIndex;
        }
    }

    // CREATE CSR
    // line;
    std::ifstream graphFile2 { fileName };

    std::getline(graphFile2, line);

    // create positions_adj
    std::vector<std::vector<uint32_t>> adj = positions_adj;
    std::vector<std::vector<uint32_t>> adj_reverse = positions_adj_reverse;
    while(std::getline(graphFile2, line)) {
        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            uint32_t i = adj[predicate][subject]++;
            uint32_t i_reverse = adj_reverse[predicate][object]++;

            IA[i] = object;
            IA_reverse[i_reverse] = subject;
        }
    }

    graphFile2.close();

    // Testing CSR selectIdLabel() & selectLabel()
    // std::pair<uint32_t, uint32_t> res = SelectLabel(1, true);
    // std::cout << "selectIdLabel first: " <<  res.first << "  second: " << res.second << std::endl;
    // res = SelectIdLabel(42, 1, false);
    // std::cout << "selectIdLabel first: " <<  res.first << "  second: " << res.second << std::endl;
}

// TODO: Make more generic
std::shared_ptr<SimpleGraph> SimpleGraph::createGraphSelectLabelSource(uint32_t source, uint32_t label, bool reverse) {
    
    // auto res = SelectIdLabel(source, label, reverse);
    // auto g = std::make_shared<SimpleGraph>();
    
    // // Use addEdge from Nikolay 
    // // adj 
    // // [source] -> [target]

    // uint32_t first = res.first; 
    // if (reverse) {
    //     for (uint32_t i = res.first; i < res.second; i++) {
    //         IA_reverse[i];
    //     }
    // } else {
    //     for (uint32_t i = res.first; i < res.second; i++) {
    //         IA[i];
    //     }
    // }

    // return g;
}

// TODO: Isn't used right? Delete if agree
/// add edge on the adjacency structure to evaluation of the query
/// can only be used if, graph is created for use of adj
// void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
//     if(from >= V || to >= V || edgeLabel >= L)
//         throw std::runtime_error(std::string("Edge data out of bounds: ") +
//                                          "(" + std::to_string(from) + "," + std::to_string(to) + "," +
//                                          std::to_string(edgeLabel) + ")");
//     adj[from].emplace_back(std::make_pair(edgeLabel, to));
// }
   // READ INITIAL INFO
    // readInitialInfo(fileName);
    // std::ifstream graphFile { fileName };

    // // parse the header (1st line)
    // std::getline(graphFile, line);
    // std::smatch matches;
    // if(std::regex_search(line, matches, headerPat)) {
    //     uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
    //     uint32_t noEdges = (uint32_t) std::stoul(matches[2]);
    //     uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

    //     setNoVertices(noNodes);
    //     setNoLabels(noLabels);

    //     LabelCount.resize(L);
    //     LabelSource.resize(L);
    //     LabelTarget.resize(L);
    //     for(int i = 0; i < L; i++) {
    //         LabelCount[i] = 0;
    //         std::vector<uint32_t> zeroes(V, 0);
    //         LabelSource[i] = zeroes;
    //         LabelTarget[i] = zeroes;
    //     }

    //     IA.resize(noEdges);
    //     IA_reverse.resize(noEdges);
    // } 
    // else {
    //     throw std::runtime_error(std::string("Invalid graph header!"));
    // }

    // // parse edge data
    // while(std::getline(graphFile, line)) {
    //     if(std::regex_search(line, matches, edgePat)) {
    //         uint32_t subject = (uint32_t) std::stoul(matches[1]);
    //         uint32_t predicate = (uint32_t) std::stoul(matches[2]);
    //         uint32_t object = (uint32_t) std::stoul(matches[3]);

    //         LabelCount[predicate]++;
    //         LabelSource[predicate][subject]++;
    //         LabelTarget[predicate][object]++;
    //     }
    // }

    // graphFile.close();


    // INITIALIZE POSITIONS
    // initialize_positions_adj();
    // positions_adj.resize(L);
    // positions_adj_reverse.resize(L);

    // for (uint32_t label = 0; label < L; label++) {
    //     positions_adj[label].resize(V+1);
    //     positions_adj_reverse[label].resize(V+1);
    // }

    // add label posisitons to adj
    // positions_adj[0][0] = 0;
    // positions_adj_reverse[0][0] = 0;

    // for (uint32_t label = 0; label < L; label++){
    //     uint32_t count = positions_adj[label][0] + LabelCount[label];

    //     if (label < L-1) {
    //         positions_adj[label+1][0] = count;
    //         positions_adj_reverse[label+1][0] = count;
    //     }

    //     positions_adj[label][V] = count;
    //     positions_adj_reverse[label][V] = count;
    // }

    // // add target positions to adj
    // for (uint32_t label = 0; label < L; label++){
    //     uint32_t sourceIndex = positions_adj[label][0];
    //     uint32_t targetIndex = positions_adj[label][0];
        
    //     for (uint32_t i = 1; i < V; i++) {
    //         sourceIndex += LabelSource[label][i-1];
    //         targetIndex += LabelTarget[label][i-1];
    //         positions_adj[label][i] = sourceIndex;
    //         positions_adj_reverse[label][i] = targetIndex;
    //     }
    // }

    // CREATE CSR
    // line;
    // std::ifstream graphFile2 { fileName };

    // std::getline(graphFile2, line);

    // create positions_adj
    // std::vector<std::vector<uint32_t>> adj = positions_adj;
    // std::vector<std::vector<uint32_t>> adj_reverse = positions_adj_reverse;
    // while(std::getline(graphFile2, line)) {
    //     if(std::regex_search(line, matches, edgePat)) {
    //         uint32_t subject = (uint32_t) std::stoul(matches[1]);
    //         uint32_t predicate = (uint32_t) std::stoul(matches[2]);
    //         uint32_t object = (uint32_t) std::stoul(matches[3]);

    //         uint32_t i = adj[predicate][subject]++;
    //         uint32_t i_reverse = adj_reverse[predicate][object]++;

    //         IA[i] = object;
    //         IA_reverse[i_reverse] = subject;
    //     }
    // }

    // graphFile2.close();

    // Testing CSR selectIdLabel() & selectLabel()
    // std::pair<uint32_t, uint32_t> res = SelectLabel(1, true);
    // std::cout << "selectIdLabel first: " <<  res.first << "  second: " << res.second << std::endl;
    // res = SelectIdLabel(42, 1, false);
    // std::cout << "selectIdLabel first: " <<  res.first << "  second: " << res.second << std::endl;
// }

// TODO: Make more generic
// std::shared_ptr<SimpleGraph> SimpleGraph::createGraphSelectLabelSource(uint32_t source, uint32_t label, bool reverse) {
    
    // auto res = SelectIdLabel(source, label, reverse);
    // auto g = std::make_shared<SimpleGraph>();
    
    // // Use addEdge from Nikolay 
    // // adj 
    // // [source] -> [target]

    // uint32_t first = res.first; 
    // if (reverse) {
    //     for (uint32_t i = res.first; i < res.second; i++) {
    //         IA_reverse[i];
    //     }
    // } else {
    //     for (uint32_t i = res.first; i < res.second; i++) {
    //         IA[i];
    //     }
    // }

    // return g;
// }

// TODO: Isn't used right? Delete if agree
/// add edge on the adjacency structure to evaluation of the query
/// can only be used if, graph is created for use of adj
// void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
//     if(from >= V || to >= V || edgeLabel >= L)
//         throw std::runtime_error(std::string("Edge data out of bounds: ") +
//                                          "(" + std::to_string(from) + "," + std::to_string(to) + "," +
//                                          std::to_string(edgeLabel) + ")");
//     adj[from].emplace_back(std::make_pair(edgeLabel, to));
// }
