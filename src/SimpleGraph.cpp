#include "SimpleGraph.h"

// struct pair_hash {
//     inline std::size_t operator()(const std::pair<uint32_t, uint32_t> & v) const {
//         return v.first*31+v.second;
//     }
// }

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


std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::SelectLabel(uint32_t label, bool reverse, bool sort_sec) {
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

    // if (sort_sec)
    //     sort(pairs.begin(), pairs.end(), sortbysec);
    // else
    //     sort(pairs.begin(), pairs.end());
    // pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::SelectIdLabel(uint32_t id, uint32_t label, bool reverse, bool isTarget) {

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (reverse) { // 1<
        // 42,1<,*
        if(!isTarget) {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj_reverse[label][id],
                positions_adj_reverse[label][id+1]
            );
            for (uint32_t i = indexes.first; i < indexes.second; i++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(id, IA_reverse[i]));  // not sure
        }
        // *,1<,42
        // 1<
        // *,1<,t => t,l,*
        else {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj[label][id],
                positions_adj[label][id+1]
            );
            for (uint32_t i = indexes.first; i < indexes.second; i++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(IA[i], id));  // not sure
        }
    } else { // 1>
        // 42,1>,*
        if(!isTarget) {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj[label][id],
                positions_adj[label][id+1]
            );
            for (uint32_t i = indexes.first; i < indexes.second; i++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(id, IA[i])); // correct way
        }
        // *,1>,42 
        else {
            auto indexes = std::pair<uint32_t, uint32_t>(
                positions_adj_reverse[label][id],
                positions_adj_reverse[label][id+1]
            );
            for (uint32_t i = indexes.first; i < indexes.second; i++)
                pairs.emplace_back(std::pair<uint32_t, uint32_t>(IA_reverse[i], id)); // not sure
        }
    }
    // sort(pairs.begin(), pairs.end());
    // pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());
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
    sort(pairs.begin(), pairs.end());
    // pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

// std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::TC(uint32_t label) {
//     std::vector<std::pair<uint32_t, uint32_t>> pairs;
    
//     uint32_t numNewAdded = 1;

//     while (numNewAdded) {
//         std::vector<std::pair<uint32_t, uint32_t>> delta = join(tc, base);
//         numNewAdded = unionDistinct(tc, delta);
//     }

//     // return tc;



//     return pairs;
// }

/**
 * A naive transitive closure (TC) computation.
 * @param label A graph edge label to compute the TC for.
 * @param in Input graph
 * @return
 */
std::vector<std::pair<uint32_t, uint32_t>> SimpleGraph::transitiveClosure(uint32_t label, int s, int t) {
    // std::unordered_set<std::pair<uint32_t, uint32_t>, pair_hash> tc_set;
    // std::vector<std::pair<uint32_t, uint32_t>> tc = SelectLabel(label, false);
    std::vector<std::pair<uint32_t, uint32_t>> base = SelectLabel(label, false);
    std::vector<std::pair<uint32_t, uint32_t>> tc;

    uint32_t numNewAdded = 1;
    uint32_t old_tc_len;
    uint32_t joinTarget = 0;
    // new
    uint32_t join_max_id = 0; //std::min(left[left.size()-1].second, right[right.size()-1].first);
    for (auto p : base) {
        if (p.second > join_max_id) {
            join_max_id = p.second;
        }
        if (p.first > join_max_id) {
            join_max_id = p.first;
        }
    }
    
    std::vector<std::vector<uint32_t>> left_adj;
    std::vector<std::vector<uint32_t>> right_adj;
    left_adj.resize(join_max_id+1);
    right_adj.resize(join_max_id+1);
    if (s == -1 && t == -1) {
        tc = base;
        for (auto p : base) {
            left_adj[p.second].emplace_back(p.first);
            right_adj[p.first].emplace_back(p.second);
        }
    } else if (s != -1) {
        for (auto p : base) {
            if (p.first == s){
                left_adj[p.second].emplace_back(s);
                tc.emplace_back(std::make_pair(s, p.second));
            }
            right_adj[p.first].emplace_back(p.second);
        }
    }
    
    for (uint32_t join_id = 0; join_id <= join_max_id; join_id++) {
        if ((!left_adj[join_id].empty()) && (!right_adj[join_id].empty())) {
            for (uint32_t lp : left_adj[join_id]) {
                for (uint32_t rp : right_adj[join_id]) {
                    tc.emplace_back(std::make_pair(lp, rp));
                }
            }
        }
    }
    std::sort(tc.begin(),tc.end());
    tc.erase(unique(tc.begin(), tc.end()), tc.end());

    while (numNewAdded) {
        old_tc_len = tc.size();
        right_adj = {};
        right_adj.resize(join_max_id+1);
        for (auto p : tc) {
            right_adj[p.first].emplace_back(p.second);
        }
        for (uint32_t join_id = 0; join_id <= join_max_id; join_id++) {
            if ((!left_adj[join_id].empty()) && (!right_adj[join_id].empty())) {
                for (uint32_t lp : left_adj[join_id]) {
                    for (uint32_t rp : right_adj[join_id]) {
                        tc.emplace_back(std::make_pair(lp, rp));
                    }
                }
            }
        }
        std::sort(tc.begin(),tc.end());
        tc.erase(unique(tc.begin(), tc.end()), tc.end());



    //     auto delta =  join(tc, base);
        numNewAdded = old_tc_len - tc.size();
    }

    return tc;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {    std::string line;
    using_csr = true;

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)");  // subject predicate object .
    std::regex headerPat (R"((\d+),(\d+),(\d+))");      // noNodes,noEdges,noLabels

    std::vector<uint32_t> LabelCount;               // number of edges for each label
    std::vector<std::vector<uint32_t>> LabelSource; // number of edges for each label of each source
    std::vector<std::vector<uint32_t>> LabelTarget; // number of edges for each label of each target

    //_____________________________________________________________________________________________________________________

    std::vector<std::vector<std::vector<uint32_t>>> temp_adj;
    uint32_t noNodes;
    uint32_t noEdges;
    uint32_t noLabels;

    std::ifstream graphFile { fileName };

    // parse the header (1st line)
    std::getline(graphFile, line);
    std::smatch matches;
    if(std::regex_search(line, matches, headerPat)) {
        noNodes = (uint32_t) std::stoul(matches[1]);
        noEdges = (uint32_t) std::stoul(matches[2]);
        noLabels = (uint32_t) std::stoul(matches[3]);
        setNoVertices(noNodes);
        setNoLabels(noLabels);
        LabelCount.resize(L);
        LabelSource.resize(L);
        LabelTarget.resize(L);
        temp_adj.resize(L);
        for(int i = 0; i < L; i++) {
            LabelCount[i] = 0;
            std::vector<uint32_t> zeroes(V, 0);
            LabelSource[i] = zeroes;
            LabelTarget[i] = zeroes;
            temp_adj[i].resize(noNodes);
        }
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

            bool exists = false;
            for (auto node : temp_adj[predicate][subject]) {
                if (node == object)
                    exists = true;
            }
            if (exists) {
                noEdges--;
            }
            else {
                temp_adj[predicate][subject].emplace_back(object);
                LabelCount[predicate]++;
                LabelSource[predicate][subject]++;
                LabelTarget[predicate][object]++;
            }
        }
    }
    IA.resize(noEdges);
    IA_reverse.resize(noEdges);

    graphFile.close();

    positions_adj.resize(L);
    positions_adj_reverse.resize(L);

    for (uint32_t label = 0; label < L; label++) {
        positions_adj[label].resize(V+1);
        positions_adj_reverse[label].resize(V+1);
    }

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

    // create positions_adj
    std::vector<std::vector<uint32_t>> adj = positions_adj;
    std::vector<std::vector<uint32_t>> adj_reverse = positions_adj_reverse;
     for (uint32_t predicate = 0; predicate < temp_adj.size(); predicate++){
        for (uint32_t subject = 0; subject < temp_adj[predicate].size(); subject++) {
            for (auto object : temp_adj[predicate][subject]) {
                uint32_t i = adj[predicate][subject]++;
                uint32_t i_reverse = adj_reverse[predicate][object]++;

                IA[i] = object;
                IA_reverse[i_reverse] = subject;
            }
        std::vector<uint32_t >().swap(temp_adj[predicate][subject]);
        }
    }
}