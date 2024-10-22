#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include <sstream>

bool sortbysec(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

SimpleEvaluator::SimpleEvaluator(std::shared_ptr<SimpleGraph> &g) {

    // works only with SimpleGraph
    graph = g;
    est = nullptr; // estimator not attached by default
}

void SimpleEvaluator::attachEstimator(std::shared_ptr<SimpleEstimator> &e) {
    est = e;
}

void SimpleEvaluator::prepare() {

    // if attached, prepare the estimator
    if(est != nullptr) est->prepare();

    // prepare other things here.., if necessary

}

/// Helper functions to parse
void inorderParseTree(PathTree *node,
                  std::vector<std::string> *query) {
    if (node == nullptr) {
        return;
    }
    inorderParseTree(node->left, query);

    if (node->data != "/") {
        query->push_back(node->data);
    }

    inorderParseTree(node->right, query);
}

// Parse the given PathTree
std::vector<std::string> parsePathToTree(PathTree *tree) {
    std::vector<std::string> query;

    if (!tree->isLeaf()) {
        inorderParseTree(tree, &query);
    } else {
        query.push_back(tree->data);
    }
    return query;
}

PathTree* reverseJoinOrder(PathTree* p) {
    std::vector<std::string> path = parsePathToTree(p);
    std::string query;

    if(path.size() >= 3) {
        std::reverse(path.begin(), path.end());

        query = "(" + path[1] + "/" + path[0] + ")";
        for(int i = 2; i < path.size(); i++) {
            query = "(" + path[2] + "/" + query + ")";
        }
    } else {
        return p;
    }

    return PathTree::strToTree(query);
}

cardStat SimpleEvaluator::computeStats(std::vector<std::pair<uint32_t, uint32_t>> &g) {
    cardStat stats {};

    std::unordered_set<uint32_t> sources = {};
    std::unordered_set<uint32_t> targets = {};
    std::vector<std::pair<uint32_t, uint32_t>> tuples(g.begin(), g.end());
    std::sort(tuples.begin(), tuples.end());

    for (auto tuple : g) {
        sources.insert(tuple.first);
        targets.insert(tuple.second);
    }

    stats.noIn = targets.size();
    stats.noPaths = std::unique(tuples.begin(), tuples.end()) - tuples.begin();
    stats.noOut = sources.size();
    
    return stats;
}

/**
 * Simple implementation of a join of two graphs.
 * @param left A graph to be joined.
 * @param right Another graph to join with.
 * @return Answer graph for a join. Note that all labels in the answer graph are "0".
 */
std::vector<std::pair<uint32_t, uint32_t>>  SimpleEvaluator::join(
            std::vector<std::pair<uint32_t, uint32_t>> &left, 
            std::vector<std::pair<uint32_t, uint32_t>>  &right) {
    std::vector<std::pair<uint32_t,uint32_t>> join;
    uint32_t joinTarget = 0;
    uint32_t source;
    uint32_t left_i = 0;
    uint32_t left_max = left.size();
    uint32_t leftJoinTarget = 0;
    uint32_t right_i = 0;
    uint32_t right_step = 0;
    uint32_t right_max = right.size();
    
    uint32_t join_max_id = 0;
    for (auto p : left) {
        if (p.second > join_max_id) {
            join_max_id = p.second;
        }
    }
    for (auto p : right) {
        if (p.first > join_max_id) {
            join_max_id = p.first;
        }
    }

    std::vector<std::vector<uint32_t>> left_adj;
    left_adj.resize(join_max_id+1);
    for (auto p : left) {
        left_adj[p.second].emplace_back(p.first);
    }
    std::vector<std::vector<uint32_t>> right_adj;
    right_adj.resize(join_max_id+1);
    for (auto p : right) {
        right_adj[p.first].emplace_back(p.second);
    }

    // remove duplicates
    for (uint32_t i = 0; i < left_adj.size(); i++) {
        sort(left_adj[i].begin(), left_adj[i].end());
        left_adj[i].erase( unique( left_adj[i].begin(), left_adj[i].end() ), left_adj[i].end() );
        sort(right_adj[i].begin(), right_adj[i].end());
        right_adj[i].erase( unique( right_adj[i].begin(), right_adj[i].end() ), right_adj[i].end() );
    }

    for (uint32_t join_id = 0; join_id <= join_max_id; join_id++) {
        if ((!left_adj[join_id].empty()) && (!right_adj[join_id].empty())) {
            for (uint32_t lp : left_adj[join_id]) {
                for (uint32_t rp : right_adj[join_id]) {
                    join.emplace_back(std::make_pair(lp, rp));
                }
            }
        }
    }
    return join;
}

/**
 * Given an AST, evaluate the query and produce an answer graph.
 * @param q Parsed AST.
 * @return Solution as a graph.
 */
std::vector<std::pair<uint32_t, uint32_t>> SimpleEvaluator::evaluatePath(PathTree *q, int s, int t) {
    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {
        // selectLabel out the label in the AST
        std::regex directLabel(R"((\d+)\>)");
        std::regex inverseLabel(R"((\d+)\<)");
        std::regex kleeneStar(R"((\d+)\+)");

        std::smatch matches;

        uint32_t label;
        bool inverse;
        
        if (std::regex_search(q->data, matches, directLabel)) {
            // Case: 1>
            label = (uint32_t) std::stoul(matches[1]);

            if(s != -1 && t != -1) return graph->SelectSTL(s, t, label, false); // 42, 1>, 43
            if(s == -1 && t == -1) return graph->SelectLabel(label, false); // *, 1>, *
            if(s != -1) return graph->SelectIdLabel(s, label, false, false); // 42, 1>, *
            if(t != -1) return graph->SelectIdLabel(t, label, false, true); // *, 1>, 42
            // return SimpleEvaluator::selectLabel(label, label, false, graph);
        } else if (std::regex_search(q->data, matches, inverseLabel)) {
            // Case: 1<
            label = (uint32_t) std::stoul(matches[1]);
            
            if(s != -1 && t != -1) return graph->SelectSTL(s, t, label, true); // 42, 1<, 43
            if(s == -1 && t == -1) return graph->SelectLabel(label, true); // *, 1<, *
            if(s != -1) return graph->SelectIdLabel(s, label, true, false); // 42, 1<, *
            if(t != -1) return graph->SelectIdLabel(t, label, true, true); // *, 1<, 42
        }

        else if(std::regex_search(q->data, matches, kleeneStar)) {
            // Case: 1+
            label = (uint32_t) std::stoul(matches[1]);

            // if(s != -1 && t != -1) return graph->TC(s, t, label, true); // 42, 1+, 43
            if(s == -1 && t == -1) return graph->transitiveClosure(label); // *, 1+, *
            if(s != -1) return graph->transitiveClosure(label, s); // 42, 1+, *
            // if(t != -1) return graph->TC(t, label); // *, 1+, 42
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
        }

        return std::vector<std::pair<uint32_t, uint32_t>> {};
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftPairs = SimpleEvaluator::evaluatePath(q->left, s, -1);
        auto rightPairs = SimpleEvaluator::evaluatePath(q->right, -1, t);

        // sort(leftPairs.begin(), leftPairs.end(), sortbysec);
        // sort(rightPairs.begin(), rightPairs.end());

        // join left with right
        return SimpleEvaluator::join(leftPairs, rightPairs);

    }
    
    return std::vector<std::pair<uint32_t, uint32_t>> {};
}


///
/// FindBestPlan
///
std::vector<std::string> parseStringToPath(const std::string &query) {
    std::stringstream ss(query);
    std::string label;
    std::vector<std::string> path;
    while (std::getline(ss, label, '/')) {
        path.push_back(label);
    }
    return path;
}

// Smaller bestPlan object to store string instead

uint32_t SimpleEvaluator::estimateQueryCost(std::string left, std::string right) {
    std::string query = "";
    if(right != "") {
        query = "(" + left + "/" + right + ")";
    } else {
        query = "(" + left + ")";
    }

    PathTree* tree = PathTree::strToTree(query);

    std::string asterisk = "*";
    PathQuery* pq = new PathQuery(asterisk, tree, asterisk);

    auto noPaths = est->estimate(pq).noPaths;
    delete pq;
    return noPaths;
}

BestPlan SimpleEvaluator::findBestPlan(std::string query) {
    std::vector<std::string> path = parseStringToPath(query);

    // If already cached
    if(cachedPlans.count(query) != 0) {
        planSpace[query].query = cachedPlans[query].query;
        planSpace[query].cost = cachedPlans[query].cost;
    } else {
        if(path.size() == 1) {
            if(planSpace.count(query) == 0) {
                std::string left = path[0];

                BestPlan ps;
                ps.query = query;
                ps.cost = estimateQueryCost(path[0], "");
                planSpace[query] = ps;
            }
        } else if(path.size() == 2) {
            if(planSpace.count(query) == 0) {
                std::string left = path[0];
                std::string right = path[1];

                BestPlan ps;
                ps.query = "(" + left + "/" + right + ")";
                ps.cost = estimateQueryCost(path[0], path[1]);
                planSpace[query] = ps;
            }
        } else {
            for(int i = 0; i < path.size()-1; ++i) {
                std::string left = "";
                std::string right = "";
                for (int j = 0; j < path.size(); ++j) {
                    if (j <= i) {
                        left += path[j] + "/";
                    } else {
                        right += path[j] + "/";
                    }
                }
                left = left.substr(0, left.size()-1);
                right = right.substr(0, right.size()-1);

                BestPlan P1 = SimpleEvaluator::findBestPlan(left);
                BestPlan P2 = SimpleEvaluator::findBestPlan(right);

                std::string leftLabel = P1.query;
                std::string rightLabel = P2.query;
                uint32_t cost = estimateQueryCost(leftLabel, rightLabel);
                
                if(planSpace.count(query) == 0) {
                    // TODO: Adjust?
                    BestPlan P;
                    P.query = ""; //query
                    P.cost = INT32_MAX;
                    planSpace[query] = P;
                }
                if(cost < planSpace[query].cost) {
                    std::string P1Query = "(" + P1.query + ")";
                    std::string P2Query = "(" + P2.query + ")";
                    planSpace[query].query = P1Query + "/" + P2Query;
                    planSpace[query].cost = cost;
                }
            }
        }
    }
    return planSpace[query];
}

std::string SimpleEvaluator::PathQueryBestPlan(PathQuery* query) {
    std::vector<std::string> path = parsePathToTree(query->path);
    std::string queryString = "" + path[0];
    for(int i = 1; i < path.size(); i++) {
        queryString += "/" + path[i];
    }
    auto plan = findBestPlan(queryString);
    planSpace.clear();
    
    query->path = PathTree::strToTree(plan.query);
    return queryString;
}

/**
 * Evaluate a path query. Produce a cardinality of the answer graph.
 * @param query Query to evaluate.
 * @return A cardinality statistics of the answer graph.
 */
cardStat SimpleEvaluator::evaluate(PathQuery *query) {
    std::vector<std::pair<uint32_t, uint32_t>> res = {};

    int s = -1;
    int t = -1;

    if(query->s != "*" && query->t != "*") {
        s = std::stoi(query->s);
        t = std::stoi(query->t);
        res = evaluatePath(query->path, s, t);
    } 
    // evaluate source id first and continue from left to right
    else if(query->s != "*") {
        s = std::stoi(query->s);
        res = evaluatePath(query->path, s, t);
    }
    // evaluate target id first and continue from right to left
    else if(query->t != "*") { 
        t = std::stoi(query->t);
        res = evaluatePath(query->path, s, t);
    } else {
        // Findbestplan DP
        // Mutates query->path in place
        std::string qString = PathQueryBestPlan(query);

        // TODO: Possibility to cache query results
        // Cache all queries for entire query execution duration
        // cachedQuery[qString] = std::make_shared<SimpleGraph>(*res);
        
        res = evaluatePath(query->path, s, t);
    }

    return SimpleEvaluator::computeStats(res);
}