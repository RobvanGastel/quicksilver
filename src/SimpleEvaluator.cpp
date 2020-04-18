#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include <sstream>

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

cardStat SimpleEvaluator::computeStats(std::vector<std::pair<uint32_t, uint32_t>> &g) {
    cardStat stats {};

    std::vector<uint32_t> sources = {};
    std::vector<uint32_t> targets = {};
    std::vector<std::pair<uint32_t, uint32_t>> paths = {};

    for(int i = 0; i < g.size(); i++) {
        sources.emplace_back(g[i].first);
        targets.emplace_back(g[i].second);

        bool exists = false;
        for(int j = 0; j < paths.size(); j++) {
            if(paths[j] == g[i]) {
                exists = true;
                break;
            }
        }
        if(!exists) {
            paths.emplace_back(g[i]);
        }
    }
    
    sort(sources.begin(), sources.end());
    sources.erase(unique(sources.begin(), sources.end()), sources.end());

    sort(targets.begin(), targets.end());
    targets.erase(unique(targets.begin(), targets.end()), targets.end());

    stats.noIn = sources.size();
    stats.noPaths = paths.size();
    stats.noOut = targets.size();

    return stats;
}

/**
 * Merges a graph into another graph.
 * @param left A graph to be merged into.
 * @param right A graph to be merged from.
 * @return A number of distinct new edges added from the "right" graph into the "left" graph.
 */
uint32_t SimpleEvaluator::unionDistinct(std::vector<std::pair<uint32_t, uint32_t>> &left, 
            std::vector<std::pair<uint32_t, uint32_t>> &right) {

    uint32_t numNewAdded = 0;
    return numNewAdded;
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

    return std::vector<std::pair<uint32_t, uint32_t>> {};
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

            if(s != -1 && t != -1) return graph->SelectSTL(s, t, label, true); // 42, 1>, 43
            if(s == -1 && t == -1) return graph->SelectLabel(label, true); // *, 1>, *
            if(s != -1) return graph->SelectIdLabel(s, label, false); // *, 1>, 42
            if(t != -1) return graph->SelectIdLabel(t, label, true); // 42, 1>, *
            // return SimpleEvaluator::selectLabel(label, label, false, graph);
        } else if (std::regex_search(q->data, matches, inverseLabel)) {
            // Case: 1<
            label = (uint32_t) std::stoul(matches[1]);

            if(s != -1 && t != -1) return graph->SelectSTL(s, t, label, false); // 42, 1>, 43
            if(s == -1 && t == -1) return graph->SelectLabel(label, false); // *, 1>, *
            if(s != -1) return graph->SelectIdLabel(s, label, false); // *, 1>, 42
            if(t != -1) return graph->SelectIdLabel(t, label, true); // 42, 1>, *
        }
        else if(std::regex_search(q->data, matches, kleeneStar)) {
            // Case: 1+

            // TODO: Implement the TC
            // label = (uint32_t) std::stoul(matches[1]);
            // return SimpleEvaluator::transitiveClosure(label, graph);
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
        }

        return std::vector<std::pair<uint32_t, uint32_t>> {};
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftGraph = SimpleEvaluator::evaluatePath(q->left, -1, -1);
        auto rightGraph = SimpleEvaluator::evaluatePath(q->right, -1, -1);

        // join left with right
        return SimpleEvaluator::join(leftGraph, rightGraph);

    }
    

    return std::vector<std::pair<uint32_t, uint32_t>> {};
}

///
/// FindBestPlan
///
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

    std::cout << "\n\nUnprocessed query:" << query;
    PathTree* tree = PathTree::strToTree(query);
    std::cout << "\nProcessing query: " << *tree;

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
        res = evaluatePath(query->path, s, t);
        
        // Findbestplan DP
        // TODO: Manually merge findBestPlan from bestPlan
    }

    return SimpleEvaluator::computeStats(res);
}