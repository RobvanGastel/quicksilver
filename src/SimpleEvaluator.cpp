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

cardStat SimpleEvaluator::computeStats(std::shared_ptr<SimpleGraph> &g) {

    cardStat stats {};

    for(int source = 0; source < g->getNoVertices(); source++) {
        if(!g->adj[source].empty()) stats.noOut++;
    }

    stats.noPaths = g->getNoDistinctEdges();

    for(int target = 0; target < g->getNoVertices(); target++) {
        if(!g->reverse_adj[target].empty()) stats.noIn++;
    }

    return stats;
}

/**
 * Select all edges from a graph with a given edge label.
 * @param projectLabel Label to select.
 * @param outLabel Label to rename the selected edge labels to (used in the TC).
 * @param inverse Follow the edges in inverse direction.
 * @param in The graph to select from.
 * @return The graph which only has edges with the specified edge label.
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::selectLabel(uint32_t projectLabel, uint32_t outLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if(!inverse) {
        // going forward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, outLabel);
            }
        }
    } else {
        // going backward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->reverse_adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, outLabel);
            }
        }
    }

    return out;
}

/**
 * A naive transitive closure (TC) computation.
 * @param label A graph edge label to compute the TC for.
 * @param in Input graph
 * @return
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::transitiveClosure(uint32_t label, std::shared_ptr<SimpleGraph> &in) {

    auto tc = selectLabel(label, 0, 0, in);
    auto base = selectLabel(label, 0, 0, in);
    uint32_t numNewAdded = 1;

    while (numNewAdded) {
        auto delta = join(tc, base);
        numNewAdded = unionDistinct(tc, delta);
    }

    return tc;
}

/**
 * Merges a graph into another graph.
 * @param left A graph to be merged into.
 * @param right A graph to be merged from.
 * @return A number of distinct new edges added from the "right" graph into the "left" graph.
 */
uint32_t SimpleEvaluator::unionDistinct(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    uint32_t numNewAdded = 0;

    for(uint32_t source = 0; source < right->getNoVertices(); source++) {
        for (auto labelTarget : right->adj[source]) {

            auto label = labelTarget.first;
            auto target = labelTarget.second;

            if(!left->edgeExists(source, target, label)) {
                left->addEdge(source, target, label);
                numNewAdded++;
            }
        }
    }

    return numNewAdded;
}

/**
 * Simple implementation of a join of two graphs.
 * @param left A graph to be joined.
 * @param right Another graph to join with.
 * @return Answer graph for a join. Note that all labels in the answer graph are "0".
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    auto out = std::make_shared<SimpleGraph>(left->getNoVertices());
    out->setNoLabels(1);

    for(uint32_t leftSource = 0; leftSource < left->getNoVertices(); leftSource++) {
        for (auto labelTarget : left->adj[leftSource]) {

            int leftTarget = labelTarget.second;
            // try to join the left target with right s
            for (auto rightLabelTarget : right->adj[leftTarget]) {

                auto rightTarget = rightLabelTarget.second;
                out->addEdge(leftSource, rightTarget, 0);

            }
        }
    }

    return out;
}

/**
 * Given an AST, evaluate the query and produce an answer graph.
 * @param q Parsed AST.
 * @return Solution as a graph.
 */
std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluatePath(PathTree *q) {

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
            label = (uint32_t) std::stoul(matches[1]);
            return SimpleEvaluator::selectLabel(label, label, false, graph);
        } else if (std::regex_search(q->data, matches, inverseLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            return SimpleEvaluator::selectLabel(label, label, true, graph);
        }
        else if(std::regex_search(q->data, matches, kleeneStar)) {
            label = (uint32_t) std::stoul(matches[1]);
            return SimpleEvaluator::transitiveClosure(label, graph);
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
        }

        return nullptr;
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftGraph = SimpleEvaluator::evaluatePath(q->left);
        auto rightGraph = SimpleEvaluator::evaluatePath(q->right);

        // join left with right
        return SimpleEvaluator::join(leftGraph, rightGraph);

    }

    return nullptr;
}

/**
 * Perform a selection on a source constant.
 * @param s A source constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<SimpleGraph> selectSource(std::string &s, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    for (auto labelTarget : in->adj[std::stoi(s)]) {

        auto label = labelTarget.first;
        auto target = labelTarget.second;

        out->addEdge((uint32_t) std::stoi(s), target, label);
    }

    return out;
}

/**
 * Perform a selection on a target constant.
 * @param s A target constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<SimpleGraph> selectTarget(std::string &t, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    for (auto labelSource : in->reverse_adj[std::stoi(t)]) {

        auto label = labelSource.first;
        auto source = labelSource.second;
        out->addEdge(source, (uint32_t) std::stoi(t), label);
    }

    return out;
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

    // Findbestplan DP
    std::vector<std::string> path = parsePathToTree(query->path);
    std::string queryString = "" + path[0];
    for(int i = 1; i < path.size(); i++) {
        queryString += "/" + path[i];
    }
    auto plan = findBestPlan(queryString);
    planSpace.clear();

    // Recreate new PathTree
    PathTree* tree = PathTree::strToTree(plan.query);

    auto res = evaluatePath(tree);
    if(query->s != "*") res = selectSource(query->s, res);
    else if(query->t != "*") res = selectTarget(query->t, res);
    return SimpleEvaluator::computeStats(res);
}