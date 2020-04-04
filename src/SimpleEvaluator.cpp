#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
// #include <bits/stdc++.h>

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

void inorderPrint(PathTree *node) {
    if (node == nullptr) {
        return;
    }
    inorderPrint(node->left);
    std::cout << node->data << std::endl;
    inorderPrint(node->right);
}

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

std::vector<std::string> parsePathToTree(PathTree *tree) {
    std::vector<std::string> query;

    if (!tree->isLeaf()) {
        inorderParseTree(tree, &query);
    } else {
        query.push_back(tree->data);
    }
    return query;
}

std::vector<std::pair<PathQuery*, std::pair<std::string, std::string>>> createSubsetJoins(PathQuery *q) 
{
    // Parse path of tree to temporary vector<int>
    std::vector<std::string> path = parsePathToTree(q->path);
    std::cout << "Parsed tree to vector of int" << std::endl;

    // All possible pairs
	std::vector<std::pair<std::string, std::string>> combinations;
    for(int i = 1; i < path.size(); i++) {
        combinations.push_back(std::make_pair(path[i-1], path[i]));
    }

    std::vector<std::pair<PathQuery*, std::pair<std::string, std::string>>> trees;
    for(auto pair : combinations) {
        std::string first = pair.first;
        PathTree* left = new PathTree(first, nullptr, nullptr);

        std::string second = pair.second;
        PathTree* right = new PathTree(second, nullptr, nullptr);

        std::string join = "/";
        PathTree* tree = new PathTree(join, left, right);

        std::string asterisk = "*";
        // TODO: Make sure source target are correct
        // Currently not implemented
        PathQuery* pq = new PathQuery(asterisk, tree, asterisk);

        auto treepair = std::make_pair(pq, pair);
        trees.push_back(treepair);
    }
    std::cout << "Created trees: " << trees.size() << std::endl;
    return trees;
} 

void inorderWalkRemove(PathTree *node, std::string remove) {
    if (node == nullptr) {
        return;
    }

    if (node->left->data == remove) {
        std::cout << "remove left " << remove << std::endl;
        std::cout << node->left->data << std::endl;
        node->left = NULL;
    } else {
        std::cout << "dont remove left " << remove <<  std::endl;
        std::cout << node->left->data << std::endl;
        inorderWalkRemove(node->left, remove);
    }
    
    if (node->right->data == remove) {
        std::cout << "remove right " << remove << std::endl;
        std::cout << node->right->data << std::endl;
        node->right = NULL;
    } else {
        std::cout << "dont remove right " << remove << std::endl;
        std::cout << node->right->data << std::endl;
        inorderWalkRemove(node->right, remove);
    }    
}

void setDifference(BestPlan *S, std::string S1) {
    inorderWalkRemove(S->plan->path, S1);
}

bool containsOneRelation(BestPlan S) {
    // return S.plan.left.isLeaf();
    std::cout << "One join: " << (parsePathToTree(S.plan->path).size() <= 2) << std::endl;
    return parsePathToTree(S.plan->path).size() <= 2;
}

// TODO: Extend to join with TC
BestPlan SimpleEvaluator::findBestPlan(BestPlan S) {
    if (S.cost != INT_MAX) {
        return S;
    }
    if (containsOneRelation(S)) {
        // Based on the best way of assessing S
        S.cost = est->estimate(S.plan).noPaths;
        S.plan = S.plan;
    } else {
        // All possible join combinations
	    auto queries = createSubsetJoins(S.plan);
        std::cout << "Subsets created: " << queries.size() << std::endl;

        // for each non-empty subset S1 of S such that S1 != S
        for (std::pair<PathQuery*, std::pair<std::string, std::string>> q : queries) {
            auto S1 = BestPlan(INT_MAX, q.first);
            auto P1 = findBestPlan(S1);

            auto SMinusS1 = S.clone();
            // The corresponding left side of the join
            setDifference(SMinusS1, q.second.first);
            // std::cout << "dsadasd"<< std::endl;
            setDifference(SMinusS1, q.second.second);
            
            // std::cout << "q: " << q.second.first << ", " << q.second.second << std::endl;
            // inorderPrint(SMinusS1->plan->path);
            std::cout << "dadasd" << std::endl;
            inorderPrint(S.plan->path);

            auto P2 = findBestPlan(*SMinusS1);

            // TODO: Estimate the join of P1 and P2 "Naive"?
            // $A=$ best algorithm for joining results of $P1$ and $P2$
            int cost = P1.cost + P2.cost; // + cost of A
            if (cost < S.cost) {
                // TODO: Set the new cost if its lower and create the tree
                S.cost = cost;
                std::cout << "cost: " << cost << std::endl;
                // S.plan = execute P1. plan; execute P2. plan; join results of $P1$ and $P2$ using A
            }
        }
        std::cout << "end for loop combinations" << std::endl;
    }
    return S;
}

/**
 * Evaluate a path query. Produce a cardinality of the answer graph.
 * @param query Query to evaluate.
 * @return A cardinality statistics of the answer graph.
 */
cardStat SimpleEvaluator::evaluate(PathQuery *query) {
    /// Find best plan
    std::cout << std::endl;
    auto S = BestPlan(INT_MAX, query);
    std::cout << "BestPlan S created" << std::endl;
    findBestPlan(S);

    auto res = evaluatePath(query->path);
    if(query->s != "*") res = selectSource(query->s, res);
    else if(query->t != "*") res = selectTarget(query->t, res);
    return SimpleEvaluator::computeStats(res);
}