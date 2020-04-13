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
    // std::cout << node->data << std::endl;
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
    // TODO keep joined tuples

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
    // std::cout << "Created trees: " << trees.size() << std::endl;
    return trees;
} 

void inorderWalkRemove(PathTree *node, std::string remove) {
    if (node == nullptr) {
        return;
    }
    if(node->left != nullptr) {
        if (node->left->data == remove) {
            node->left = NULL;
        } else {
            inorderWalkRemove(node->left, remove);
        }
    }

    if(node->right != nullptr) {
        if (node->right->data == remove) {
            node->right = NULL;
        } else {
            inorderWalkRemove(node->right, remove);
        }
    }
}

void setDifference(BestPlan *S, std::string S1) {
    inorderWalkRemove(S->plan->path, S1);
}

bool containsOneRelation(BestPlan S) {
    // return S.plan.left.isLeaf();
    // std::cout << "One or less joins: " << (parsePathToTree(S.plan->path).size() <= 2) << std::endl;
    return parsePathToTree(S.plan->path).size() <= 2;
}

PathQuery* constructNewTree(BestPlan S, BestPlan P1, std::pair<std::string, std::string> combination) {
    auto path = parsePathToTree(S.plan->path);

    // TODO: Make sure source target are correct
    // Currently not implemented
    std::string asterisk = "*";
    std::string join = "/";

    PathTree* left;
    PathTree* right;
    for(int i = 0; i < path.size(); i++) {

        // Right node
        if(i % 2 == 1) {
            if(path[i] == combination.first) {
                right = P1.plan->path;
                continue;
            }
            right = new PathTree(path[i], nullptr, nullptr);
        }

        // Left node
        if (i % 2 == 0) {
            if(path[i] == combination.first) {
                left = P1.plan->path;
                continue;
            }
            left = new PathTree(path[i], left, right);
        }
    }
    return new PathQuery(asterisk, left, asterisk);
}

void determineOrder(std::vector<std::string> paths) {

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
        std::cout << "\n";

        // for each non-empty subset S1 of S such that S1 != S
        for (std::pair<PathQuery*, std::pair<std::string, std::string>> q : queries) {
            auto S1 = BestPlan(INT_MAX, q.first);
            auto P1 = findBestPlan(S1);

            // TODO: Split up recursive look up for up down search in
            // determineOrder.

            // Loop over all possibilities except the combination join
            std::vector<std::string> path = parsePathToTree(S.plan->path);
            std::vector<std::string> joinsDown;
            std::vector<std::string> joinsUp;
            for(int i = 0; i < path.size(); i++) {
                // Find amount to go down
                if(path[i] == q.second.first && path[i+1] == q.second.second) {
                    // Check if there exist an element
                    if (i > 0) {
                        for(int j = i-1; j > -1; j--) {
                            joinsDown.push_back(path[j]);
                        }
                    }
                }
                // Find amount to go up
                if(path[i] == q.second.second && path[i-1] == q.second.first) {
                    if (i < path.size() - 1) {
                        for(int j = i+1; j < path.size(); j++) {   
                            joinsUp.push_back(path[j]);
                        }                    
                    }
                }
            }

            std::cout << "combination: " << q.second.first << ", " << q.second.second << std::endl;

            std::cout << "joinsDown: ";
            for(int i = 0; i < joinsDown.size(); i++) {
                std::cout << joinsDown[i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "joinsUp: ";
            for(int i = 0; i < joinsUp.size(); i++) {
                std::cout << joinsUp[i] << ", ";
            }
            std::cout << "\n" << std::endl;

            std::string join = "/";
            std::string asterisk = "*";

            // JoinDown tree
            BestPlan S2 = BestPlan(INT_MAX, S.plan);

            // Store root node
            PathTree* root = q.first->path;
            PathTree* tempNode = q.first->path;
            
            if(joinsDown.size() != 0) {
                for(int i = 0; i < joinsDown.size(); i++) {
                    
                    // Need to create a join of case
                    //    join
                    //   /    \
                    //  1   combi-join
                    //        /    \
                    //       2      3
                    if(joinsDown.size()-i == 1) {
                        if(i == 0) { // set root node pointer I want to keep
                            PathTree* node = new PathTree(joinsDown[i], nullptr, nullptr);
                            root = new PathTree(join, node, tempNode);
                        } else {
                            PathTree* node = new PathTree(joinsDown[i], nullptr, nullptr);
                            root = new PathTree(join, node, tempNode);
                        }
                    // Need to create a join of case
                    //          join
                    //         /    \
                    //     join    combi-join
                    //    /   \      /    \
                    //   1     4    2      3
                    } else if(joinsDown.size()-i == 2) {
                        PathTree* node1 = new PathTree(joinsDown[i], nullptr, nullptr);
                        i++;
                        PathTree* node4 = new PathTree(joinsDown[i], nullptr, nullptr);
                        PathTree* nodeJoin = new PathTree(join, node1, node4);

                        root = new PathTree(join, nodeJoin, tempNode);

                    // At least 3 joins
                    //          join
                    //         /    \
                    //     join    combi-join
                    //    /   \      /    \
                    //  join  ...   2      3
                    //  /  \
                    // 1    4
                    } else {

                    }
                }

                PathQuery* queryDown = new PathQuery(asterisk, root, asterisk);
                std::cout << "JoinsDown tree: " << "\n";
                std::cout << *queryDown->path << "\n";
                auto pathsdown = parsePathToTree(queryDown->path);
                for(int k = 0; k < pathsdown.size(); k++) {
                    std::cout << pathsdown[k] << ", ";
                }
                std::cout << std::endl;

                // Estimate cost
                int cost = est->estimate(queryDown).noPaths + P1.cost;
                std::cout << "cost down: " << cost << "\n";
                S2 = BestPlan(cost, queryDown);
            }

            // JoinUp tree
            // if(joinsUp.size() != 0) {
            //     // Combine up tree with down tree if the down tree exists
            //     PathTree* combined = q.first->path;
            //     if(S2.cost != INT_MAX) {
            //         combined = S2.plan->path;
            //     } 
                
            //     // Give temp value to start off with
            //     PathTree* joinUp = q.first->path;
            //     PathTree* unfinishedFirstNode = q.first->path;
            //     PathTree* secondNode = q.first->path;
            //     for(int i = 0; i < joinsUp.size(); i++) {
            //         PathTree* joinLeaf = new PathTree(joinsUp[i], nullptr, nullptr);

            //         if(i == 0) {
            //             // After the evaluation the nullptr will be filled with
            //             // joinDown.
            //             unfinishedFirstNode = new PathTree(join, nullptr, joinLeaf);
            //         } else {
            //             secondNode = new PathTree(join, unfinishedFirstNode, joinLeaf);
            //         }
            //     }

            //     PathQuery* queryUpUnfinished = new PathQuery(asterisk, secondNode, asterisk);
            //     std::cout << *secondNode << "\n";

            //     // Estimate cost
            //     int cost = 0;
            //     if(S2.cost != INT_MAX) {
            //         cost = est->estimate(queryUpUnfinished).noPaths + S2.cost + P1.cost;
            //     } else {
            //         cost = est->estimate(queryUpUnfinished).noPaths + P1.cost;
            //     }

            //     // Finish the actual tree after calculating seperate estimations
            //     unfinishedFirstNode->left = S2.plan->path;
            //     queryUpUnfinished = new PathQuery(asterisk, secondNode, asterisk);
            //     std::cout << *secondNode << "\n";

            //     // The tree is finished in above line!
            //     std::cout << "cost up: " << cost << "\n";
            //     S2 = BestPlan(cost, queryUpUnfinished);
            // }

            // // Determine which combination is better
            // if(S.cost > S2.cost) {
            //     std::cout << "new cost: " << S2.cost << ", previous cost: " << S.cost << std::endl;
            //     std::cout << "\n";
            //     S = S2;
            // }
        }
    }
    return S;
}

/**
 * Evaluate a path query. Produce a cardinality of the answer graph.
 * @param query Query to evaluate.
 * @return A cardinality statistics of the answer graph.
 */
cardStat SimpleEvaluator::evaluate(PathQuery *query) {

    // Start normal evaluation
    auto start = std::chrono::steady_clock::now();
    PathQuery* planQuery = query;
    auto res = evaluatePath(planQuery->path);
    if(query->s != "*") res = selectSource(planQuery->s, res);
    else if(query->t != "*") res = selectTarget(planQuery->t, res);
    auto end = std::chrono::steady_clock::now();
    std::cout << " " << std::endl;
    std::cout << "Time to evaluate: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

    // Start find best plan and evaluate
    std::cout << " " << std::endl;
    std::cout << "Start evaluating BestPlan" << std::endl;
    auto S = BestPlan(INT_MAX, query);
    auto SPlan = findBestPlan(S);
    std::cout << "end BestPlan, final cost: " << SPlan.cost << std::endl;

    start = std::chrono::steady_clock::now();
    planQuery = SPlan.plan;
    res = evaluatePath(planQuery->path);
    if(query->s != "*") res = selectSource(planQuery->s, res);
    else if(query->t != "*") res = selectTarget(planQuery->t, res);
    end = std::chrono::steady_clock::now();
    std::cout << "Time to evaluate BestPlan: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;
    std::cout << " " << std::endl;

    return SimpleEvaluator::computeStats(res);
}