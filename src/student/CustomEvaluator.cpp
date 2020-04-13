#include "SimpleEstimator.h"
#include "CustomEvaluator.h"

CustomEvaluator::CustomEvaluator(std::shared_ptr<csr> &g) {

    // works only with csr
    graph = g;
    est = nullptr; // estimator not attached by default
}

void CustomEvaluator::attachEstimator(std::shared_ptr<SimpleEstimator> &e) {
    est = e;
}

void CustomEvaluator::prepare() {

    // if attached, prepare the estimator
    if(est != nullptr) est->prepare();

    // prepare other things here.., if necessary

}

cardStat CustomEvaluator::computeStats(std::shared_ptr<csr> &g) {

    cardStat stats {};

    // for(int source = 0; source < g->getNoVertices(); source++) {
    //     if(!g->adj[source].empty()) stats.noOut++;
    // }

    // stats.noPaths = g->getNoDistinctEdges();

    // for(int target = 0; target < g->getNoVertices(); target++) {
    //     if(!g->reverse_adj[target].empty()) stats.noIn++;
    // }

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
std::shared_ptr<csr> CustomEvaluator::selectLabel(uint32_t projectLabel, uint32_t outLabel, bool inverse, std::shared_ptr<csr> &in) {

    auto out = std::make_shared<csr>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if(!inverse) {
        // going forward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            // for (auto labelTarget : in->adj[source]) {

            //     auto label = labelTarget.first;
            //     auto target = labelTarget.second;

            //     if (label == projectLabel)
            //         out->addEdge(source, target, outLabel);
            // }
            int i=0;
        }
    } else {
        // going backward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            // for (auto labelTarget : in->reverse_adj[source]) {

            //     auto label = labelTarget.first;
            //     auto target = labelTarget.second;

            //     if (label == projectLabel)
            //         out->addEdge(source, target, outLabel);
            // }
            int i=0;
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
std::shared_ptr<csr> CustomEvaluator::transitiveClosure(uint32_t label, std::shared_ptr<csr> &in) {

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
uint32_t CustomEvaluator::unionDistinct(std::shared_ptr<csr> &left, std::shared_ptr<csr> &right) {

    uint32_t numNewAdded = 0;

    for(uint32_t source = 0; source < right->getNoVertices(); source++) {
        // for (auto labelTarget : right->adj[source]) {

        //     auto label = labelTarget.first;
        //     auto target = labelTarget.second;

        //     if(!left->edgeExists(source, target, label)) {
        //         left->addEdge(source, target, label);
        //         numNewAdded++;
        //     }
        // }
        int i=0;
    }

    return numNewAdded;
}

/**
 * Simple implementation of a join of two graphs.
 * @param left A graph to be joined.
 * @param right Another graph to join with.
 * @return Answer graph for a join. Note that all labels in the answer graph are "0".
 */
std::shared_ptr<csr> CustomEvaluator::join(std::shared_ptr<csr> &left, std::shared_ptr<csr> &right) {

    auto out = std::make_shared<csr>(left->getNoVertices());
    out->setNoLabels(1);

    for(uint32_t leftSource = 0; leftSource < left->getNoVertices(); leftSource++) {
        // for (auto labelTarget : left->adj[leftSource]) {

        //     int leftTarget = labelTarget.second;
        //     // try to join the left target with right s
        //     for (auto rightLabelTarget : right->adj[leftTarget]) {

        //         auto rightTarget = rightLabelTarget.second;
        //         out->addEdge(leftSource, rightTarget, 0);

        //     }
        // }
        int i=0;
    }

    return out;
}

/**
 * Given an AST, evaluate the query and produce an answer graph.
 * @param q Parsed AST.
 * @return Solution as a graph.
 */
std::shared_ptr<csr> CustomEvaluator::evaluatePath(PathTree *q) {

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
            return CustomEvaluator::selectLabel(label, label, false, graph);
        } else if (std::regex_search(q->data, matches, inverseLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            return CustomEvaluator::selectLabel(label, label, true, graph);
        }
        else if(std::regex_search(q->data, matches, kleeneStar)) {
            label = (uint32_t) std::stoul(matches[1]);
            return CustomEvaluator::transitiveClosure(label, graph);
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
        }

        return nullptr;
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftGraph = CustomEvaluator::evaluatePath(q->left);
        auto rightGraph = CustomEvaluator::evaluatePath(q->right);

        // join left with right
        return CustomEvaluator::join(leftGraph, rightGraph);

    }

    return nullptr;
}

/**
 * Perform a selection on a source constant.
 * @param s A source constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<csr> selectSource(std::string &s, std::shared_ptr<csr> &in) {

    auto out = std::make_shared<csr>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    // for (auto labelTarget : in->adj[std::stoi(s)]) {

    //     auto label = labelTarget.first;
    //     auto target = labelTarget.second;

    //     out->addEdge((uint32_t) std::stoi(s), target, label);
    // }

    return out;
}

/**
 * Perform a selection on a target constant.
 * @param s A target constant.
 * @param in A graph to select from.
 * @return An answer graph as a result of the given selection.
 */
std::shared_ptr<csr> selectTarget(std::string &t, std::shared_ptr<csr> &in) {

    auto out = std::make_shared<csr>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    // for (auto labelSource : in->reverse_adj[std::stoi(t)]) {

    //     auto label = labelSource.first;
    //     auto source = labelSource.second;
    //     out->addEdge(source, (uint32_t) std::stoi(t), label);
    // }

    return out;
}

void CustomEvaluator::findBestPlan(PathQuery *query) {
    
}

/**
 * Evaluate a path query. Produce a cardinality of the answer graph.
 * @param query Query to evaluate.
 * @return A cardinality statistics of the answer graph.
 */
cardStat CustomEvaluator::evaluate(PathQuery *query) {

    /// Find best plan
    findBestPlan(query);

    auto res = evaluatePath(query->path);
    if(query->s != "*") res = selectSource(query->s, res);
    else if(query->t != "*") res = selectTarget(query->t, res);
    return CustomEvaluator::computeStats(res);
}