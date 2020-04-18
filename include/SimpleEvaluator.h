#ifndef QS_SIMPLEEVALUATOR_H
#define QS_SIMPLEEVALUATOR_H


#include <memory>
#include <cmath>
#include "SimpleGraph.h"
#include "PathTree.h"
#include "Evaluator.h"
#include "Graph.h"
#include <map>

struct BestPlan {
    std::string query;
    uint32_t cost;
};

class SimpleEvaluator : public Evaluator {

    std::shared_ptr<SimpleGraph> graph;
    std::shared_ptr<SimpleEstimator> est;
    std::map<std::string, BestPlan> planSpace;
    std::map<std::string, BestPlan> cachedPlans;

public:

    explicit SimpleEvaluator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEvaluator() = default;

    void prepare() override ;
    BestPlan findBestPlan(std::string query);
    uint32_t estimateQueryCost(std::string left, std::string right);

    cardStat evaluate(PathQuery *query) override ;
    void attachEstimator(std::shared_ptr<SimpleEstimator> &e);

    std::vector<std::pair<uint32_t, uint32_t>> evaluatePath(PathTree *q, int source, int target);
    static std::vector<std::pair<uint32_t, uint32_t>> join(std::vector<std::pair<uint32_t, uint32_t>> &left, std::vector<std::pair<uint32_t, uint32_t>> &right);
    static std::vector<std::pair<uint32_t, uint32_t>> transitiveClosure(uint32_t label, std::vector<std::pair<uint32_t, uint32_t>> &in);
    static uint32_t unionDistinct(std::vector<std::pair<uint32_t, uint32_t>> &left, std::vector<std::pair<uint32_t, uint32_t>> &right);

    static cardStat computeStats(std::vector<std::pair<uint32_t, uint32_t>> &g);

};


#endif //QS_SIMPLEEVALUATOR_H