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

    std::map<std::string, std::shared_ptr<SimpleGraph>> cachedQuery;

public:

    explicit SimpleEvaluator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEvaluator() = default;

    void prepare() override ;
    BestPlan findBestPlan(std::string query);
    std::string PathQueryBestPlan(PathQuery* query);
    uint32_t estimateQueryCost(std::string left, std::string right);

    cardStat evaluate(PathQuery *query) override ;

    void attachEstimator(std::shared_ptr<SimpleEstimator> &e);

    std::shared_ptr<SimpleGraph> evaluatePath(PathTree *q);
    static std::shared_ptr<SimpleGraph> selectLabel(uint32_t projectLabel, uint32_t outLabel, bool inverse, std::shared_ptr<SimpleGraph> &in);
    static std::shared_ptr<SimpleGraph> join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right);
    static std::shared_ptr<SimpleGraph> transitiveClosure(uint32_t label, std::shared_ptr<SimpleGraph> &in);
    static uint32_t unionDistinct(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right);

    static cardStat computeStats(std::shared_ptr<SimpleGraph> &g);

};


#endif //QS_SIMPLEEVALUATOR_H
