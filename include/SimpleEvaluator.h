#ifndef QS_SIMPLEEVALUATOR_H
#define QS_SIMPLEEVALUATOR_H


#include <memory>
#include <cmath>
#include "SimpleGraph.h"
#include "PathTree.h"
#include "Evaluator.h"
#include "Graph.h"

struct BestPlan {
    int cost;
    PathQuery* plan;

    public:
    BestPlan(int c, PathQuery* p) : plan(p), cost(c) {};
    BestPlan* clone() const { return new BestPlan(*this); }
};

class SimpleEvaluator : public Evaluator {

    std::shared_ptr<SimpleGraph> graph;
    std::shared_ptr<SimpleEstimator> est;

public:

    explicit SimpleEvaluator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEvaluator() = default;

    void prepare() override ;
    BestPlan findBestPlan(BestPlan S);
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
