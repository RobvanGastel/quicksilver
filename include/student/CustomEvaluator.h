#ifndef QS_CUSTOMEVALUATOR_H
#define QS_CUSTOMEVALUATOR_H


#include <memory>
#include <cmath>
#include "csr.h"
#include "PathTree.h"
#include "Evaluator.h"
#include "Graph.h"

class CustomEvaluator : public Evaluator {

    std::shared_ptr<csr> graph;
    std::shared_ptr<SimpleEstimator> est;

public:

    CustomEvaluator(std::shared_ptr<csr> &g);
    ~CustomEvaluator() = default;

    void prepare() override ;
    void findBestPlan(PathQuery *query);
    cardStat evaluate(PathQuery *query) override ;

    void attachEstimator(std::shared_ptr<SimpleEstimator> &e);

    std::shared_ptr<csr> evaluatePath(PathTree *q);
    static std::shared_ptr<csr> selectLabel(uint32_t projectLabel, uint32_t outLabel, bool inverse, std::shared_ptr<csr> &in);
    static std::shared_ptr<csr> join(std::shared_ptr<csr> &left, std::shared_ptr<csr> &right);
    static std::shared_ptr<csr> transitiveClosure(uint32_t label, std::shared_ptr<csr> &in);
    static uint32_t unionDistinct(std::shared_ptr<csr> &left, std::shared_ptr<csr> &right);

    static cardStat computeStats(std::shared_ptr<csr> &g);

};


#endif //QS_CUSTOMEVALUATOR_H
