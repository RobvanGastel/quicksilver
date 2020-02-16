#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;
    std::vector<int> sampleVertices;
//    Data structure to store vector of tuples for each relation
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> relation_vector;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(PathQuery *q) override ;
};


#endif //QS_SIMPLEESTIMATOR_H
