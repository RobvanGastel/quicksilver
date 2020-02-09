#include "SimpleGraph.h"
#include "SimpleEstimator.h"

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // do your prep here

}

cardStat SimpleEstimator::estimate(PathQuery *q) {

    // perform your estimation here

    return cardStat {0, 0, 0};
}