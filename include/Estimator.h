#ifndef QS_ESTIMATOR_H
#define QS_ESTIMATOR_H

#include "PathTree.h"
#include "PathQuery.h"
#include <iostream>

struct cardStat {
    uint32_t noOut;
    uint32_t noPaths;
    uint32_t noIn;

    void print() {
        std::cout << "(" << noOut << ", " << noPaths << ", " << noIn << ")" << std::endl;
    }
};

class Estimator {

public:

    virtual void prepare() = 0;
    virtual cardStat estimate(PathQuery *q) = 0;

};


#endif //QS_ESTIMATOR_H
