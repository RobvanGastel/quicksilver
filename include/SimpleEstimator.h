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

/// Histogram class
class Histogram {

public:
    uint32_t labels;
    uint32_t vertices;
    uint32_t depth;
    uint32_t width_size;
    uint32_t total_memory;
    uint32_t bucket_memory;
    uint32_t noBuckets;
    uint32_t histogram_type;
    std::vector<uint32_t> total_relations;
    
    std::vector<uint32_t> distinct_source_relations;
    std::vector<uint32_t> distinct_target_relations;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> relation_pairs;
    std::vector<std::vector<uint32_t>> source_relations_count;
    std::vector<std::vector<uint32_t>> target_relations_count;
    std::vector<std::vector<std::vector<uint32_t>>> source_buckets;
    std::vector<std::vector<std::vector<uint32_t>>> target_buckets;

public:
    Histogram(std::string &type_of_histogram, uint32_t noLabels, uint32_t noVertices, uint32_t u_depth);
    ~Histogram();

    void create_histograms(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

    void create_voptimal_histograms();

    void create_equidepth_histograms();

    void create_equiwidth_histograms();

    // void create_voptimal_histograms(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj);

    void create_frequency_vectors(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

    void print_histogram(uint32_t query_var, uint32_t relation);

    uint32_t get_query_results(uint32_t nodeID, uint32_t query_var, uint32_t relation);

};

#endif //QS_SIMPLEESTIMATOR_H
