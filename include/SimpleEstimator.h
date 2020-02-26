#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

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
        // Array of ints with the total counts for each relation
        std::vector<uint32_t> total_relations;
        // Array of ints with the distinct counts for each relation, e.g (1 0 2) and (1 0 3)
        // would return 1 for source
        std::vector<uint32_t> distinct_source_relations;
        // would return 2 for target
        std::vector<uint32_t> distinct_target_relations;
        // Array of tuples for each relation
        std::vector<std::vector<std::vector<std::pair<uint32_t, uint32_t>>>> relation_pairs;
        std::vector<std::vector<std::vector<std::pair<uint32_t, uint32_t>>>> reverse_relation_pairs;
        // Array of ints with the relation count of each node for each relation,
        // e.g. source_relations_count[0][5] for relation 0 and node 5
        std::vector<std::vector<uint32_t>> source_relations_count;
        std::vector<std::vector<uint32_t>> target_relations_count;
        // Array of tuples for each relation with the starting node, ending node and bucket size
        std::vector<std::vector<std::vector<uint32_t>>> source_buckets;
        std::vector<std::vector<std::vector<uint32_t>>> target_buckets;
        std::vector<std::vector<std::vector<std::vector<std::vector<uint32_t>>>>> multidimensional_matrix;

    public:
        Histogram() = default;
        Histogram(std::string &type_of_histogram, uint32_t noLabels, uint32_t noVertices);
        ~Histogram();

        void create_histograms(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

        void create_voptimal_histograms();

        void create_equidepth_histograms();

        void create_equiwidth_histograms();

        void create_frequency_vectors(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

        void print_histogram(uint32_t query_var, uint32_t relation);

        uint32_t get_query_results(uint32_t nodeID, uint32_t query_var, uint32_t relation);
};

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;
    Histogram histogram;
    std::shared_ptr<SimpleGraph> SampleTransitiveClosure(int T, float sample);
    std::shared_ptr<SimpleGraph> SampleTransitiveClosure(int T, int source, bool reverse);
    // Data structure to store vector of tuples for each relation
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> relation_vector;

    std::vector<std::vector<int>> charsets;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(PathQuery *q) override ;
};

#endif //QS_SIMPLEESTIMATOR_H
