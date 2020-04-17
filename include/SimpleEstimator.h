#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

class Stats {

public:
    uint32_t labels;
    uint32_t vertices;
    // Array of ints with the total counts for each relation
    std::vector<uint32_t> total_relations;
    std::vector<std::vector<std::vector<uint32_t>>>  relation_pairs;
    std::vector<std::vector<std::vector<uint32_t>>>  reverse_relation_pairs;
    // std::vector<std::vector<std::vector<std::pair<uint32_t, uint32_t>>>>  relation_pairs;
    // std::vector<std::vector<std::vector<std::pair<uint32_t, uint32_t>>>>  reverse_relation_pairs;
    // Array of ints with the distinct counts for each relation, e.g (1 0 2) and (1 0 3)
    // would return 1 for source
    std::vector<uint32_t> distinct_source_relations;
    // would return 2 for target
    std::vector<uint32_t> distinct_target_relations;
    // Array of ints with the relation count of each node for each relation,
    // e.g. source_relations_count[0][5] for relation 0 and node 5
    std::vector<std::vector<uint32_t>> source_relations_count;
    std::vector<std::vector<uint32_t>> target_relations_count;
    // matrix[rel_x][rel_y][x_normal][y_normal]{tuples, s, m, f}
    std::vector<std::vector<std::vector<std::vector<std::vector<uint32_t>>>>> multidimensional_matrix;

public:
    Stats() = default;
    Stats(uint32_t noLabels, uint32_t noVertices);

    // void create_stats(std::vector<std::vector<uint32_t>> *positions_adj, std::vector<std::vector<uint32_t>> *positions_adj_reverse,
    //     std::vector<uint32_t> *IA, std::vector<uint32_t> *IA_reverse);
    void create_stats(std::shared_ptr<SimpleGraph> *g);
};

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;
    Stats stats;
    std::shared_ptr<SimpleGraph> SampleTransitiveClosure(int T, float sample);
    std::shared_ptr<SimpleGraph> SampleTransitiveClosure(int T, int source, bool reverse);
    uint32_t get_in(std::vector<uint32_t> relation_info);
    // Data structure to store vector of tuples for each relation
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> relation_vector;

    // Characteristic sets
    std::vector<int> countNumber;
    std::vector<std::vector<int>> countRelations;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(PathQuery *q) override ;
};

#endif //QS_SIMPLEESTIMATOR_H