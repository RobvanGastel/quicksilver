#ifndef QUICKSILVER_HISTOGRAM_H
#define QUICKSILVER_HISTOGRAM_H

#include <iostream>
#include <vector>

class Histogram {

private:
    uint32_t labels;
    uint32_t vertices;
    uint32_t depth;
    uint32_t width_size;
    std::string histogram_type;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> relations;
    std::vector<std::vector<uint32_t>> source_relations_count;
    std::vector<std::vector<uint32_t>> target_relations_count;
    std::vector<std::vector<std::vector<uint32_t>>> source_buckets;
    std::vector<std::vector<std::vector<uint32_t>>> target_buckets;

public:
    Histogram(std::string &type_of_histogram, uint32_t noLabels, uint32_t noVertices, uint32_t u_depth, uint32_t u_width_size);
    ~Histogram();

    void create_histograms(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

    void create_equidepth_histograms();

    void create_equiwidth_histograms();

    void create_voptimal_histograms();

    void create_frequency_vectors(std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj);

    void print_histogram(std::string query_var, uint32_t relation);

    uint32_t get_query_results(uint32_t nodeID, std::string query_var, uint32_t relation);

};


#endif //QUICKSILVER_HISTOGRAM_H
