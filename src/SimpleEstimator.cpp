#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include <set>
#include <cmath>
#include <cfloat>
#include <queue>

/////
///// Histogram class
/////
Histogram::Histogram(std::string &type_of_histogram, uint32_t noLabels,
                     uint32_t noVertices) {
    labels = noLabels;
    vertices = noVertices;
    total_memory = 1000000;
    bucket_memory = 3 * 32;
    noBuckets = 200;

    if (type_of_histogram == "equidepth")
        histogram_type = 0;
    else if (type_of_histogram == "equiwidth")
        histogram_type = 1;
    else if (type_of_histogram == "voptimal")
        histogram_type = 2;
    else {
        // std::cout << "Incorrect type of histogram specified" << std::endl;
        exit (EXIT_FAILURE);
    }
    source_buckets.push_back({});
    total_relations.push_back({});
    distinct_source_relations.push_back({});
    distinct_target_relations.push_back({});
}

Histogram::~Histogram() {
//    TODO: Free all memory
}

void Histogram::create_histograms(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj) {
    create_frequency_vectors(adj);
    if (histogram_type == 0)
        create_equidepth_histograms();
    else if (histogram_type == 1)
        create_equiwidth_histograms();
    else if (histogram_type == 2)
        create_voptimal_histograms();
}

void Histogram::create_equidepth_histograms() {
    for (int i = 0; i < labels; i++) {
        depth = total_relations[i]/ noBuckets;
        uint32_t n = 0;
        source_buckets.push_back({});
        source_buckets[i].push_back({});
        for (int j = 0; j < source_relations_count[i].size(); j++) {
            uint32_t start = j;
            uint32_t end = j;
            uint32_t sum = 0;
            while (sum < depth) {
                if ((sum > 0) && (source_relations_count[i][j] > depth)) {
                    break;
                } else if (j == source_relations_count[i].size()) {
                    break;
                } else {
                    sum += source_relations_count[i][j];
                    end = j;
                    j++;
                }
            }
            source_buckets[i][n].push_back(start);
            source_buckets[i][n].push_back(end);
            source_buckets[i][n].push_back(sum);
            n++;
            if (j < source_relations_count[i].size()) {
                source_buckets[i].push_back({});
                j--;
            }
        }
        n = 0;
        target_buckets.push_back({});
        target_buckets[i].push_back({});
        for (int j = 0; j < target_relations_count[i].size(); j++) {
            uint32_t start = j;
            uint32_t end = j;
            uint32_t sum = 0;
            while (sum < depth) {
                if ((sum > 0) && (target_relations_count[i][j] > depth)) {
                    break;
                } else if (j == target_relations_count[i].size()) {
                    break;
                } else {
                    sum += target_relations_count[i][j];
                    end = j;
                    j++;
                }
            }
            target_buckets[i][n].push_back(start);
            target_buckets[i][n].push_back(end);
            target_buckets[i][n].push_back(sum);
            n++;
            if (j < target_relations_count[i].size()) {
                target_buckets[i].push_back({});
                j--;
            }
        }
    }
}

void Histogram::create_equiwidth_histograms() {
    for (int i = 0; i < labels; i++) {
        width_size = source_relations_count[i].size()/ noBuckets;
        uint32_t n = 0;
        source_buckets.push_back({});
        source_buckets[i].push_back({});
        for (int j = 0; j < source_relations_count[i].size(); j++) {
            uint32_t count_it = 0;
            uint32_t start = j;
            uint32_t sum = 0;
            while (count_it < width_size) {
                if (j == source_relations_count[i].size()) {
                    break;
                } else {
                    sum += source_relations_count[i][j];
                    j++;
                    count_it++;
                }
            }
            source_buckets[i][n].push_back(start);
            source_buckets[i][n].push_back(j - 1);
            source_buckets[i][n].push_back(sum);
            n++;
            if (j < source_relations_count[i].size()) {
                source_buckets[i].push_back({});
                j--;
            }
        }
        n = 0;
        target_buckets.push_back({});
        target_buckets[i].push_back({});
        for (int j = 0; j < target_relations_count[i].size(); j++) {
            uint32_t count_it = 0;
            uint32_t start = j;
            uint32_t sum = 0;
            while (count_it < width_size) {
                count_it++;
                if (j == target_relations_count[i].size()) {
                    break;
                } else {
                    sum += target_relations_count[i][j];
                    j++;
                }
            }
            target_buckets[i][n].push_back(start);
            target_buckets[i][n].push_back(j - 1);
            target_buckets[i][n].push_back(sum);
            n++;
            if (j < target_relations_count[i].size()) {
                target_buckets[i].push_back({});
                j--;
            }
        }
    }
}

void Histogram::create_voptimal_histograms() {
    // std::cout << "Creating V-Optimal Histogram" << std::endl;
    for (int i = 0; i < labels; i++) {
        // std::cout << "Relation " << i << std::endl;
        uint32_t n = 0;
        source_buckets.push_back({});
        for (uint32_t j = 0; j < source_relations_count[i].size(); j++) {
            source_buckets[i].push_back({j, j, source_relations_count[i][j]});
            n++;
        }
        uint32_t v_error = 0;
        while (source_buckets[i].size() > noBuckets) {
            double_t min = DBL_MAX;
            uint32_t k = 0;
            for (uint32_t j = 0; j < source_buckets[i].size() - 1; j++) {
                std::vector<uint32_t> temp_bucket = {source_buckets[i][j][0], source_buckets[i][j+1][1], source_buckets[i][j][2] + source_buckets[i][j+1][2]};
                double_t error1 = source_buckets[i][j][2] * ((double_t)(source_buckets[i][j][1] - source_buckets[i][j][0]) / 12);
                double_t error2 = source_buckets[i][j+1][2] * ((double_t)(source_buckets[i][j+1][1] - source_buckets[i][j+1][0]) / 12);
                double_t error3 = temp_bucket[2] * ((double_t)(temp_bucket[1] - temp_bucket[0]) / 12);
                double_t new_error = v_error - error1 - error2 + error3;
                if (new_error < min) {
                    min = new_error;
                    k = j;
                }
            }
            source_buckets[i][k][2] += source_buckets[i][k+1][2];
            source_buckets[i][k][1] = source_buckets[i][k+1][1];
            source_buckets[i].erase(source_buckets[i].begin() + k + 1);
            v_error = min;
        }
        n = 0;
        target_buckets.push_back({});
        for (uint32_t j = 0; j < target_relations_count[i].size(); j++) {
            target_buckets[i].push_back({j, j, target_relations_count[i][j]});
            n++;
        }
        v_error = 0;
        while (target_buckets[i].size() > noBuckets) {
            double_t min = DBL_MAX;
            uint32_t k = 0;
            for (uint32_t j = 0; j < target_buckets[i].size() - 1; j++) {
                std::vector<uint32_t> temp_bucket = {target_buckets[i][j][0], target_buckets[i][j+1][1], target_buckets[i][j][2] + target_buckets[i][j+1][2]};
                double_t error1 = target_buckets[i][j][2] * ((double_t)(target_buckets[i][j][1] - target_buckets[i][j][0]) / 12);
                double_t error2 = target_buckets[i][j+1][2] * ((double_t)(target_buckets[i][j+1][1] - target_buckets[i][j+1][0]) / 12);
                double_t error3 = temp_bucket[2] * ((double_t)(temp_bucket[1] - temp_bucket[0]) / 12);
                double_t new_error = v_error - error1 - error2 + error3;
                if (new_error < min) {
                    min = new_error;
                    k = j;
                }
            }
            target_buckets[i][k][2] += target_buckets[i][k+1][2];
            target_buckets[i][k][1] = target_buckets[i][k+1][1];
            target_buckets[i].erase(target_buckets[i].begin() + k + 1);
            v_error = min;
        }
    }
}

void Histogram::create_frequency_vectors(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj) {
    for (int i = 0; i < labels; i++) {
        relation_pairs.push_back({});
        reverse_relation_pairs.push_back({});
        total_relations.push_back({0});
        source_relations_count.push_back({});
        target_relations_count.push_back({});
        distinct_source_relations.push_back({0});
        distinct_target_relations.push_back({0});
        for (int k = 0; k < vertices; k++) {
            relation_pairs[i].push_back({});
            reverse_relation_pairs[i].push_back({});
            source_relations_count[i].push_back({0});
            target_relations_count[i].push_back({0});
        }
    }

    for (uint32_t i = 0; i < adj.size(); i++){
        for (uint32_t j = 0; j < adj[i].size() ; j++) {
            uint32_t rel_type = adj[i][j].first;
            uint32_t rel_target = adj[i][j].second;
            source_relations_count[rel_type][i]++;
            target_relations_count[rel_type][rel_target]++;
            total_relations[rel_type]++;
            relation_pairs[rel_type][i].push_back(std::make_pair(i, rel_target));
            reverse_relation_pairs[rel_type][rel_target].push_back(std::make_pair(rel_target, i));
        }
    }
    for (int rel = 0; rel < labels; rel++) {
        for (int i = 0; i < vertices; i++) {
             if (source_relations_count[rel][i] > 0) {
                 distinct_source_relations[rel] += 1;
             }
             if (target_relations_count[rel][i] > 0)
                distinct_target_relations[rel] += 1;
        }
    }

    multidimensional_matrix.push_back({});
    for (int rel_x = 0; rel_x < labels; rel_x++) {
        multidimensional_matrix.push_back({});
        for (int rel_y = 0; rel_y < labels; rel_y++) {
            multidimensional_matrix[rel_x].push_back({});
            for (int x_normal = 0; x_normal < 2; x_normal++) {
                multidimensional_matrix[rel_x][rel_y].push_back({});
                std::vector<std::vector<std::pair<uint32_t, uint32_t>>> x_pairs;
                if (x_normal == 0)
                    x_pairs = relation_pairs[rel_x];
                else
                    x_pairs = reverse_relation_pairs[rel_x];

                for (int y_normal = 0; y_normal < 2; y_normal++) {
                    multidimensional_matrix[rel_x][rel_y][x_normal].push_back({});
                    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> y_pairs;
                    if (y_normal == 0)
                        y_pairs = relation_pairs[rel_y];
                    else
                        y_pairs = reverse_relation_pairs[rel_y];
                    uint32_t tuples = 0;
                    std::vector<uint32_t> source_answers = {};
                    std::vector<uint32_t> middle_answers = {};
                    std::vector<uint32_t> final_answers = {};
                    for (int source_x = 0; source_x < x_pairs.size(); source_x++) {
                        for (int k = 0; k < x_pairs[source_x].size(); k++) {
                            uint32_t target = x_pairs[source_x][k].second;
                            if (y_pairs[target].size() > 0) {
                                source_answers.push_back(source_x);
                                middle_answers.push_back(target);
                            }
                            for (int source_y = 0; source_y < y_pairs[target].size(); source_y++) {
                                uint32_t final_target = y_pairs[target][source_y].second;
                                final_answers.push_back(final_target);
                                tuples += 1;
                            }
                        }
                    }
                    uint32_t source_count = std::distance(source_answers.begin(),
                                                    std::unique(source_answers.begin(), source_answers.end()));
                    uint32_t middle_count = std::distance(middle_answers.begin(),
                                                     std::unique(middle_answers.begin(), middle_answers.end()));
                    uint32_t final_count = std::distance(final_answers.begin(),
                                                    std::unique(final_answers.begin(), final_answers.end()));
                    multidimensional_matrix[rel_x][rel_y][x_normal][y_normal] = {tuples, source_count, middle_count, final_count};
                    //            if (final_count > 0 ){
                    //                std::cout << std::endl;
                    //                std::cout << "For " << rel_x << " and " << rel_y << std::endl;
                    //                std::cout << tuples << std::endl;
                    //                std::cout << middle_count << std::endl;
                    //                std::cout << final_count << std::endl;
                    //            }
                }
            }
        }
    }


//    uint32_t bla = 0;
//    for (int i = 0; i < labels; i++)
//        bla += total_relations[i];
//    printf("Total relations: %d\n", bla);

//    Print size of each relation
//    for (int i = 0; i < labels; i++) {
//        uint32_t relation_sum = 0;
//        for (int k = 0; k < vertices; k++) {
//            relation_sum += source_relations_count[i][k];
//            std::cout << source_relations_count[i][k] << std::endl;
//        }
//        std::cout << "Relation " << i << ": " << relation_sum << std::endl;
//    }
}

void Histogram::print_histogram(uint32_t query_var, uint32_t relation) {
    std::cout << "Histogram " << histogram_type << " for " << query_var << " relation " << relation << std::endl;
    if (query_var == 0) {
        for (int i = 0; i < source_buckets[relation].size(); i++) {
            std::cout << source_buckets[relation][i][0] << "\t" << source_buckets[relation][i][1] << "\t"
                      << source_buckets[relation][i][2] << std::endl;
        }
    }
    else if (query_var == 1) {
        for (int i = 0; i < target_buckets[relation].size(); i++) {
            std::cout << target_buckets[relation][i][0] << "\t" << target_buckets[relation][i][1] << "\t"
                      << target_buckets[relation][i][2] << std::endl;
        }
    }
    std::cout << std::endl;
}

uint32_t Histogram::get_query_results(uint32_t nodeID, uint32_t query_var, uint32_t relation) {
    if (nodeID > vertices)
        return -1;
    int i = 0;
    if (query_var == 0) {
        int noBuckets = source_buckets[relation].size();
        while (nodeID > source_buckets[relation][i][1]) {
            i++;
            if (i == noBuckets)
                return -1;
        }
        return source_buckets[relation][i][2];
    }
    else if (query_var == 1) {
        int noBuckets = target_buckets[relation].size();
        while (nodeID > target_buckets[relation][i][1]) {
            i++;
            if (i == noBuckets)
                return -1;
        }
        return target_buckets[relation][i][2];
    }
    else
        return -1;
}


/////
///// Stats class
/////
Stats::Stats(uint32_t noLabels, uint32_t noVertices) {
    labels = noLabels;
    vertices = noVertices;
    
    total_relations.push_back({});
    source_buckets.push_back({});
    total_relations.push_back({});
    distinct_source_relations.push_back({});
    distinct_target_relations.push_back({});
}

Stats::~Stats() {
    // TODO: Free all memory
}

void Stats::create_stats(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj) {
    relation_pairs = std::vector<uint32_t> (labels);
    reverse_relation_pairs = std::vector<std::vector<uint32_t>> (labels);
    distinct_source_relations = std::vector<uint32_t> (labels);
    distinct_target_relations = std::vector<uint32_t> (labels);

    for (int rel_type = 0; rel_type < labels; rel_type++) {
        total_relations.push_back({0});
        source_relations_count.push_back(std::vector<uint32_t> (vertices));
        target_relations_count.push_back(std::vector<uint32_t> (vertices));
    }

    for (uint32_t rel_source = 0; rel_source < adj.size(); rel_source++){
        for (uint32_t j = 0; j < adj[rel_source].size() ; j++) {
            uint32_t rel_type = adj[rel_source][j].first;
            uint32_t rel_target = adj[rel_source][j].second;
            source_relations_count[rel_type][rel_source]++;
            target_relations_count[rel_type][rel_target]++;
            total_relations[rel_type]++;
        }
    }

    for (uint32_t rel_type = 0; rel_type < labels; rel_type++) {
        for (uint32_t i = 0; i < vertices; i++) {
             if (source_relations_count[rel_type][i] > 0) {
                 distinct_source_relations[rel_type] += 1;
             }
             if (target_relations_count[rel_type][i] > 0)
                distinct_target_relations[rel_type] += 1;
        }
    }
    // for (int i = 0; i < labels; i++) {
    //     relation_pairs.push_back({});
    //     reverse_relation_pairs.push_back({});
    //     total_relations.push_back({0});
    //     source_relations_count.push_back({});
    //     target_relations_count.push_back({});
    //     distinct_source_relations.push_back({0});
    //     distinct_target_relations.push_back({0});
    //     for (int k = 0; k < vertices; k++) {
    //         relation_pairs[i].push_back({});
    //         reverse_relation_pairs[i].push_back({});
    //         source_relations_count[i].push_back({0});
    //         target_relations_count[i].push_back({0});
    //     }
    // }

    for (uint32_t i = 0; i < adj.size(); i++){
        for (uint32_t j = 0; j < adj[i].size() ; j++) {
            uint32_t rel_type = adj[i][j].first;
            uint32_t rel_target = adj[i][j].second;
            source_relations_count[rel_type][i]++;
            target_relations_count[rel_type][rel_target]++;
            total_relations[rel_type]++;
            relation_pairs[rel_type][i].push_back(std::make_pair(i, rel_target));
            reverse_relation_pairs[rel_type][rel_target].push_back(std::make_pair(rel_target, i));
        }
    }
    for (int rel = 0; rel < labels; rel++) {
        for (int i = 0; i < vertices; i++) {
             if (source_relations_count[rel][i] > 0) {
                 distinct_source_relations[rel] += 1;
             }
             if (target_relations_count[rel][i] > 0)
                distinct_target_relations[rel] += 1;
        }
    }


    // LOW MEMORY IMPLEMENTATION
    // std::vector<std::vector<std::vector<std::vector<uint32_t>>>> tuples = std::vector<std::vector<std::vector<std::vector<uint32_t>>>> (labels,
    //     std::vector<std::vector<std::vector<uint32_t>>> (labels,
    //         std::vector<std::vector<uint32_t>> (2,
    //             std::vector<uint32_t> (2, (uint32_t)0))));
    // std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> source_dist = std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> (labels,
    //     std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>> (labels,
    //         std::vector<std::vector<std::unordered_set<uint32_t>>> (2,
    //             std::vector<std::unordered_set<uint32_t>> (2, std::unordered_set<uint32_t> ()))));
    // std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> middle_dist = std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> (labels,
    //     std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>> (labels,
    //         std::vector<std::vector<std::unordered_set<uint32_t>>> (2,
    //             std::vector<std::unordered_set<uint32_t>> (2, std::unordered_set<uint32_t> ()))));
    // std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> final_dist = std::vector<std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>>> (labels,
    //     std::vector<std::vector<std::vector<std::unordered_set<uint32_t>>>> (labels,
    //         std::vector<std::vector<std::unordered_set<uint32_t>>> (2,
    //             std::vector<std::unordered_set<uint32_t>> (2, std::unordered_set<uint32_t> ()))));

    // // matrix[rel_label_i][rel_label_j][rel_type_i][rel_type_j] = {tuples, source_dist, middle_dist, final_dist}
    // multidimensional_matrix = std::vector<std::vector<std::vector<std::vector<std::vector<uint32_t>>>>> (labels,
    //     std::vector<std::vector<std::vector<std::vector<uint32_t>>>> (labels,
    //         std::vector<std::vector<std::vector<uint32_t>>> (2,
    //             std::vector<std::vector<uint32_t>> (2, 
    //                 std::vector<uint32_t> (4, 0)))));

    // for (uint32_t rel_1_source = 0; rel_1_source < adj.size(); rel_1_source++){
    //     for (uint32_t i = 0; i < adj[rel_1_source].size(); i++) {
    //         uint32_t rel_1_type = adj[rel_1_source][i].first;
    //         uint32_t rel_1_target = adj[rel_1_source][i].second;

    //         for (uint32_t rel_2_source = 0; rel_2_source < adj.size(); rel_2_source++){

    //             for (uint32_t j = 0; j < adj[rel_2_source].size(); j++) {
    //                 uint32_t rel_2_type = adj[rel_2_source][j].first;
    //                 uint32_t rel_2_target = adj[rel_2_source][j].second;

    //                 if (rel_1_target == rel_2_source) { // l1>/l2>   l2</l1<
    //                     tuples[rel_1_type][rel_2_type][0][0]++;
    //                 }
    //                 if (rel_1_target == rel_2_target) { // l1>/l2<   l2</l1>
    //                     tuples[rel_1_type][rel_2_type][0][1]++;

    //                 }
    //                 if (rel_1_source == rel_2_source) { // l1</l2>   l2>/l1<
    //                     tuples[rel_1_type][rel_2_type][1][0]++;

    //                 }
    //                 if (rel_1_source == rel_2_target) { // l1</l2<   l2>/l1>
    //                     tuples[rel_1_type][rel_2_type][1][2]++;

    //                 }
    //             }
    //         }
    //     }
    // }

    // std::cout << rel_source << "   " << rel_type << "   " << rel_target << std::endl;

    // multidimensional_matrix.push_back({});
    // for (int rel_x = 0; rel_x < labels; rel_x++) {
    //     multidimensional_matrix.push_back({});
    //     for (int rel_y = 0; rel_y < labels; rel_y++) {
    //         multidimensional_matrix[rel_x].push_back({});
    //         for (int x_normal = 0; x_normal < 2; x_normal++) {
    //             multidimensional_matrix[rel_x][rel_y].push_back({});
    //             std::vector<std::vector<std::pair<uint32_t, uint32_t>>> x_pairs;
    //             if (x_normal == 0)
    //                 x_pairs = relation_pairs[rel_x];
    //             else
    //                 x_pairs = reverse_relation_pairs[rel_x];

    //             for (int y_normal = 0; y_normal < 2; y_normal++) {
    //                 multidimensional_matrix[rel_x][rel_y][x_normal].push_back({});
    //                 std::vector<std::vector<std::pair<uint32_t, uint32_t>>> y_pairs;
    //                 if (y_normal == 0)
    //                     y_pairs = relation_pairs[rel_y];
    //                 else
    //                     y_pairs = reverse_relation_pairs[rel_y];
    //                 uint32_t tuples = 0;
    //                 std::vector<uint32_t> source_answers = {};
    //                 std::vector<uint32_t> middle_answers = {};
    //                 std::vector<uint32_t> final_answers = {};
    //                 for (int source_x = 0; source_x < x_pairs.size(); source_x++) {
    //                     for (int k = 0; k < x_pairs[source_x].size(); k++) {
    //                         uint32_t target = x_pairs[source_x][k].second;
    //                         if (y_pairs[target].size() > 0) {
    //                             source_answers.push_back(source_x);
    //                             middle_answers.push_back(target);
    //                         }
    //                         for (int source_y = 0; source_y < y_pairs[target].size(); source_y++) {
    //                             uint32_t final_target = y_pairs[target][source_y].second;
    //                             final_answers.push_back(final_target);
    //                             tuples += 1;
    //                         }
    //                     }
    //                 }
    //                 uint32_t source_count = std::distance(source_answers.begin(),
    //                                                 std::unique(source_answers.begin(), source_answers.end()));
    //                 uint32_t middle_count = std::distance(middle_answers.begin(),
    //                                                  std::unique(middle_answers.begin(), middle_answers.end()));
    //                 uint32_t final_count = std::distance(final_answers.begin(),
    //                                                 std::unique(final_answers.begin(), final_answers.end()));
    //                 multidimensional_matrix[rel_x][rel_y][x_normal][y_normal] = {tuples, source_count, middle_count, final_count};
    //                 //            if (final_count > 0 ){
    //                 //                std::cout << std::endl;
    //                 //                std::cout << "For " << rel_x << " and " << rel_y << std::endl;
    //                 //                std::cout << tuples << std::endl;
    //                 //                std::cout << middle_count << std::endl;
    //                 //                std::cout << final_count << std::endl;
    //                 //            }
    //             }
    //         }
    //     }
    // }
}


std::vector<uint32_t> get_relation_info(std::string relation) { // path[0]
    std::vector<uint32_t> relation_info;

    // relation label
    uint32_t T = std::stoi(relation.substr(0, relation.size()-1));
    relation_info.push_back(T);

    // relation direction
    std::string dir = relation.substr(relation.size()-1, 1);
    relation_info.push_back(uint32_t(dir != ">") + uint32_t(dir == "*")); // 0:>; 1:<; 2:+;

    return relation_info;
}


uint32_t get_in(std::vector<uint32_t> relation_info, Histogram histogram) {
    if (relation_info[1] == 0) {
        return histogram.distinct_target_relations[relation_info[0]];
    } else if (relation_info[1] == 1) {
        return histogram.distinct_source_relations[relation_info[0]];
    } else {
        return (uint32_t)-1;
    }
}


/////
///// SimpleEstimator class
/////
SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {
//    std::vector<std::vector<uint32_t >> cardinalities;
//    for (int i = 0; i < graph->getNoLabels(); i++) {
//        cardinalities.push_back({});
//        for (int j = 0; j < graph->getNoVertices(); j++)
//            cardinalities[i].push_back({0});
//    }
//    for (uint32_t i = 0; i < graph->reverse_adj.size(); i++){
//        for (uint32_t j = 0; j < graph->reverse_adj[i].size() ; j++) {
//            uint32_t rel_type = graph->reverse_adj[i][j].first;
//            uint32_t rel_source = graph->reverse_adj[i][j].second;
//            cardinalities
//        }
//    }




    int noLabels = graph->getNoLabels();
    int noVertices = graph->getNoVertices();

    /// Creation histograms
    std::string histogram_type = "equiwidth";
    histogram = Histogram(histogram_type, noLabels, noVertices);
    histogram.create_histograms(graph->adj);
    // histogram.print_histogram(0, 0);
    // std::cout << histogram.get_query_results(985, 0, 0) << std::endl;
    stats = Stats(noLabels, noVertices);
    stats.create_stats(graph->adj);
}

/// Parse tree to Vector of queries
void inorderParse(PathTree *node,
        std::vector<std::string> *query) {
    if (node == nullptr) {
        return;
    }
    inorderParse(node->left, query);

    if (node->data != "/") {
        query->push_back(node->data);
    }

    inorderParse(node->right, query);
}

std::vector<std::string> parsePathTree(PathTree *tree) {
    std::vector<std::string> query;

    if (!tree->isLeaf()) {
        inorderParse(tree, &query);
    } else {
        query.push_back(tree->data);
    }
    return query;
}

/// Sample transitive closure queries
std::shared_ptr<SimpleGraph> SimpleEstimator::SampleTransitiveClosure(int T, float sample) {
    auto se = SimpleEvaluator(graph);
    
    int sampleSize = ceil(sample * graph->getNoVertices());
    int numNewAdded = 1;

    /// Create sample graph (TC)
    // Use max upperbound for labels
    auto sampleGraph = std::make_shared<SimpleGraph>(graph->getNoVertices());
    sampleGraph->setNoLabels(sampleSize); 
    
    // Use max upperbound for labels
    auto base = std::make_shared<SimpleGraph>(graph->getNoVertices());
    base->setNoLabels(sampleSize); 

    for(uint32_t source = 0; source < sampleSize; source++) {
        int index = rand() % graph->getNoVertices();

        for (auto labelTarget : graph->adj[index]) {

            auto label = labelTarget.first;
            auto target = labelTarget.second;

            if (label == T) {
                if(source < sampleSize) {
                    sampleGraph->addEdge(source, target, 0);
                    base->addEdge(source, target, 0);
                } else {
                    break;
                }
            }
        }
    }

    while (numNewAdded) {
        auto delta = se.join(sampleGraph, base);
        numNewAdded = se.unionDistinct(sampleGraph, delta);
    }

    return sampleGraph;
}

/// Sample transitive closure for 1 source or target
std::shared_ptr<SimpleGraph> SimpleEstimator::SampleTransitiveClosure(int T, int node, bool reverse) {
    auto se = SimpleEvaluator(graph);    
    int numNewAdded = 1;

    /// Create sample graph (TC)
    auto sampleGraph = std::make_shared<SimpleGraph>(graph->getNoVertices());
    sampleGraph->setNoLabels(1); 

    auto base = std::make_shared<SimpleGraph>(graph->getNoVertices());
    base->setNoLabels(1);
    
    for (auto labelTarget : graph->adj[node]) {

        auto label = labelTarget.first;
        auto target = labelTarget.second;
        sampleGraph->addEdge(0, target, 0);
        base->addEdge(0, target, 0);
    }

    while (numNewAdded) {
        auto delta = se.join(sampleGraph, base);
        numNewAdded = se.unionDistinct(sampleGraph, delta);
    }

    return sampleGraph;
}

cardStat SimpleEstimator::estimate(PathQuery *q) {
    int32_t T = -1; /// Current Tuple "Table"
    auto path = parsePathTree(q->path);

    uint32_t noSources = 1;
    uint32_t noPaths = 1;
    uint32_t noTargets = 1;
    // std::cout << "\n\npath size: " << path.size() << std::endl;

    // Either there are no joins (e.g. just 1 relation/table) 
    // or it's a transitive closure (TC).
    if (path.size() == 1) {
        T = std::stoi(path[0].substr(0, path[0].size()-1));
        std::string relation = path[0].substr(path[0].size()-1, 1);
        
        // Cases (estimates): 
        // std::cout << histogram.get_query_results(29,1,0) << std::endl;
        // std::cout << histogram.target_relations_count[0][29] << std::endl;
        if (relation == ">") { // (s,t) such that (s, l, t)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                    noSources = histogram.distinct_source_relations[T];
                    noPaths = histogram.total_relations[T];
                    noTargets = histogram.distinct_target_relations[T];
                } else { // source: *, target: i
                    int t_i = std::stoi(q->t);
                    int result = histogram.target_relations_count[T][t_i];
                    noSources = result;
                    noPaths = result;
                    noTargets = 1;                   
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = histogram.source_relations_count[T][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.target_relations_count[T][t_i], histogram.source_relations_count[T][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        } else if(relation == "<") { // (s,t) such that (t, l, s)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                noSources = histogram.distinct_target_relations[T];
                noPaths = histogram.total_relations[T];
                noTargets = histogram.distinct_source_relations[T];
                } else { // source: *, target: j
                    int t_i = std::stoi(q->t);
                    int result = histogram.source_relations_count[T][t_i];
                    noSources = result; 
                    noPaths = result; 
                    noTargets = 1;                  
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = histogram.target_relations_count[T][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.source_relations_count[T][t_i], histogram.target_relations_count[T][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        }
            // - Source: *, Target: * (TC)
        else if(relation == "+") {
            /// TODO: Paper to improve, uses GRIPP data structure: 
            /// Estimating Result Size and Execution Times for Graph Queries
            /// Silke Trißl and Ulf Leser
            
            float sample = 0.05;
            if (q->s == "*") { 
                if (q->t =="*") { // - Source: *, Target: *
                    auto out = SampleTransitiveClosure(T, sample);
                    noSources = histogram.distinct_source_relations[T];
                    // To retrieve 100% value estimate
                    noPaths = out->getNoDistinctEdges()*1/sample; 
                    noTargets = histogram.distinct_source_relations[T];
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    auto out = SampleTransitiveClosure(T, t_i, true);

                    noSources = out->getNoDistinctEdges(); 
                    noPaths = out->getNoDistinctEdges(); 
                    noTargets = 1;
                
                }
            } else {
                int s_i = std::stoi(q->s);
                if (q->t =="*") { // - Source: i, Target: *
                    auto out = SampleTransitiveClosure(T, s_i, false);

                    noSources = 1;
                    noPaths = out->getNoDistinctEdges(); 
                    noTargets = out->getNoDistinctEdges();
                } else { // - Source: i, Target: j
                    int t_j = std::stoi(q->t);
                }
            }  
        }
    } else if (path.size() > 1) {
        // if (q->t != "*") {
        //     std::reverse(path.begin(), path.end());
        // }

        
        // std::vector<int> relation_0 = get_relation_info(path[0]);
        std::vector<uint32_t> relation_i;
        std::vector<uint32_t> relation_j;

        // basic info
        int in;     // l1.in
        int d_oi;   // d(o, T_{r/l1})

        // multidimensional matrix
        std::vector<uint32_t> join_stats;
        float T_j;        // |T_{l1/l2}| -> |T_{j-1/j}|
        int d_sj;       // d(s, T_{l1/l2})
        int middle_j;   // l1/l2.middle
        int d_oj;       // d(o, T_{l1/l2})

        // past join iterations
        int T;      // |T_{r/l1/l2}|
        int d_s;    // d(s, T_{r/l1/l2})
        int d_o;    // d(o, T_{r/l1/l2})

        if (q->s == "*") { // source: *
            if (q->t =="*") { // source: *, target: *
                relation_j = get_relation_info(path[0]);
                std::cout << "\n    path_0: " << path[0] << "  relation_0: " << relation_j[0] << " " << relation_j[1] << std::endl;
                if (relation_j[1] == 2) {
                    std::cout << "TC FOUND" << std::endl;
                }

                T = histogram.total_relations[relation_j[0]];
                // if (relation_j[0] == 0) {
                //     d_oi = histogram.distinct_target_relations[relation_j[0]];
                // } else {
                //     d_oi = histogram.distinct_target_relations[relation_j[0]];
                // }

                for (int j = 1; j < path.size(); j++) {
                    relation_i = relation_j;
                    relation_j = get_relation_info(path[j]);
                    std::cout << "    path_j: " << path[j] << "  relation_j: " << relation_j[0] << " " << relation_j[1] << std::endl;

                    in = get_in(relation_i, histogram);
                    std::cout << "        in: " << in << std::endl;

                    join_stats = histogram.multidimensional_matrix[relation_i[0]][relation_j[0]][relation_i[1]][relation_j[1]];
                    T_j = join_stats[0];
                    d_sj = join_stats[1];
                    middle_j = join_stats[2];
                    d_oj = join_stats[3];
                    std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;
                    join_stats = histogram.multidimensional_matrix[relation_j[0]][relation_i[0]][1-relation_j[1]][1-relation_i[1]];
                    T_j = join_stats[0];
                    d_sj = join_stats[1];
                    middle_j = join_stats[2];
                    d_oj = join_stats[3];
                    std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;
                    std::cout << "        T: " << T << "  d_oi: " << d_oi << std::endl;

                    // calculations
                    T_j = ((float)middle_j / (float)in) * (float)T * ((float)T_j / (float)d_sj);
                    std::cout << "        results T_j: " << T_j << std::endl;
                    
                }
            } else { // source: *, target: j
                // should never happen -> reverse path
            }
        } else {
            if (q->t =="*") { // source: i, target: *

            } else { // source: i, target: j

            }
        }
    }
    

    /// Cases of joins:
    /// Order doesn't matter => s = "*" and t = "*"
    /// Order right to left => s = "*" and t = 1
    /// Order left to right => s = 1 and t = "*", so reverse
    // if (q->s != "*") {
    //     std::reverse(path.begin(), path.end());
    // }


    // Causes; segmentation fault
    // int j = path.size()-2;
    // int fullSize = path.size()-1;
    // if(path.size() > 1 && q->s == "*" && q->t == "*") { // - Source: *, Target: *
    //     int Trs;

    //     while (path.size() > 0) {
    //         int Tr = std::stoi(path[j].substr(0, path[j].size()-1));
    //         std::string relation = path[j].substr(path[j].size()-1, 1);

    //         // std::cout << histogram.total_relations[Tr]/histogram.distinct_target_relations[Tr] << std::endl;

    //         if (relation == ">" || relation == "<") {
    //             int Tr_count = histogram.total_relations[Tr];
    //             int Ts_count;
    //             int v;

    //             if (j+1 == fullSize) {
    //                 int Ts = std::stoi(path[j+1].substr(0, path[j+1].size()-1));
    //                 Ts_count = histogram.total_relations[Ts];
    //                 v = std::max(
    //                     Tr_count / histogram.distinct_source_relations[Tr],
    //                     Ts_count / histogram.distinct_target_relations[Ts]
    //                 );
    //                 std::cout << "v: " << v << std::endl;
    //                 path.pop_back();
    //             } else {
    //                 v = 100;
    //                 Ts_count = Trs;
    //             }
                
    //             int Trs = ceil(Tr_count * Ts_count / v);
    //             std::cout << Trs << std::endl;
    //         }

    //         path.pop_back();
    //         j--;
    //     }
    // }
    
    // To prevent 0 predictions
    // noSources = std::max(noSources, (uint32_t)1);
    // noPaths = std::max(noPaths, (uint32_t)1);
    // noTargets = std::max(noTargets, (uint32_t)1);

    return cardStat {noSources, noPaths, noTargets};
}
