#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include <set>
#include <cmath>
#include <cfloat>

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
    // if (histogram_type == 0)
    //     create_equidepth_histograms();
    // else if (histogram_type == 1)
    //     create_equiwidth_histograms();
    // else if (histogram_type == 2)
    //     create_voptimal_histograms();
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
        return target_buckets[relation][i][2]/(target_buckets[relation][i][1]-target_buckets[relation][i][0]);
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
    distinct_source_relations.push_back({});
    distinct_target_relations.push_back({});
}

Stats::~Stats() {
    // TODO: Free all memory
}

void Stats::create_stats(std::vector<std::vector<std::pair<uint32_t, uint32_t>>> adj) {
    total_relations = std::vector<uint32_t> (labels, 0);
    source_relations_count = std::vector<std::vector<uint32_t>> (labels, std::vector<uint32_t> (vertices));
    target_relations_count = std::vector<std::vector<uint32_t>> (labels, std::vector<uint32_t> (vertices));
    relation_pairs = std::vector<std::vector<std::vector<uint32_t>>> (labels,
        std::vector<std::vector<uint32_t>> (vertices,
            std::vector<uint32_t> ()));
    reverse_relation_pairs = std::vector<std::vector<std::vector<uint32_t>>> (labels,
        std::vector<std::vector<uint32_t>> (vertices,
            std::vector<uint32_t> ()));
    distinct_source_relations = std::vector<uint32_t> (labels);
    distinct_target_relations = std::vector<uint32_t> (labels);

    for (uint32_t rel_source = 0; rel_source < adj.size(); rel_source++){
        for (uint32_t i = 0; i < adj[rel_source].size() ; i++) {
            uint32_t rel_type = adj[rel_source][i].first;
            uint32_t rel_target = adj[rel_source][i].second;
            source_relations_count[rel_type][rel_source]++;
            target_relations_count[rel_type][rel_target]++;
            total_relations[rel_type]++;
            relation_pairs[rel_type][rel_source].push_back(rel_target);
            reverse_relation_pairs[rel_type][rel_target].push_back(rel_source);
        }
    }

    for (uint32_t rel_type = 0; rel_type < labels; rel_type++) {
        for (uint32_t vertice = 0; vertice < vertices; vertice++) {
             if (source_relations_count[rel_type][vertice] > 0) {
                 distinct_source_relations[rel_type]++;
             }
             if (target_relations_count[rel_type][vertice] > 0)
                distinct_target_relations[rel_type]++;
        }
    }

    // // matrix[rel_label_i][rel_label_j][rel_type_i][rel_type_j] = {tuples, source_dist, middle_dist, final_dist}
    multidimensional_matrix = std::vector<std::vector<std::vector<std::vector<std::vector<uint32_t>>>>> (labels,
        std::vector<std::vector<std::vector<std::vector<uint32_t>>>> (labels,
            std::vector<std::vector<std::vector<uint32_t>>> (2,
                std::vector<std::vector<uint32_t>> (2, 
                    std::vector<uint32_t> (4, 0)))));
    std::vector<std::vector<uint32_t>> x_pairs;
    std::vector<std::vector<uint32_t>> y_pairs;
    
    for (uint32_t rel_x = 0; rel_x < labels; rel_x++) {
        for (uint32_t rel_y = 0; rel_y < rel_x; rel_y++) { // < labels
            for (uint32_t x_normal = 0; x_normal < (uint32_t)2; x_normal++) {
                if (x_normal == (uint32_t)0)
                    x_pairs = relation_pairs[rel_x];
                else
                    x_pairs = reverse_relation_pairs[rel_x];
                
                for (uint32_t y_normal = 0; y_normal < (uint32_t)2; y_normal++) {
                    uint32_t tuples = 0;
                    std::unordered_set<uint32_t> source_answers = {};
                    std::unordered_set<uint32_t> middle_answers = {};
                    std::unordered_set<uint32_t> final_answers = {};
                    if (y_normal == (uint32_t)0)
                        y_pairs = relation_pairs[rel_y];
                    else
                        y_pairs = reverse_relation_pairs[rel_y];

                    for (uint32_t source_x = 0; source_x < (uint32_t)x_pairs.size(); source_x++) {
                        for (uint32_t i = 0; i < x_pairs[source_x].size(); i++) {
                            uint32_t x_target = x_pairs[source_x][i];

                            if (y_pairs[x_target].size() > (uint32_t)0) {
                                source_answers.insert(source_x);
                                middle_answers.insert(x_target);
                                final_answers.insert(y_pairs[x_target].begin(), y_pairs[x_target].end());
                                tuples += y_pairs[x_target].size();
                            }
                        }
                    }
                    // std::cout << rel_x << "  " << 
                    //     rel_y << "  " << 
                    //     x_normal << "  " <<
                    //     tuples << "  " << 
                    //     source_answers.size() << "  " << 
                    //     middle_answers.size() << "  " << 
                    //     final_answers.size() << std::endl;

                    multidimensional_matrix[rel_x][rel_y][x_normal][y_normal] = {
                        tuples,
                        (uint32_t)source_answers.size(),
                        (uint32_t)middle_answers.size(), 
                        (uint32_t)final_answers.size()
                    };
                    multidimensional_matrix[rel_y][rel_x][1-y_normal][1-x_normal] = {
                        tuples,
                        (uint32_t)final_answers.size(),
                        (uint32_t)middle_answers.size(),
                        (uint32_t)source_answers.size()
                    };     
                }              
            }
        }

        /// for rel_y == rel_x

        /// rel_x 0 / rel_x 0 && rel_x 1 / rel_x 1
        uint32_t tuples = 0;
        std::unordered_set<uint32_t> source_answers = {};
        std::unordered_set<uint32_t> middle_answers = {};
        std::unordered_set<uint32_t> final_answers = {};
        x_pairs = relation_pairs[rel_x];
        for (uint32_t source_x = 0; source_x < (uint32_t)x_pairs.size(); source_x++) {
            for (uint32_t i = 0; i < x_pairs[source_x].size(); i++) {
                uint32_t x_target = x_pairs[source_x][i];

                if (x_pairs[x_target].size() > (uint32_t)0) {
                    source_answers.insert(source_x);
                    middle_answers.insert(x_target);
                    tuples += x_pairs[x_target].size();
                }
            }
        }
        uint32_t source_answers_size = source_answers.size();
        std::vector<uint32_t> result = {
            tuples,
            source_answers_size,
            (uint32_t)middle_answers.size(), 
            source_answers_size
        };
        multidimensional_matrix[rel_x][rel_x][0][0] = result;
        multidimensional_matrix[rel_x][rel_x][1][1] = result;

        // rel_x 0 / rel_x 1 && rel_x 1 / rel_x 0
        for (uint32_t x_normal = 0; x_normal < (uint32_t)2; x_normal++) {
            if (x_normal == (uint32_t)0) {
                x_pairs = relation_pairs[rel_x];
                y_pairs = reverse_relation_pairs[rel_x];
            }
            else {
                x_pairs = reverse_relation_pairs[rel_x];
                y_pairs = relation_pairs[rel_x];
            }
            
            uint32_t tuples = 0;
            std::unordered_set<uint32_t> source_answers = {};
            std::unordered_set<uint32_t> middle_answers = {};
            std::unordered_set<uint32_t> final_answers = {};

            for (uint32_t source_x = 0; source_x < (uint32_t)x_pairs.size(); source_x++) {
                for (uint32_t i = 0; i < x_pairs[source_x].size(); i++) {
                    uint32_t x_target = x_pairs[source_x][i];

                    if (y_pairs[x_target].size() > (uint32_t)0) {
                        source_answers.insert(source_x);
                        middle_answers.insert(x_target);
                        final_answers.insert(y_pairs[x_target].begin(), y_pairs[x_target].end());
                        tuples += y_pairs[x_target].size();
                    }
                }
            }

            multidimensional_matrix[rel_x][rel_x][x_normal][1-x_normal] = {
                tuples,
                (uint32_t)source_answers.size(),
                (uint32_t)middle_answers.size(), 
                (uint32_t)final_answers.size()
            };
        }
    }    
    // std::cout << "|T|=" << multidimensional_matrix[3][2][1][0][0] <<
    //     "   s=" << multidimensional_matrix[3][2][1][0][1] <<
    //     "   m=" << multidimensional_matrix[3][2][1][0][2] <<
    //     "   o=" << multidimensional_matrix[3][2][1][0][3] << std::endl;
    // std::cout << "|T|=" << multidimensional_matrix[2][3][1][0][0] <<
    //     "   s=" << multidimensional_matrix[2][3][1][0][1] <<
    //     "   m=" << multidimensional_matrix[2][3][1][0][2] <<
    //     "   o=" << multidimensional_matrix[2][3][1][0][3] << std::endl;
    // std::cout << "|T|=" << multidimensional_matrix[1][1][1][0][0] <<
    //     "   s=" << multidimensional_matrix[1][1][1][0][1] <<
    //     "   m=" << multidimensional_matrix[1][1][1][0][2] <<
    //     "   o=" << multidimensional_matrix[1][1][1][0][3] << std::endl;
    // std::cout << "|T|=" << multidimensional_matrix[1][1][0][1][0] <<
    //     "   s=" << multidimensional_matrix[1][1][0][1][1] <<
    //     "   m=" << multidimensional_matrix[1][1][0][1][2] <<
    //     "   o=" << multidimensional_matrix[1][1][0][1][3] << std::endl;
    // std::cout << "|T|=" << multidimensional_matrix[2][2][0][0][0] <<
    //     "   s=" << multidimensional_matrix[2][2][0][0][1] <<
    //     "   m=" << multidimensional_matrix[2][2][0][0][2] <<
    //     "   o=" << multidimensional_matrix[2][2][0][0][3] << std::endl;
    // std::cout << "|T|=" << multidimensional_matrix[2][2][1][1][0] <<
    //     "   s=" << multidimensional_matrix[2][2][1][1][1] <<
    //     "   m=" << multidimensional_matrix[2][2][1][1][2] <<
    //     "   o=" << multidimensional_matrix[2][2][1][1][3] << std::endl;

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


uint32_t get_in(std::vector<uint32_t> relation_info, Stats stats) {
    if (relation_info[1] == 0) {
        return stats.distinct_target_relations[relation_info[0]];
    } else if (relation_info[1] == 1) {
        return stats.distinct_source_relations[relation_info[0]];
    } else {
        return (uint32_t)-1;
    }
}

std::vector<std::__cxx11::string> reverse_path(std::vector<std::__cxx11::string> path) {
    std::vector<std::__cxx11::string> newPath;
    for (int i = path.size()-1; i >= 0 ; i--) {
        if (path[i].substr(1, 2) == ">")
            newPath.push_back(path[i].substr(0, 1)+(std::__cxx11::string)"<");
        else
            newPath.push_back(path[i].substr(0, 1)+(std::__cxx11::string)">");
    }
    return newPath;
}


/////
///// SimpleEstimator class
/////
SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {
    int noLabels = graph->getNoLabels();
    int noVertices = graph->getNoVertices();

    /// Creation histograms
    // std::string histogram_type = "equidepth";
    // histogram = Histogram(histogram_type, noLabels, noVertices);
    // histogram.create_histograms(graph->adj);
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
    
    if (reverse) {
        for (auto labelTarget : graph->reverse_adj[node]) {

            auto label = labelTarget.first;
            auto target = labelTarget.second;
            sampleGraph->addEdge(0, target, 0);
            base->addEdge(0, target, 0);
        }
    } else {
        for (auto labelTarget : graph->adj[node]) {

            auto label = labelTarget.first;
            auto target = labelTarget.second;
            sampleGraph->addEdge(0, target, 0);
            base->addEdge(0, target, 0);
        }    
    }

    while (numNewAdded) {
        auto delta = se.join(sampleGraph, base);
        numNewAdded = se.unionDistinct(sampleGraph, delta);
    }

    return sampleGraph;
}

cardStat SimpleEstimator::estimate(PathQuery *q) {
    int32_t rel_type = -1; /// Current Tuple "Table"
    auto path = parsePathTree(q->path);

    uint32_t noSources = 1;
    uint32_t noPaths = 1;
    uint32_t noTargets = 1;
    // std::cout << "\n\npath size: " << path.size() << std::endl;

    // Either there are no joins (e.g. just 1 relation/table) 
    // or it's a transitive closure (TC).
    if (path.size() == 1) {
        rel_type = std::stoi(path[0].substr(0, path[0].size()-1));
        std::string relation = path[0].substr(path[0].size()-1, 1);
        
        if (relation == ">") { // (s,t) such that (s, l, t)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                    noSources = histogram.distinct_source_relations[rel_type];
                    noPaths = histogram.total_relations[rel_type];
                    noTargets = histogram.distinct_target_relations[rel_type];
                } else { // source: *, target: i
                    int t_i = std::stoi(q->t);
                    int result = histogram.target_relations_count[rel_type][t_i];
                    noSources = result;
                    noPaths = result;
                    noTargets = 1;                   
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = histogram.source_relations_count[rel_type][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.target_relations_count[rel_type][t_i], histogram.source_relations_count[rel_type][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        } else if(relation == "<") { // (s,t) such that (t, l, s)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                noSources = histogram.distinct_target_relations[rel_type];
                noPaths = histogram.total_relations[rel_type];
                noTargets = histogram.distinct_source_relations[rel_type];
                } else { // source: *, target: j
                    int t_i = std::stoi(q->t);
                    int result = histogram.source_relations_count[rel_type][t_i];
                    noSources = result; 
                    noPaths = result; 
                    noTargets = 1;                  
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = histogram.target_relations_count[rel_type][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.source_relations_count[rel_type][t_i], histogram.target_relations_count[rel_type][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        }
            // - Source: *, Target: * (TC)
        else if(relation == "+") {
            
            float sample = 0.05;
            if (q->s == "*") { 
                if (q->t =="*") { // - Source: *, Target: *
                    noSources = stats.distinct_source_relations[rel_type];
                    noTargets = stats.distinct_target_relations[rel_type];
                    noPaths = stats.total_relations[rel_type] + stats.multidimensional_matrix[rel_type][rel_type][0][0][0];
                    // auto out = SampleTransitiveClosure(rel_type, sample);
                    
                    // int count = 0;
                    // for (int i = 0; i < out->adj.size(); i++) {
                    //     if(out->adj[i].size() > 0) {
                    //         count += 1;
                    //     }
                    // }
                    // noSources = count*1/sample;
                    // // To retrieve 100% value estimate
                    // noPaths = out->getNoDistinctEdges()*1/sample; 
                    // noTargets = count*1/sample;
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    auto out = SampleTransitiveClosure(rel_type, t_i, true);

                    noSources = out->getNoDistinctEdges(); 
                    noPaths = out->getNoDistinctEdges(); 
                    noTargets = 1;
                }
            } else {
                int s_i = std::stoi(q->s);
                if (q->t =="*") { // - Source: i, Target: *
                    auto out = SampleTransitiveClosure(rel_type, s_i, false);

                    noSources = 1;
                    noPaths = out->getNoDistinctEdges(); 
                    noTargets = out->getNoDistinctEdges();
                } else { // - Source: i, Target: j
                    /// TODO: Implement
                    int t_j = std::stoi(q->t);
                }
            }  
        }
    } else if(path.size() > 1) {

        /// Cases of joins:
        /// Order doesn't matter => s = "*" and t = "*"
        /// Order right to left => s = "*" and t = 1
        /// Order left to right => s = 1 and t = "*", so reverse
        if (q->t != "*") {
            // std::reverse(path.begin(), path.end());
            path = reverse_path(path);
        }        


        if ((q->s == "*") && (q->t == "*")) { // source: *
            if (q->t == "*") { // source: *, target: *
                std::vector<uint32_t> relation_i;
                std::vector<uint32_t> relation_j;

                // basic info
                float in;     // l1.in
                float T_i;
                float d_si;
                // uint32_t middle_i;
                float d_oi;   // d(o, T_{r/l1})

                float part1;

                // multidimensional matrix
                std::vector<uint32_t> join_stats;
                float T_j;        // |T_{l1/l2}| -> |T_{j-1/j}|
                float d_sj;       // d(s, T_{l1/l2})
                float middle_j;   // l1/l2.middle
                float d_oj;       // d(o, T_{l1/l2})

                // past join iterations
                // uint32_t T;      // |T_{r/l1/l2}|
                // uint32_t d_s;    // d(s, T_{r/l1/l2})
                // uint32_t d_o;    // d(o, T_{r/l1/l2})

                relation_i = get_relation_info(path[0]);
                relation_j = get_relation_info(path[1]);
                // std::cout << "\n    path_0: " << path[0] << "  relation_0: " << relation_i[0] << " " << relation_i[1] << std::endl;
                if ((relation_i[1] == 2) || (relation_j[1] == 2))
                    std::cout << "TC FOUND!" << std::endl;
                
                // T = stats.total_relations[relation_i[0]];
                // if (relation_j[0] == 0) {
                //     d_oi = histogram.distinct_target_relations[relation_j[0]];
                // } else {
                //     d_oi = histogram.distinct_target_relations[relation_j[0]];
                // }

                join_stats = stats.multidimensional_matrix[relation_i[0]][relation_j[0]][relation_i[1]][relation_j[1]];
                T_i = join_stats[0];
                d_si = join_stats[1];
                // int middle_i = join_stats[2];
                d_oi = join_stats[3];
                // std::cout << "        T_i: " << T_i << "  d_si: " << d_si << "  middle_i: " << middle_i << "  d_oi: " << d_oi << std::endl;
                
                for (int j = 2; j < path.size(); j++) {
                    std::cout << "\n        results d_si: " << d_si;
                    std::cout << "        results T_i: " << T_i;
                    std::cout << "        results d_oi: " << d_oi << std::endl;
                    relation_i = relation_j;
                    relation_j = get_relation_info(path[j]);
                    // std::cout << "    relation_j: " << relation_j[0] << " " << relation_j[1] << std::endl;
                    // std::cout << "    path_j: " << path[j] << "  relation_j: " << relation_j[0] << " " << relation_j[1] << std::endl;

                    in = get_in(relation_i, stats);
                    // std::cout << "        in: " << in << std::endl;

                    join_stats = stats.multidimensional_matrix[relation_i[0]][relation_j[0]][relation_i[1]][relation_j[1]];
                    T_j = join_stats[0];
                    d_sj = join_stats[1];
                    middle_j = join_stats[2];
                    d_oj = join_stats[3];
                    std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;

                    // calculations
                    part1 = middle_j / in;
                    // part1 = middle_j / d_oi;
                    std::cout << "        results part1: " << part1 << "  " << in << "  " << d_oi;
                    d_si = d_si * part1;
                    std::cout << "  " << middle_j << "  " << in << std::endl;
                    T_i = T_i * part1 * (T_j / d_sj)/4;
                    d_oi = d_oi * d_oj / in;
                }

                // middle_j = join_stats[2];
                // part1 = (float)middle_j / (float)stats.distinct_target_relations[relation_i[0]];
                // d_si = (float)stats.distinct_target_relations[relation_i[0]] * part1;
                // T_i = (float)stats.total_relations[relation_i[0]] * (float)part1 * ((float)stats.total_relations[relation_j[0]]/(float)stats.distinct_source_relations[relation_j[1]]);
                // d_oi = (float)stats.distinct_source_relations[relation_j[0]] * (float)middle_j / (float)stats.distinct_target_relations[relation_i[0]];
                // std::cout << "        " << part1 << "  T_i: " << T_i << "  d_si: " << d_si << "  d_oi: " << d_oi << std::endl;

                noSources = d_si;
                noPaths = T_i;
                noTargets = d_oi;
            } else { // source: *, target: j
                // should never happen -> reverse path
                std::cout << "WRONG: source: *, target: j" << std::endl;
            }
        } else {
            if ((q->t == "*") || q->s == "*") { // source: i, target: *
                std::cout << "source: i, target: *" << std::endl;
                uint32_t source;
                // basic info
                float in;     // l1.in
                float T_i;
                float d_si;
                // uint32_t middle_i;
                float d_oi;   // d(o, T_{r/l1})

                std::vector<uint32_t> relation_i;
                std::vector<uint32_t> relation_j;
                float part1;

                relation_i = get_relation_info(path[0]);
                if ((q->t == "*")) {
                    source = std::stoi(q->s);
                    if (relation_i[0] == 0)
                        T_i = stats.source_relations_count[relation_i[0]][source];
                    else
                        T_i = stats.target_relations_count[relation_i[0]][source];
                }
                else {
                    source = std::stoi(q->t);
                    if (relation_i[0] == 0)
                        T_i = stats.source_relations_count[relation_i[0]][source];
                    else
                        T_i = stats.target_relations_count[relation_i[0]][source];
                }


                // multidimensional matrix
                std::vector<uint32_t> join_stats;
                float T_j;        // |T_{l1/l2}| -> |T_{j-1/j}|
                float d_sj;       // d(s, T_{l1/l2})
                float middle_j;   // l1/l2.middle
                float d_oj;       // d(o, T_{l1/l2})
                
                
                if ((relation_i[1] == 2))
                    std::cout << "TC FOUND!" << std::endl;
                
                d_oi = T_i;
                std::cout << "    source:" << source << " T_i:" << T_i << std::endl;
                
                for (int j = 1; j < path.size(); j++) {
                    // std::cout << "\n        results d_si: " << d_si;
                    // std::cout << "        results T_i: " << T_i;
                    // std::cout << "        results d_oi: " << d_oi << std::endl;
                    relation_j = get_relation_info(path[j]);
                    // std::cout << "    relation_j: " << relation_j[0] << " " << relation_j[1] << std::endl;
                    // std::cout << "    path_j: " << path[j] << "  relation_j: " << relation_j[0] << " " << relation_j[1] << std::endl;

                    in = get_in(relation_i, stats);
                    // std::cout << "        in: " << in << std::endl;

                    join_stats = stats.multidimensional_matrix[relation_i[0]][relation_j[0]][relation_i[1]][relation_j[1]];
                    T_j = join_stats[0];
                    d_sj = join_stats[1];
                    middle_j = join_stats[2];
                    d_oj = join_stats[3];
                    std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;

                    // calculations
                    part1 = middle_j / in;
                    d_si = d_si * part1;
                    std::cout << "  " << middle_j << "  " << in << std::endl;
                    T_i = T_i * part1 * (T_j / d_sj);
                    d_oi = d_oi * d_oj / in;


                    relation_i = relation_j;
                    T_i = (T_i+d_oi)/2;
                    d_oi = T_i;
                } 

                // middle_j = join_stats[2];
                // part1 = (float)middle_j / (float)stats.distinct_target_relations[relation_i[0]];
                // d_si = (float)stats.distinct_target_relations[relation_i[0]] * part1;
                // T_i = (float)stats.total_relations[relation_i[0]] * (float)part1 * ((float)stats.total_relations[relation_j[0]]/(float)stats.distinct_source_relations[relation_j[1]]);
                // d_oi = (float)stats.distinct_source_relations[relation_j[0]] * (float)middle_j / (float)stats.distinct_target_relations[relation_i[0]];
                // std::cout << "        " << part1 << "  T_i: " << T_i << "  d_si: " << d_si << "  d_oi: " << d_oi << std::endl;
                
                noPaths = T_i;
                if ((q->t == "*")) {
                    noSources = T_i > 0;
                    noTargets = d_oi;
                }
                else {
                    noSources = d_oi;
                    noTargets = T_i > 0;
                }

            } else { // source: i, target: j

            }
        }
    
    

        /// Cases of joins:
        /// Order doesn't matter => s = "*" and t = "*"
        /// Order right to left => s = "*" and t = 1
        /// Order left to right => s = 1 and t = "*", so reverse
        // if (q->s != "*") {
        //     std::reverse(path.begin(), path.end());
        // }

        /// Source and target tuple/"Table" 
        // int T_s = std::stoi(path[0].substr(0, path[0].size()-1));
        // int T_t = std::stoi(path[path.size()-1].substr(0, path[0].size()-1));
        // float card = 1;

        // bool containsTC = false;
        // for(int i = 0; i < path.size(); i++) {
        //     std::string relation = path[i].substr(path[i].size()-1, 1);
        //     if(relation == "+") {
        //         containsTC = true;
        //         break;
        //     }
        // }

        // if (!containsTC) {
        //     if (q->s == "*") {
        //         if (q->t =="*") { // - Source: *, Target: *
        //             /// Basic approach to joins from the characteristic set paper,
        //             /// Of a star join.
        //             int T = 0;
        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
        //             }
        //             card = card * graph->getNoVertices();
                    
        //             noSources = histogram.distinct_target_relations[T_s];
        //             noPaths = card;
        //             noTargets = histogram.distinct_source_relations[T_t];
                
        //         } else { // - Source: *, Target: i
        //             int t_i = std::stoi(q->t);

        //             int T = 0;
        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
        //             }
        //             card = card * graph->getNoVertices();

        //             noSources = card; 
        //             noPaths = card; 
        //             noTargets = 1;                  
        //         }
        //     } else {
        //         int s_i = std::stoi(q->s);
        //         if (q->t =="*") { // - Source: i, Target: *
        //             int T = 0;
        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
        //             }
        //             card = card * graph->getNoVertices();

        //             noSources = 1; 
        //             noPaths = card; 
        //             noTargets = card;

        //         } else { // - Source: i, Target: j
        //             /// TODO: Implement
        //             int t_i = std::stoi(q->t);
        //             int result = std::min(histogram.source_relations_count[T_s][t_i], 
        //                 histogram.target_relations_count[T_t][s_i]);
        //             noSources = result;
        //             noPaths = result;
        //             noTargets = result;
        //         }
        //     }
        // } else { /// The path contains TC
        //     float sample = 0.05;

        //     if (q->s == "*") { 
        //         if (q->t =="*") { // - Source: *, Target: *
        //             float result = 1;
        //             int count = 0;

        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 auto out = SampleTransitiveClosure(T, sample);

        //                 if(i == 0) {
        //                     for (int i = 0; i < out->adj.size(); i++) {
        //                        if(out->adj[i].size() > 0) {
        //                             count += 1;
        //                         }
        //                     }
        //                 }

        //                 result = out->getNoDistinctEdges()*1/sample;
        //                 if(i != path.size()-1) {
        //                     card = card * (float)result/(float)graph->getNoVertices();
        //                 }
        //             }
        //             card = card * result;

        //             if(card != 0) {
        //                 noSources = count*1/sample;
        //                 noPaths = card;
        //                 noTargets = card;
        //             } else {
        //                 noSources = 0;
        //                 noPaths = 0;
        //                 noTargets = 0;
        //             }
                    
        //         } else { // - Source: *, Target: i
        //             int t_i = std::stoi(q->t);                    
        //             float result = 1;
        //             int count = 0;

        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 auto out =  SampleTransitiveClosure(T, t_i, true);
                        
        //                 if(i == 0) {
        //                     for (int i = 0; i < out->adj.size(); i++) {
        //                        if(out->adj[i].size() > 0) {
        //                             count += 1;
        //                         }
        //                     }
        //                 }

        //                 result = out->getNoDistinctEdges()*1/sample;
        //                 if(i != path.size()-1) {
        //                     card = card * (float)result/(float)graph->getNoVertices();
        //                 }
        //             }
        //             card = card * result;

        //             noSources = count*1/sample; 
        //             noPaths = card; 
        //             noTargets = 1;    
        //         }
        //     } else { // - Source: i
        //         int s_i = std::stoi(q->s);
        //         if (q->t =="*") { // - Source: i, Target: *             
        //             float result = 1;
        //             int count = 0;

        //             for (int i = 0; i < path.size(); i++) {
        //                 T = std::stoi(path[i].substr(0, path[i].size()-1));
        //                 auto out =  SampleTransitiveClosure(T, s_i, false);
                        
        //                 if(i == 0) {
        //                     for (int i = 0; i < out->adj.size(); i++) {
        //                        if(out->adj[i].size() > 0) {
        //                             count += 1;
        //                         }
        //                     }
        //                 }

        //                 result = out->getNoDistinctEdges()*1/sample;
        //                 if(i != path.size()-1) {
        //                     card = card * (float)result/(float)graph->getNoVertices();
        //                 }
        //             }
        //             card = card * result;

        //             noSources = 1; 
        //             noPaths = card; 
        //             noTargets = count*1/sample;  

        //         } else { // - Source: i, Target: j
        //             int t_j = std::stoi(q->t);
        //         }
        //     }  
        // }
    }

    return cardStat {noSources, noPaths, noTargets};
}
