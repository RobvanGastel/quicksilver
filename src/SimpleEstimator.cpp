#include "SimpleGraph.h"
#include "SimpleEstimator.h"
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
//    noBuckets = total_memory / bucket_memory;
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
        total_relations.push_back({0});
        source_relations_count.push_back({});
        target_relations_count.push_back({});
        distinct_source_relations.push_back({0});
        distinct_target_relations.push_back({0});
        for (int k = 0; k < vertices; k++) {
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
            relation_pairs[rel_type].push_back(std::make_pair(i, rel_target));
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
///// SimpleEstimator class
/////
SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;

}

void SimpleEstimator::prepare() {

    int noLabels = graph->getNoLabels();
    int noVertices = graph->getNoVertices();

    /// Sample tuple sizes
    /// TODO: Find good ratio, sample 5% of the vertices
    int samples = ceil(0.05 * noVertices);

    std::vector<int> sampleCount(noLabels, 0);
    for(int i = 0; i < samples; i++) {
        int index = rand() % noVertices;

        if (graph->adj[index].size() > 0) {
            int adjIndex = rand() % graph->adj[index].size();
            int ct = graph->adj[index][adjIndex].first;
            sampleCount[ct] = sampleCount[ct] + 1;
        } else { // No relation edge found,
            i--;
        }
    }
    // Multiply sample by 20
    for(int i = 0; i < sampleCount.size(); i++) {
        sampleCount[i] = sampleCount[i] * 20;
    }

    sampleVertices = sampleCount;

    /// Creation histograms
    std::string histogram_type = "equiwidth";
    histogram = Histogram(histogram_type, noLabels, noVertices);
    histogram.create_histograms(graph->adj);
    // histogram.print_histogram(0, 0);
    // std::cout << histogram.get_query_results(985, 0, 0) << std::endl;
}


// Calculate the V(R, A) values based on the histogram
int distinctValuesFor(int relation, std::string attribute) {
    // Default is 10
    return 10;
}

/// Parse tree to 2D vector
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
    
    // std::cout << std::endl;
    // std::cout << "Pairs: ";
    for (int i = 0; i < query.size(); i++) {
        // std::cout << query.at(i) << ", ";
    }
    // std::cout << std::endl;
    return query;
}

/// Input param; <index vertex, adjacency list index>
// int findTransitiveClosure(int index, std::vector<int>> adjIndices) {
    // std::unordered_set<int> passedVertices;
    // int tcCount = adjIndices.size();
    // passedVertices.insert(graph->adj[index]);

    // std::queue<int> Q;
    // for (int i = 0; i < adjIndices.size(); i++) {
    //     vertices.push(graph->adj[index][adjIndices[i]].second);
    // }

    // while (!queue.empty()) {
    //     s = queue.front(); 
    //     cout << s << " "; 
    //     queue.pop_front();
    // }
    // return 0;
// }

cardStat SimpleEstimator::estimate(PathQuery *q) {
    int32_t T = -1; /// Current Tuple "Table"
    auto path = parsePathTree(q->path);

    uint32_t noSources = 1;
    uint32_t noPaths = 1;
    uint32_t noTargets = 1;
    
    Either there are no joins (e.g. just 1 relation/table) 
    or it's a transitive closure (TC).
    if (path.size() == 1) {
        T = std::stoi(path[0].substr(0, path[0].size()-1));
        std::string relation = path[0].substr(path[0].size()-1, 1);

        // std::cout << histogram.get_query_results(33,0,0) << std::endl;
        // std::cout << histogram.source_relations_count[0][33] << std::endl;
        // std::cout << histogram.get_query_results(29,1,0) << std::endl;
        // std::cout << histogram.target_relations_count[0][29] << std::endl;

        // // Cases (precise): 
        // if (relation == ">") { // (s,t) such that (s, l, t)
        //     if (q->s == "*") { 
        //         if (q->t =="*") { // - Source: *, Target: *
        //         noSources = histogram.distinct_source_relations[T];
        //         noPaths = histogram.total_relations[T];
        //         noTargets = histogram.distinct_target_relations[T];
        //         } else { // - Source: *, Target: i
        //             int t_i = std::stoi(q->t);
        //             noSources = histogram.target_relations_count[T][t_i];
        //             noPaths = histogram.target_relations_count[T][t_i];
        //             noTargets = 1;                    
        //         }
        //     } else {
        //         int s_i = std::stoi(q->s);

        //         if (q->t =="*") { // - Source: i, Target: *
        //             noSources = 1;
        //             noPaths = histogram.source_relations_count[T][s_i];
        //             noTargets = histogram.source_relations_count[T][s_i];
        //         } else { // - Source: i, Target: i
        //             int t_i = std::stoi(q->t);
        //             int result = std::min(histogram.target_relations_count[T][t_i], histogram.source_relations_count[T][s_i]);
        //             noSources = result;
        //             noPaths = result;
        //             noTargets = result;
        //         }
        //     }
        // } else if(relation == "<") { // (s,t) such that (t, l, s)
        //     if (q->s == "*") { 
        //         if (q->t =="*") { // - Source: *, Target: *
        //         noSources = histogram.distinct_source_relations[T];
        //         noPaths = histogram.total_relations[T];
        //         noTargets = histogram.distinct_target_relations[T];
        //         } else { // - Source: *, Target: i
        //             int t_i = std::stoi(q->t);
        //             noSources = histogram.source_relations_count[T][t_i];
        //             noPaths = histogram.source_relations_count[T][t_i];
        //             noTargets = 1;                    
        //         }
        //     } else {
        //         int s_i = std::stoi(q->s);

        //         if (q->t =="*") { // - Source: i, Target: *
        //             noSources = 1;
        //             noPaths = histogram.target_relations_count[T][s_i];
        //             noTargets = histogram.target_relations_count[T][s_i];
        //         } else { // - Source: i, Target: i
        //             int t_i = std::stoi(q->t);
        //             int result = std::min(histogram.source_relations_count[T][t_i], histogram.target_relations_count[T][s_i]);
        //             noSources = result;
        //             noPaths = result;
        //             noTargets = result;
        //         }
        //     }
        // }
        
        // Cases (estimates): 
        // std::cout << histogram.get_query_results(29,1,0) << std::endl;
        // std::cout << histogram.target_relations_count[0][29] << std::endl;
        if (relation == ">") { // (s,t) such that (s, l, t)
            if (q->s == "*") { 
                if (q->t =="*") { // - Source: *, Target: *
                noSources = histogram.distinct_source_relations[T];
                noPaths = histogram.total_relations[T];
                noTargets = histogram.distinct_target_relations[T];
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    noSources = histogram.target_relations_count[T][t_i];
                    noPaths = histogram.target_relations_count[T][t_i];
                    noTargets = 1;                    
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // - Source: i, Target: *
                    noSources = 1;
                    noPaths = histogram.source_relations_count[T][s_i];
                    noTargets = histogram.source_relations_count[T][s_i];
                } else { // - Source: i, Target: i
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.target_relations_count[T][t_i], histogram.source_relations_count[T][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        } else if(relation == "<") { // (s,t) such that (t, l, s)
            if (q->s == "*") { 
                if (q->t =="*") { // - Source: *, Target: *
                noSources = histogram.distinct_source_relations[T];
                noPaths = histogram.total_relations[T];
                noTargets = histogram.distinct_target_relations[T];
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    noSources = histogram.source_relations_count[T][t_i];
                    noPaths = histogram.source_relations_count[T][t_i];
                    noTargets = 1;                    
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // - Source: i, Target: *
                    noSources = 1;
                    noPaths = histogram.target_relations_count[T][s_i];
                    noTargets = histogram.target_relations_count[T][s_i];
                } else { // - Source: i, Target: i
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.source_relations_count[T][t_i], histogram.target_relations_count[T][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        }
        /// - Source: *, Target: * (TC)
        else if(relation == "+") {
            /// TODO: Paper to improve: 
            /// Estimating Result Size and Execution Times for Graph Queries
            /// Silke Trißl and Ulf Leser

            /// Current approach follow a set of paths:
            // int tupleCount = histogram.total_relations(T);
            // int tcCount = 0;
            // /// 5% sample
            // for (int j = 0; j < ceil(tupleCount * 0,05); j++) {
            //     int index = rand() % graph->getNoVertices();
            //     std::vector<int>> AdjIndices;

            //     for (int i = 0; i < graph->adj[index].size(); i++) {
            //         if (graph->adj[index][i].first == T) {
            //             indices.insert(adjIndex));
            //         }
            //     }
            //     if (indices.size() == 0) {
            //         j--;
            //     }

            //     // Calculate the transitive closure
            //     tcCount += findTransitiveClosure(index, adjIndices);
            // }

            // tcCount = tcCount * 20;
        }
    }

    /// Cases of joins:
    /// Order doesn't matter => s = "*" and t = "*"
    /// Order right to left => s = "*" and t = 1
    /// Order left to right => s = 1 and t = "*", so reverse
    if (q->s != "*") {
        std::reverse(path.begin(), path.end());
    }

    int j = path.size()-2;
    int fullSize = path.size()-1;
    if(path.size() > 1 && q->s == "*" && q->t == "*") { // - Source: *, Target: *
        int Trs;

        while (path.size() > 0) {
            int Tr = std::stoi(path[j].substr(0, path[j].size()-1));
            std::string relation = path[j].substr(path[j].size()-1, 1);

            // std::cout << histogram.total_relations[Tr]/histogram.distinct_target_relations[Tr] << std::endl;

            if (relation == ">" || relation == "<") {
                int Tr_count = histogram.total_relations[Tr];
                int Ts_count;
                int v;

                if (j+1 == fullSize) {
                    int Ts = std::stoi(path[j+1].substr(0, path[j+1].size()-1));
                    Ts_count = histogram.total_relations[Ts];
                    v = std::max(
                        Tr_count / histogram.distinct_source_relations[Tr],
                        Ts_count / histogram.distinct_target_relations[Ts]
                    );
                    std::cout << "v: " << v << std::endl;
                    path.pop_back();
                } else {
                    v = 100;
                    Ts_count = Trs;
                }
                
                int Trs = ceil(Tr_count * Ts_count / v);
                std::cout << Trs << std::endl;
            }

            path.pop_back();
            j--;
        }



    }
    
    noSources = std::max(noSources, (uint32_t)1);
    noPaths = std::max(noPaths, (uint32_t)1);
    noTargets = std::max(noTargets, (uint32_t)1);

 
    return cardStat {noSources, noPaths, noTargets};
}