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
                    std::vector<uint32_t> middle_answers = {};
                    std::vector<uint32_t> final_answers = {};
                    for (int source_x = 0; source_x < x_pairs.size(); source_x++) {
                        for (int k = 0; k < x_pairs[source_x].size(); k++) {
                            uint32_t target = x_pairs[source_x][k].second;
                            if (y_pairs[target].size() > 0)
                                middle_answers.push_back(target);
                            for (int source_y = 0; source_y < y_pairs[target].size(); source_y++) {
                                uint32_t final_target = y_pairs[target][source_y].second;
                                final_answers.push_back(final_target);
                                tuples += 1;
                            }
                        }
                    }
                    uint32_t middle_count = std::distance(middle_answers.begin(),
                                                     std::unique(middle_answers.begin(), middle_answers.end()));
                    uint32_t final_count = std::distance(final_answers.begin(),
                                                    std::unique(final_answers.begin(), final_answers.end()));
                    multidimensional_matrix[rel_x][rel_y][x_normal].push_back({tuples, middle_count, final_count});
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
    std::string histogram_type = "equidepth";
    histogram = Histogram(histogram_type, noLabels, noVertices);
    histogram.create_histograms(graph->adj);
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
    int32_t T = -1; /// Current Tuple "Table"
    auto path = parsePathTree(q->path);

    uint32_t noSources = 1;
    uint32_t noPaths = 1;
    uint32_t noTargets = 1;

    // Either there are no joins (e.g. just 1 relation/table) 
    // or it's a transitive closure (TC).
    if (path.size() == 1) {
        T = std::stoi(path[0].substr(0, path[0].size()-1));
        std::string relation = path[0].substr(path[0].size()-1, 1);
        
        if (relation == ">") { // (s,t) such that (s, l, t)
            if (q->s == "*") {
                if (q->t =="*") { // - Source: *, Target: *
                    noSources = histogram.distinct_source_relations[T];
                    noPaths = histogram.total_relations[T];
                    noTargets = histogram.distinct_target_relations[T];
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    int result = histogram.target_relations_count[T][t_i];
                    noSources = result;
                    noPaths = result;
                    noTargets = 1;                   
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // - Source: i, Target: *
                    int result = histogram.source_relations_count[T][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // - Source: i, Target: j
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
                noSources = histogram.distinct_target_relations[T];
                noPaths = histogram.total_relations[T];
                noTargets = histogram.distinct_source_relations[T];
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);
                    int result = histogram.source_relations_count[T][t_i];
                    noSources = result; 
                    noPaths = result; 
                    noTargets = 1;                  
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // - Source: i, Target: *
                    int result = histogram.target_relations_count[T][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // - Source: i, Target: j
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
    } else if(path.size() > 1) {

        /// Cases of joins:
        /// Order doesn't matter => s = "*" and t = "*"
        /// Order right to left => s = "*" and t = 1
        /// Order left to right => s = 1 and t = "*", so reverse
        if (q->s != "*") {
            std::reverse(path.begin(), path.end());
        }

        /// Source and target tuple/"Table" 
        int T_s = std::stoi(path[0].substr(0, path[0].size()-1));
        int T_t = std::stoi(path[path.size()-1].substr(0, path[0].size()-1));
        float card = 1;

        bool containsTC;
        for(int = 0; i < path.size(); i++) {
            std::string relation = path[i].substr(path[i].size()-1, 1);
            if(relation == "+") {
                containsTC = true;
                break;
            }
        }

        if (!containsTC) {
            if (q->s == "*") {
                if (q->t =="*") { // - Source: *, Target: *
                /// Basic approach to joins from the characteristic set paper,
                /// Of a star join.
                int T = 0;
                for (int i = 0; i < path.size(); i++) {
                    T = std::stoi(path[i].substr(0, path[i].size()-1));
                    card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
                }
                card = card * graph->getNoVertices();
                
                noSources = histogram.distinct_target_relations[T_s];
                noPaths = card;
                noTargets = histogram.distinct_source_relations[T_t];
                
                } else { // - Source: *, Target: i
                    int t_i = std::stoi(q->t);

                    int T = 0;
                    for (int i = 0; i < path.size(); i++) {
                        T = std::stoi(path[i].substr(0, path[i].size()-1));
                        card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
                    }
                    card = card * graph->getNoVertices();

                    noSources = card; 
                    noPaths = card; 
                    noTargets = 1;                  
                }
            } else {
                int s_i = std::stoi(q->s);
                if (q->t =="*") { // - Source: i, Target: *
                    int T = 0;
                    for (int i = 0; i < path.size(); i++) {
                        T = std::stoi(path[i].substr(0, path[i].size()-1));
                        card = card * (float)histogram.total_relations[T]/(float)graph->getNoVertices();
                    }
                    card = card * graph->getNoVertices() * 1/histogram.total_relations[T_s];

                    noSources = 1; 
                    noPaths = card; 
                    noTargets = card;

                } else { // - Source: i, Target: j
                    /// TODO: Bad estimates
                    int t_i = std::stoi(q->t);
                    int result = std::min(histogram.source_relations_count[T_s][t_i], 
                        histogram.target_relations_count[T_t][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        }
    } else { /// The path contains TC
        /// TODO: Solve TC join
    }

    return cardStat {noSources, noPaths, noTargets};
}
