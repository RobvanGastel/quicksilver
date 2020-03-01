#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"
#include <set>
#include <cmath>
#include <cfloat>

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

std::vector<std::string> reverse_path(std::vector<std::string> path) {
    std::vector<std::string> newPath;
    for (int i = path.size()-1; i >= 0 ; i--) {
        if (path[i].substr(1, 2) == ">")
            newPath.push_back(path[i].substr(0, 1)+(std::string)"<");
        else
            newPath.push_back(path[i].substr(0, 1)+(std::string)">");
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

    /// Defaults to 1, unless we know there are no tuples.
    uint32_t noSources = 1;
    uint32_t noPaths = 1;
    uint32_t noTargets = 1;

    /// Either there are no joins (e.g. just 1 relation/table) 
    /// or it's a transitive closure (TC).
    if (path.size() == 1) {
        rel_type = std::stoi(path[0].substr(0, path[0].size()-1));
        std::string relation = path[0].substr(path[0].size()-1, 1);
        
        if (relation == ">") { // forward relation, (s,t) such that (s, l, t)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                    noSources = stats.distinct_source_relations[rel_type];
                    noPaths = stats.total_relations[rel_type];
                    noTargets = stats.distinct_target_relations[rel_type];
                } else { // source: *, target: i
                    int t_i = std::stoi(q->t);
                    int result = stats.target_relations_count[rel_type][t_i];
                    noSources = result;
                    noPaths = result;
                    noTargets = 1;                   
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = stats.source_relations_count[rel_type][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(stats.target_relations_count[rel_type][t_i], 
                        stats.source_relations_count[rel_type][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        } else if(relation == "<") { // backward relation, (s,t) such that (t, l, s)
            if (q->s == "*") {
                if (q->t =="*") { // source: *, target: *
                    noSources = stats.distinct_target_relations[rel_type];
                    noPaths = stats.total_relations[rel_type];
                    noTargets = stats.distinct_source_relations[rel_type];
                } else { // source: *, target: j
                    int t_i = std::stoi(q->t);
                    int result = stats.source_relations_count[rel_type][t_i];
                    noSources = result; 
                    noPaths = result; 
                    noTargets = 1;                  
                }
            } else {
                int s_i = std::stoi(q->s);

                if (q->t =="*") { // source: i, target: *
                    int result = stats.target_relations_count[rel_type][s_i];
                    noSources = 1;
                    noPaths = result;
                    noTargets = result;
                } else { // source: i, target: j
                    int t_i = std::stoi(q->t);
                    int result = std::min(stats.source_relations_count[rel_type][t_i], 
                        stats.target_relations_count[rel_type][s_i]);
                    noSources = result;
                    noPaths = result;
                    noTargets = result;
                }
            }
        }
        else if(relation == "+") { // Transitive Closure relation

            // The sample size for sampling on the entire graph
            float sample = 0.05;

            if (q->s == "*") { // - Source: *
                if (q->t =="*") { // - Source: *, Target: *
                    noSources = stats.distinct_source_relations[rel_type];
                    noTargets = stats.distinct_target_relations[rel_type];
                    noPaths = stats.total_relations[rel_type] + stats.multidimensional_matrix[rel_type][rel_type][0][0][0];
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
                    /// TODO: Implement or ignore?
                    int t_j = std::stoi(q->t);
                }
            }  
        }
    } else if(path.size() > 1) {
        /// There is atleast 1 join in the query,
        /// Cases of joins:
        /// Order doesn't matter => s = "*" and t = "*"
        /// Order right to left => s = "*" and t = 1
        /// Order left to right => s = 1 and t = "*", so reverse
        if (q->t != "*") {
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

                relation_i = get_relation_info(path[0]);
                relation_j = get_relation_info(path[1]);
                // std::cout << "\n    path_0: " << path[0] << "  relation_0: " << relation_i[0] << " " << relation_i[1] << std::endl;

                join_stats = stats.multidimensional_matrix[relation_i[0]][relation_j[0]][relation_i[1]][relation_j[1]];
                T_i = join_stats[0];
                d_si = join_stats[1];
                // int middle_i = join_stats[2];
                d_oi = join_stats[3];
                // std::cout << "        T_i: " << T_i << "  d_si: " << d_si << "  middle_i: " << middle_i << "  d_oi: " << d_oi << std::endl;
                
                for (int j = 2; j < path.size(); j++) {
                    // std::cout << "\n        results d_si: " << d_si;
                    // std::cout << "        results T_i: " << T_i;
                    // std::cout << "        results d_oi: " << d_oi << std::endl;
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
                    // std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;

                    // calculations
                    part1 = middle_j / in;
                    // part1 = middle_j / d_oi;
                    // std::cout << "        results part1: " << part1 << "  " << in << "  " << d_oi;
                    d_si = d_si * part1;
                    // std::cout << "  " << middle_j << "  " << in << std::endl;
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
            } 
            // else { // source: *, target: j
            //     // should never happen -> reverse path
            //     std::cout << "WRONG: source: *, target: j" << std::endl;
            // }
        } else {
            if (q->t == "*" || q->s == "*") { // source: i, target: *
                // std::cout << "source: i, target: *" << std::endl;
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
                
                d_oi = T_i;
                // std::cout << "    source:" << source << " T_i:" << T_i << std::endl;
                
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
                    // std::cout << "        T_j: " << T_j << "  d_sj: " << d_sj << "  middle_j: " << middle_j << "  d_oj: " << d_oj << std::endl;

                    // calculations
                    part1 = middle_j / in;
                    d_si = d_si * part1;
                    // std::cout << "  " << middle_j << "  " << in << std::endl;
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
                if (q->t == "*") {
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
