#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include <set>
#include <cmath>

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // Preparation
    std::cout << "No Edges: " << graph->getNoEdges() << std::endl;
    std::cout << "No Labels: "  << graph->getNoLabels() << std::endl;
    std::cout << "No Vertices: "  << graph->getNoVertices() << std::endl;
    std::cout << "No Distinct Edges: "  << graph->getNoDistinctEdges() << std::endl;

    int noLabels = graph->getNoLabels();
    int noVertices = graph->getNoVertices();

    // Sample tuple sizes
    // TODO: Find good ratio, sample 5% of the vertices
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
    sampleVertices = sampleCount;


    std::cout << "tuple 0: " << 20 * sampleCount[0]
        << "\ntuple 1: " << 20 * sampleCount[1] 
        << "\ntuple 2: " << 20 * sampleCount[2] 
        << "\ntuple 3: " << 20 * sampleCount[3] 
        << "\nCount: " << 20 * (sampleCount[0] + sampleCount[1] 
        + sampleCount[2] + sampleCount[3]) << std::endl;


    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> relations = {};
    std::vector<std::vector<uint32_t>> relations_size = {};
    for (int i = 0; i < noLabels; i++) {
        relations.push_back({});
        relations_size.push_back({});
        for (int k = 0; k < noVertices; k++) {
            relations_size[i].push_back({0});
        }
    }

    for (uint32_t i = 0; i < graph->adj.size(); i++){
        for (uint32_t j = 0; j < graph->adj[i].size() ; j++) {
            uint32_t rel_type = graph->adj[i][j].first;
            relations_size[rel_type][i] += 1;
            relations[rel_type].push_back(std::make_pair(i, graph->adj[i][j].second));
        }
    }

    for (int i = 0; i < noLabels; i++) {
        uint32_t bla = 0;
        for (int k = 0; k < noVertices; k++) {
            bla += relations_size[i][k];
//            std::cout << relations_size[i][k] << std::endl;
        }
        std::cout << "Relation " << i << ": " << bla << std::endl;
    }

//    for (uint32_t i = 0; i < relations.size(); i++){
//        for (uint32_t j = 0; j < relations[i].size() ; j++){
//            std::cout << relations[i][j].first << "  " << relations[i][j].second << std::endl;
//        }
//        std::cout << std::endl;
//    }
    relation_vector = relations;


//______________________________________________________________________________________
//  Equi depth histogram creation
    uint32_t noBuckets = 10;
    uint32_t depth = relations.max_size()/10;
    std::vector<std::vector<std::vector<uint32_t>>> equidepth_buckets = {};
    for (int i = 0; i < noLabels; i++) {
        uint32_t n = 0;
        equidepth_buckets.push_back({});
        equidepth_buckets[i].push_back({});
        for (int j = 0; j < relations_size[i].size(); j++) {
            uint32_t start = j;
            uint32_t end = j;
            uint32_t sum = 0;
            while (sum < depth) {
                if ((sum > 0) && (relations_size[i][j] > depth)) {
                    break;
                }
                else {
                    j += 1;
                    end = j;
                    sum += relations_size[i][j];
                }
            }
            equidepth_buckets[i][n].push_back(start);
            equidepth_buckets[i][n].push_back(end);
            equidepth_buckets[i][n].push_back(sum);
            n += 1;
            if (j < relations_size[i].size() - 1) {
                equidepth_buckets[i].push_back({});
            }
        }
    }

    for (int i = 0; i < equidepth_buckets[0].size(); i++) {
        std::cout << equidepth_buckets[0][i][0] << "  " << equidepth_buckets[0][i][1] << "  " << equidepth_buckets[0][i][2] << std::endl;
    }

//______________________________________________________________________________________


}

// TODO: Needs some kind of Histogram to work properly
// V(R, A)
int DistinctValuesFor(int relation, std::string attribute) {
    // hardcoded default is 10
    return 10;
}


cardStat SimpleEstimator::estimate(PathQuery *q) {

    // perform your estimation here
    std::set<int> sources = {};
    std::set<int> targets = {};
    std::vector<std::pair<uint32_t,uint32_t>> results = {};
    std::string relation_direction = q->path->data.substr(q->path->data.size()-1,1);
    std::string relation_type = q->path->data.substr(0,q->path->data.size()-1);
    std::string sourceVertex;
    std::string targetVertex;
    if (relation_direction == ">") {
        sourceVertex = q->s;
        targetVertex = q->t;
    }
    else if (relation_direction == "<") {
        sourceVertex = q->t;
        targetVertex = q->s;
    }

    if (sourceVertex != "*") {
        uint32_t int_source = std::stoul(sourceVertex,0);
        if (targetVertex != "*") {
            uint32_t int_target = std::stoul(targetVertex,0);
            for (uint32_t j = 0; j < graph->adj.at(int_source).size() ; j++){
                if ((graph->adj.at(int_source).at(j).first == std::stoul(relation_type,0)) && (graph->adj.at(int_source).at(j).second == int_target)) {
                    results.push_back(std::make_pair(int_source, int_target));
                    sources.insert(int_source);
                    targets.insert(int_target);
                }
            }
        }
        else {
            for (uint32_t j = 0; j < graph->adj.at(int_source).size() ; j++){
                if (graph->adj.at(int_source).at(j).first == std::stoul(relation_type,0)) {
                    results.push_back(std::make_pair(int_source, graph->adj.at(int_source).at(j).second));
                    sources.insert(int_source);
                    targets.insert(graph->adj.at(int_source).at(j).second);
                }
            }
        }
    }
    else {
        if (targetVertex != "*") {
            uint32_t int_target = std::stoul(targetVertex,0);
            for (uint32_t i = 0; i < graph->adj.size(); i++) {
                for (uint32_t j = 0; j < graph->adj.at(i).size() ; j++){
                    if ((graph->adj.at(i).at(j).first == std::stoul(relation_type,0)) && (graph->adj.at(i).at(j).second == int_target)) {
                        results.push_back(std::make_pair(i, int_target));
                        sources.insert(i);
                        targets.insert(int_target);
                    }
                }
            }
        }
        else {
            for (uint32_t i = 0; i < graph->adj.size(); i++) {
                for (uint32_t j = 0; j < graph->adj.at(i).size() ; j++){
                    if (graph->adj.at(i).at(j).first == std::stoul(relation_type,0)) {
                        results.push_back(std::make_pair(i, graph->adj.at(i).at(j).second));
                        sources.insert(i);
                        targets.insert(graph->adj.at(i).at(j).second);
                    }
                }
            }
        }
    }

    uint32_t noSources = sources.size();
    uint32_t noPaths = results.size();
    uint32_t noTargets = targets.size();

    return cardStat {noSources, noPaths, noTargets};
}
