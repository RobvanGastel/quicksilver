#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include <set>

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // do your prep here

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

//    Print adjacency list
//    for (uint32_t i = 0; i < graph->adj.size(); i++){
//        for (uint32_t j = 0; j < graph->adj.at(i).size() ; j++){
//            std::cout << i << " " << graph->adj.at(i).at(j).first << " " << graph->adj.at(i).at(j).second << " || ";
//        }
//        std::cout << std::endl;
//    }

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