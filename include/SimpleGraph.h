#ifndef QS_SIMPLEGRAPH_H
#define QS_SIMPLEGRAPH_H

#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"
#include <map>

class SimpleGraph : public Graph {
public:
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> adj;
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> reverse_adj; // vertex adjacency list
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> sample_adj;
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> sample_reverse_adj;
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> comp_adj;
protected:
    uint32_t V;
    uint32_t L;
    std::vector<uint32_t> act_sources;
    std::vector<uint32_t> act_targets;

public:

    SimpleGraph() : V(0), L(0) {};
    ~SimpleGraph() = default;
    explicit SimpleGraph(uint32_t n);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override ;
    uint32_t getNoLabels() const override ;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override ;
    void readFromContiguousFile(const std::string &fileName) override ;

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

};

#endif //QS_SIMPLEGRAPH_H
