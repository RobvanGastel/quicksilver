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

class SimpleGraph {
public:
    std::vector<uint32_t> IA;
    std::vector<std::vector<uint32_t>> positions_adj; // positions_adj[label][source] -> starting index in IA; positions_adj[label][source+1] -> ending index in IA
    std::vector<uint32_t> IA_reverse;
    std::vector<std::vector<uint32_t>> positions_adj_reverse;
protected:
    uint32_t V;
    uint32_t L;
    std::vector<uint32_t> LabelCount;
    std::vector<std::vector<uint32_t>> LabelSource;
    std::vector<std::vector<uint32_t>> LabelTarget;
    std::vector<std::vector<uint32_t>> subjects;
    std::vector<std::vector<uint32_t>> objects;

public:

    SimpleGraph() : V(0), L(0) {};
    ~SimpleGraph() = default;
    explicit SimpleGraph(uint32_t n);

    uint32_t getNoVertices() const  ;
    uint32_t getNoEdges() const  ;
    uint32_t getNoDistinctEdges() const  ;
    uint32_t getNoLabels() const  ;

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);
    void readFromContiguousFile(const std::string &fileName) ;
    void initialize_positions_adj();
    void readInitialInfo(const std::string &fileName);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);
    std::vector<uint32_t> findNeighbours(uint32_t id, uint32_t label, bool reverse);

};

#endif //QS_SIMPLEGRAPH_H
