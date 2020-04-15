#ifndef QUICKSILVER_ADJMATRIX_H
#define QUICKSILVER_ADJMATRIX_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"

class AdjMatrix {
public:
    std::vector<std::vector<std::vector<uint8_t>>> adj_matrix;
protected:
    uint32_t V;
    uint32_t L;
    uint32_t E;
    uint32_t dist_E;

public:

    AdjMatrix() : V(0), L(0), E(0), dist_E(0) {};
    ~AdjMatrix() = default;
    explicit AdjMatrix(uint32_t n);

    uint32_t getNoLabels() const;
    uint32_t getNoVertices() const;
    uint32_t getNoEdges() const;
    uint32_t getNoDistinctEdges() const;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel);
    void readFromContiguousFile(const std::string &fileName);

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);
    void printMatrix(int label);

};


#endif //QUICKSILVER_ADJMATRIX_H
