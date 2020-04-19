#ifndef QUICKSILVER_K2TREE_H
#define QUICKSILVER_K2TREE_H


#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"

class k2Tree : public Graph {
private:
    std::vector<std::vector<std::vector<uint32_t>>> adj;
    std::vector<std::vector<uint8_t>> trees;
    std::vector<std::vector<uint8_t>> leaves;

protected:
    uint32_t V;
    uint32_t L;
    uint8_t k;
    uint8_t h;

public:

    k2Tree() : V(0), L(0), k(2), h(1) {};
    ~k2Tree() = default;
    explicit k2Tree(uint32_t n);

    uint32_t getNoVertices() const;
    uint32_t getNoEdges() const;
    uint32_t getNoDistinctEdges() const;
    uint32_t getNoLabels() const;
    uint8_t getK();
    uint8_t getH();

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel);
    void readFromContiguousFile(const std::string &fileName);
    void createTrees();

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);
    void setK(uint8_t k_);
    void setH(uint8_t h_);

    uint8_t build (std::vector<std::vector<uint32_t>> adj, std::vector<typename std::vector<uint32_t>::iterator> &cursors, std::vector<std::vector<uint8_t>> &levels, uint32_t n, uint32_t l, uint32_t p, uint32_t q);

};

#endif //QUICKSILVER_K2TREE_H
