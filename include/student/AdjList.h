#ifndef QUICKSILVER_ADJLIST_H
#define QUICKSILVER_ADJLIST_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"

class AdjList : public Graph {
private:
    std::vector<std::vector<std::vector<uint32_t>>> adj;
protected:
    uint32_t V;
    uint32_t L;

public:

    AdjList() : V(0), L(0) {};
    ~AdjList() = default;
    explicit AdjList(uint32_t n);

    uint32_t getNoVertices() const;
    uint32_t getNoEdges() const;
    uint32_t getNoDistinctEdges() const;
    uint32_t getNoLabels() const;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel);
    void readFromContiguousFile(const std::string &fileName);

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

};

#endif //QUICKSILVER_ADJLIST_H
