#ifndef QUICKSILVER_CSR_H
#define QUICKSILVER_CSR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"

class csr : public Graph {
public:
    std::vector<uint32_t> IA;
    std::vector<std::vector<uint32_t>> positions_adj; // positions_adj[label][source] -> starting index in IA; positions_adj[label][source+1] -> ending index in IA


protected:
    uint32_t V;
    uint32_t L;
    std::vector<int> LabelCount;
    std::vector<std::vector<int>> LabelSource;


public:

    csr() : V(0), L(0) {};
    ~csr() = default;
    explicit csr(uint32_t n);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override ;
    uint32_t getNoLabels() const override ;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override ;
    void readFromContiguousFile(const std::string &fileName) override ;
    void initialize_positions_adj();
    void readInitialInfoFromContiguousFile(const std::string &fileName) ;

    bool edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel);

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

};

#endif //QUICKSILVER_CSR_H
