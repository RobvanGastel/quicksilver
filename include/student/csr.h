#ifndef QUICKSILVER_CSR_H
#define QUICKSILVER_CSR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include <set>
// #include "Graph.h"

class csr {// : public Graph {
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

    csr() : V(0), L(0) {};
    ~csr() = default;
    explicit csr(uint32_t n);

    uint32_t getNoVertices() const ;
    uint32_t getNoEdges() const;
    uint32_t getNoDistinctEdges() const;
    uint32_t getNoLabels() const;

    void readFromContiguousFile(const std::string &fileName);
    void initialize_positions_adj();
    void readInitialInfoFromContiguousFile(const std::string &fileName) ;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);
    std::vector<uint32_t> findNeighbours(uint32_t id, uint32_t label, bool reverse);
};

#endif //QUICKSILVER_CSR_H