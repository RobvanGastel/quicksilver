#ifndef QUICKSILVER_CSR_H
#define QUICKSILVER_CSR_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
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
    std::vector<uint32_t> LabelCount;                // TODO: throw after use
    std::vector<std::vector<uint32_t>> LabelSource;  // TODO: throw after use
    std::vector<std::vector<uint32_t>> LabelTarget;  // TODO: throw after use


public:

    csr() : V(0), L(0) {};
    ~csr() = default;
    explicit csr(uint32_t n);

    uint32_t getNoVertices() const ;
    uint32_t getNoEdges() const;
    uint32_t getNoDistinctEdges() const;
    uint32_t getNoLabels() const;

    // Returns starting index of specified label and ending index+1 (so second is not included in the label)
    std::pair<uint32_t, uint32_t> SelectLabel(uint32_t label, bool reverse);

    void readFromContiguousFile(const std::string &fileName);
    void initialize_positions_adj();
    void readInitialInfoFromContiguousFile(const std::string &fileName) ;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

};

#endif //QUICKSILVER_CSR_H
