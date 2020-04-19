#include "AdjMatrix.h"

AdjMatrix::AdjMatrix(uint32_t n)   {
    setNoVertices(n);
}

uint32_t AdjMatrix::getNoLabels() const {
    return L;
}

void AdjMatrix::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    adj_matrix.resize(L);
}

uint32_t AdjMatrix::getNoVertices() const {
    return V;
}

void AdjMatrix::setNoVertices(uint32_t n) {
    V = n;
    for (int i = 0; i < L; i++) {
        adj_matrix[i].resize(V);
        for (int j = 0; j < V; j ++)
            adj_matrix[i][j].resize(V);
    }
}

uint32_t AdjMatrix::getNoEdges() const {
    return E;
}

uint32_t AdjMatrix::getNoDistinctEdges() const {
    return dist_E;
}

void AdjMatrix::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                 "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                 std::to_string(edgeLabel) + ")");
    if ((int)adj_matrix[edgeLabel][from][to] == 0) {
        adj_matrix[edgeLabel][from][to] = 1;
        dist_E += 1;
    }
    E += 1;
}

bool AdjMatrix::edgeExists(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    return (adj_matrix[edgeLabel][from][to] == 1);
}

void AdjMatrix::readFromContiguousFile(const std::string &fileName) {

    std::string line;
    std::ifstream graphFile { fileName };

    std::regex edgePat (R"((\d+)\s(\d+)\s(\d+)\s\.)"); // subject predicate object .
    std::regex headerPat (R"((\d+),(\d+),(\d+))"); // noNodes,noEdges,noLabels

    // parse the header (1st line)
    std::getline(graphFile, line);
    std::smatch matches;
    if(std::regex_search(line, matches, headerPat)) {
        uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
        uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

        setNoLabels(noLabels);
        setNoVertices(noNodes);
    } else {
        throw std::runtime_error(std::string("Invalid graph header!"));
    }
    // parse edge data
    while(std::getline(graphFile, line)) {
        if(std::regex_search(line, matches, edgePat)) {
            uint32_t subject = (uint32_t) std::stoul(matches[1]);
            uint32_t predicate = (uint32_t) std::stoul(matches[2]);
            uint32_t object = (uint32_t) std::stoul(matches[3]);

            addEdge(subject, object, predicate);
        }
    }
//    std::cout << E << "  " << dist_E <<std::endl;
//    printMatrix(0);
    graphFile.close();
}

void AdjMatrix::printMatrix(int label) {
    std::cout << "Relation " << label << std::endl;
    for (int j = 0; j < V; j++) {
        for (int z = 0; z < V; z++) {
            if ((int)adj_matrix[label][j][z] == 1)
                std::cout << (int)adj_matrix[label][j][z] << " ";
        }
        std::cout << std::endl;
    }
}