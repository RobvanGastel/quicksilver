#ifndef QUICKSILVER_PATHQUERY_H
#define QUICKSILVER_PATHQUERY_H


#include <cstdint>
#include "PathTree.h"
#include <iostream>

class PathQuery {

public:
    std::string s;
    std::string t;
    PathTree *path;

    PathQuery(std::string &source, PathTree *tree, std::string &target) : s(source), path(tree), t(target) {}
    ~PathQuery();

};

std::ostream & operator << (std::ostream & out, const PathQuery & pq);

#endif //QUICKSILVER_PATHQUERY_H
