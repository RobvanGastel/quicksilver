#ifndef QS_RPQTREE_H
#define QS_RPQTREE_H

#include <string>
#include <algorithm>

class PathTree {

public:
    PathTree *left;
    PathTree *right;
    std::string data;

    PathTree(std::string &payload, PathTree *left, PathTree *right) : data(payload), left(left), right(right) {}
    ~PathTree();

    static PathTree* strToTree(std::string str);

    bool isConcat();
    bool isKleene();

    bool isLeaf();
    bool isUnary();
    bool isBinary();

};

std::ostream & operator << (std::ostream & out, const PathTree & pt);

#endif //QS_RPQTREE_H
