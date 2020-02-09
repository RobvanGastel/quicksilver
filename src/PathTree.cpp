#include <iostream>
#include "PathTree.h"

PathTree::~PathTree() {
    delete(left);
    delete(right);
}

PathTree* PathTree::strToTree(std::string str) {

    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end()); // remove spaces

    int level = 0; // inside parentheses check

    // case /
    // most right '/' (but not inside '()') search and split
    for(int i=(int) str.size()-1;i>=0;--i){
        char c = str[i];
        if(c == ')'){
            ++level;
            continue;
        }
        if(c == '('){
            --level;
            continue;
        }
        if(level>0) continue;
        if(c == '/'){
            std::string left(str.substr(0,i));
            std::string right(str.substr(i+1));
            std::string payload(1, c);
            return new PathTree(payload, strToTree(left), strToTree(right));
        }
    }

    if(str[0]=='('){
        //case ()
        //pull out inside and to strToTree
        for(int i=0;i<str.size();++i){
            if(str[i]=='('){
                ++level;
                continue;
            }
            if(str[i]==')'){
                --level;
                if(level==0){
                    std::string exp(str.substr(1, i-1));
                    return strToTree(exp);
                }
                continue;
            }
        }
    } else
        // case value
        return new PathTree(str, nullptr, nullptr);

    std::cerr << "Error: parsing RPQ failed." << std::endl;
    return nullptr;
}

std::ostream & operator << (std::ostream & out, const PathTree & pq) {

    if(pq.left == nullptr && pq.right == nullptr) {
        out << ' ' << pq.data << ' ';
    } else {
        out << '(' << pq.data << ' ';
        if(pq.left != nullptr) out << *pq.left;
        if(pq.right!= nullptr) out << *pq.right;
        out << ')';
    }

    return out;
}

bool PathTree::isConcat() {
    return (data == "/") && isBinary();
}

bool PathTree::isBinary() {
    return left != nullptr && right != nullptr;
}

bool PathTree::isLeaf() {
    return left == nullptr && right == nullptr;
}