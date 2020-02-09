#include "PathQuery.h"

PathQuery::~PathQuery() { delete(path); }

std::ostream & operator << (std::ostream &out, const PathQuery &pq) {
    out << pq.s << ", " << *pq.path << ", " << pq.t << std::endl;
    return out;
}

