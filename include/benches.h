#ifndef QS_BENCHES_H
#define QS_BENCHES_H
#include <string>

struct benchresult_t {
	long prepTime, evalTime, loadTime;
};

struct benchresult_t evaluatorBench(std::string &graphFile, std::string &queriesFile);
int estimatorBench(std::string &graphFile, std::string &queriesFile);

#endif //QS_BENCHES_H
