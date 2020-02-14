#include "SimpleGraph.h"
#include "SimpleEstimator.h"

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // Preperation
    std::cout << "No Edges: " << graph->getNoEdges() << std::endl;
    std::cout << "No Labels: "  << graph->getNoLabels() << std::endl;
    std::cout << "No Vertices: "  << graph->getNoVertices() << std::endl;
    std::cout << "No Distinct Edges: "  << graph->getNoDistinctEdges() << std::endl;

    int noLabels = graph->getNoLabels();
    int noVertices = graph->getNoVertices();

    // Sample tuple sizes
    // TODO: Find good ratio, sample 5% of the vertices
    int samples = ceil(0.05 * noVertices);

    std::vector<int> sampleCount(noLabels, 0);
    for(int i = 0; i < samples; i++) {
        int index = rand() % noVertices;

        if (graph->adj[index].size() > 0) {
            int adjIndex = rand() % graph->adj[index].size();
            int ct = graph->adj[index][adjIndex].first;
            sampleCount[ct] = sampleCount[ct] + 1;
        } else { // No relation edge found,
            i--;
        }
    }
    sampleVertices = sampleCount;


    std::cout << "tuple 0: " << 20 * sampleCount[0]
        << "\ntuple 1: " << 20 * sampleCount[1] 
        << "\ntuple 2: " << 20 * sampleCount[2] 
        << "\ntuple 3: " << 20 * sampleCount[3] 
        << "\nCount: " << 20 * (sampleCount[0] + sampleCount[1] 
        + sampleCount[2] + sampleCount[3]) << std::endl;
}

// TODO: Needs some kind of Histogram to work properly
// V(R, A)
int DistinctValuesFor(int relation, std::string attribute) {
    // hardcoded default is 10
    return 10;
}


cardStat SimpleEstimator::estimate(PathQuery *q) {

    uint32_t costNoOut, costNoPaths, costNoIn = 0;

    // Source defined
    // Selection: 1/V(R, A) * T_R
    std::cout << q->s << std::endl;
    if(q->s != "*") {
        costNoOut += DistinctValuesFor(1, q->s);
    }

    // Target defined
    if(q->t != "*") {
        costNoIn += DistinctValuesFor(1, q->t);
    }
    // Source and Target defined

    // Joins

    // Transitive Closure


    return cardStat {costNoOut, costNoPaths, costNoIn};
}

