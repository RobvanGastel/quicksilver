#include <benches.h>
#include <iostream>
#include <rss.h>

int main(int argc, char *argv[]) {

    if(argc < 3) {
        std::cout << "Usage: quicksilver <graphFile> <queriesFile>" << std::endl;
        return 0;
    }

    // args
    std::string graphFile {argv[1]};
    std::string queriesFile {argv[2]};

    // estimatorBench(graphFile, queriesFile);
    auto result = evaluatorBench(graphFile, queriesFile);
	
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "Total load time: " << result.loadTime << " ms" << std::endl;
    std::cout << "Total prep time: " << result.prepTime << " ms" << std::endl;
    std::cout << "Total eval time: " << result.evalTime << " ms" << std::endl;
    double memoryUsage = double(getPeakRSS()) / 1024.0 / 1024.0;
    std::cout << "Peak memory usage (for all workloads): " << memoryUsage << " MiB" << std::endl;
    
    
    return 0;
}



