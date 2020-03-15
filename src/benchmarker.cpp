#include <benches.h>
#include <iostream>
#include <rss.h>
#include <vector>
#include <filesystem>

void findBenchmarks(std::vector<std::pair<std::string, std::string>>& benchmarks, const std::string& directory) {
	namespace fs = std::filesystem;
	std::string graphPath, queriesPath;
	
	for (const auto & entry : fs::directory_iterator(directory)) {
		if (entry.status().type() == fs::file_type::directory) {
			findBenchmarks(benchmarks, entry.path());
		}
		if (entry.status().type() == fs::file_type::regular) {
			if (entry.path().filename() == "queries.csv") {
				queriesPath = entry.path();
			}
			if (entry.path().filename() == "graph.nt") {
				graphPath = entry.path();
			}
		}
	}
	
	if (!graphPath.empty() && !queriesPath.empty()) {
		benchmarks.emplace_back(graphPath, queriesPath);
	}
}

int main(int argc, char *argv[]) {

    if(argc < 2) {
        std::cout << "Usage: benchmarker <workloaddir>" << std::endl;
        return 0;
    }

    // args
    std::string workloadDir {argv[1]};
    
    std::vector<std::pair<std::string, std::string>> benchmarks;
    
    findBenchmarks(benchmarks, workloadDir);
	
    struct benchresult_t result = {};
    for (auto& benchmark : benchmarks) {
    	std::cout << "\n=== Benchmark " << benchmark.first << " + " << benchmark.second << "===" << std::endl;
    	
    	auto result1 = evaluatorBench(benchmark.first, benchmark.second);
    	
    	std::cout << std::endl << std::endl;
    	result.loadTime += result1.loadTime;
    	result.evalTime += result1.evalTime;
    	result.prepTime += result1.prepTime;
    	
    	std::cout << "Total load time (for this benchmark): " << result1.loadTime << " ms" << std::endl;
    	std::cout << "Total prep time (for this benchmark): " << result1.prepTime << " ms" << std::endl;
    	std::cout << "Total eval time (for this benchmark): " << result1.evalTime << " ms" << std::endl;
    }
    

    std::cout << std::endl << std::endl << std::endl;
    std::cout << "Total load time: " << result.loadTime << " ms" << std::endl;
    std::cout << "Total prep time: " << result.prepTime << " ms" << std::endl;
    std::cout << "Total eval time: " << result.evalTime << " ms" << std::endl;
    double memoryUsage = double(getPeakRSS()) / 1024.0 / 1024.0;
    std::cout << "Peak memory usage (for all workloads): " << memoryUsage << " MiB" << std::endl;
    
    return 0;
}



