#include "Results_Comp.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

vector<tuple<string, double, int>> load_pe_counts_csv(const std::string& filepath) {
    vector<tuple<string, double, int>> data;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filepath << std::endl;
        return data;
    }
    string line;
    while (std::getline(file, line)) {
        if () continue;
        std::istringstream ss(line);
        string segment;
        double time;
        int pe;
        char comma;
        if (ss >> segment >> comma >> time >> comma >> pe) {
            segment.erase(0, segment.find_first_not_of(" \t"));
            segment.erase(segment.find_last_not_of(" \t") + 1);
            data.emplace_back(segment, time, pe);
        }
    }
    return data;
}

int main(int argc, char **argv) {

	Config cfg;
	try {
		cfg = load_config(argc, argv); // parse CLI/config, load runinfo + good runs
		std::cout << "====================================" << std::endl;
		std::cout << "Data folder: "   << cfg.data_folder   << "\n";
        std::cout << "Output folder: " << cfg.output_folder << "\n";
		std::cout << "Runinfo path: "  << cfg.runinfo_path  << "\n";
		std::cout << "Good runs path: "<< cfg.good_runs_path<< "\n";
        std::cout << "Start run: "     << cfg.start_run     << "\n";
        std::cout << "End run: "       << cfg.end_run       << "\n";
        std::cout << "Save to txt: "   << (cfg.save_to_txt ? "true" : "false") << "\n";
        std::cout << "Good runs loaded: " << cfg.good_runs_set.size() << " entries\n";
		std::cout << "====================================" << std::endl;
	} catch (const std::exception& e) {
		cerr << "Error starting program: " << e.what() << endl;
		return 1;
	}

	init_global_pdf(cfg);

	std::string data_folder   = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;
	vector<EventList> run_data;
	
	for (int z =startrun; z<endrun;z++){

		string run = std::to_string(z);
		// skip runs not in good runs list
		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

		if (params.contains(run) && params[run]["run_type"] == "production") {
			
		} else {
			cerr << "Run " << run << " not found or not a production run. Skipping." << endl;
		}

	}	

	return 0;
}