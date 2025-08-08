#include "Pulse_Analysis.h"
#include "Pulse_Fitting.h"
#include "File_Loader.h"
#include <json.hpp>
#include <iostream>
#include <fstream>
#include <set>


using namespace std;

void analysis_setup(const vector<EventList> run_data, json params, string output_folder) { // Event format: <time (us), PE #, event #, window width, # of events in window>
	double start = (double)params["fill_time"] + (double)params["hold_time"] + (double)params["clean_time"] + 40;
	double stop = start + 60;
	double bg_start = stop + 50;
	
	vector<string> segment_labels = {"12", "34", "56", "78"};
	string output_file = output_folder + "results/PulseAnalysis_" + to_string(params["run_number"]) + ".csv";

	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "Error opening output file: " << output_file << endl;
		return;
	}

	out << "Segment, Time (us), PE, Event\n";

	for (size_t seg = 0; seg < run_data.size(); ++seg) {
		cout << "Segment: " << segment_labels[seg] << endl;
		Pulse_Fitting fitter(run_data[seg]);
		fitter.setWindow(start * 1e6, stop * 1e6);
		fitter.setBackgroundWindow(bg_start * 1e6);
		fitter.analyze();
		auto signalPulses = fitter.getSignalPulses();
		auto backgroundPulses = fitter.getBackgroundPulses();

		for (const auto& event : signalPulses) {
			out << segment_labels[seg] << ", "
				<< get<0>(event)/1e6 << ", "
				<< get<1>(event) << ", "
				<< "1 \n";
		}

		for (const auto& event : backgroundPulses) {
			out << segment_labels[seg] << ", "
				<< get<0>(event)/1e6 << ", "
				<< get<1>(event) << ", "
				<< "0 \n";
		}
	}
	
	out.close();
}
	
int main(int argc, char **argv) {
	Config cfg;
	try {
		cfg = load_config(argc, argv);
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

	std::string data_folder   = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
    bool        save_to_txt   = cfg.save_to_txt;
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;
	vector<EventList> run_data;
	
	if (save_to_txt) {
		cout << "** Note: converting data to text, no analysis will be performed **" << endl;
	}

	for (int z =startrun; z<endrun;z++){

		string run = std::to_string(z);
		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

		if (params.contains(run) && params[run]["run_type"] == "production") {
			if (save_to_txt) {
				processfile(data_folder, output_folder, run);
			} else {
				run_data = processfile(data_folder, run);
				if (run_data.empty()) {
					cerr << "No data found for run " << run << ". Skipping analysis." << endl;
					continue;
				}
				analysis_setup(run_data, params[run], output_folder);
			}
		} else {
			cerr << "Run " << run << " not found or not a production run. Skipping." << endl;
		}

	}	

	return 0;
}