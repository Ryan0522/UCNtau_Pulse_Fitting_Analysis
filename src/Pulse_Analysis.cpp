#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include "Pulse_Analysis.h"
#include "Pulse_Fitting.h"

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

    if (argc < 7) {
        cerr << "Error: Expected at least 6 arguemnts.\n";
        cerr << "Usage: " << argv[0] << " ./data_folder/ ./output_folder/ ./runinfo.json\n";
        cerr << " <startRun #> <endRun #> <save_to_txt (true/false)>\n";
        return 1;
    }
   	
    string data_folder = ensureTrailingSlash(argv[1]);
    string output_folder = ensureTrailingSlash(argv[2]);
    string json_filename = argv[3]; // For run info.

    std::ifstream i(json_filename);
    json params;
	i >> params;

   	int startrun = atoi(argv[4]);
    int endrun = atoi(argv[5]);
	bool save_to_txt = (string(argv[6]) == "true");
	string run;
	string hold;
	vector<EventList> run_data;

	if (save_to_txt) {
		cout << "** Note: converting data to text, no analysis will be performed **" << endl;
	}

	for (int z =startrun; z<endrun;z++){

		string run = std::to_string(z);
		if (params.contains(run) && params[run]["run_type"] == "production"){
			int holdtime = params[run]["hold_time"];
			int filltime = params[run]["fill_time"];
			int cleantime = params[run]["clean_time"];
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
		}

	};	

	return 0;
}