#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include "Pulse_Analysis.h"
#include "Pulse_Fitting.h"

using namespace std;

string ensureTrailingSlash(const string& folder) {
    
    if (!folder.empty() && folder.back() != '/') {
        return folder + "/";
    }

    return folder;

}

vector<EventList> processfile(string data_folder, string runnum) {

	string part1 = "processed_output_";
	string part3 = ".root";
	string filename = data_folder+part1+runnum+part3;
	cout << "Processing file: " << filename << endl;
	
	//check if file exists and if it is empty
	if (FILE *file = fopen(filename.c_str(), "r")) {
		if (file == NULL){
			cout << filename;
			fclose(file);
			return {};
		} else {
			fclose(file);
		}
    } else {
        return {};
    } 
    
    TFile* fin = TFile::Open(filename.c_str());
	
    TTree* tmcs_0 = (TTree*)fin->Get("tmcs_0");
    TTree* tmcs_1 = (TTree*)fin->Get("tmcs_1");
    TTree* tems = (TTree*)fin->Get("tems");
    
    if (tems == NULL) {
		return {};
	}
	
	double run_duration = tems->GetMaximum("time");
	double backup_run_duration = tmcs_0->GetMaximum("realtime");
    if(run_duration < 0 and backup_run_duration >0){
		run_duration = backup_run_duration;
	}
	else if(run_duration < 0){
		run_duration = 2200;
	}
	else if(run_duration > 3000 and backup_run_duration < 3000){
		run_duration = backup_run_duration;
	};
	cout << runnum << " final run duration = " << run_duration << endl;

    EventList PMT12, PMT34, PMT1112, PMT1314;
    event evt0;
    tmcs_0->SetBranchAddress("channel",&evt0.channel);
    tmcs_0->SetBranchAddress("edge",&evt0.edge);
    tmcs_0->SetBranchAddress("tag",&evt0.tag);
    tmcs_0->SetBranchAddress("full",&evt0.full);
    tmcs_0->SetBranchAddress("time",&evt0.time);
    tmcs_0->SetBranchAddress("realtime",&evt0.realtime);
    
    event evt1;
    tmcs_1->SetBranchAddress("channel",&evt1.channel);
    tmcs_1->SetBranchAddress("edge",&evt1.edge);
    tmcs_1->SetBranchAddress("tag",&evt1.tag);
    tmcs_1->SetBranchAddress("full",&evt1.full);
    tmcs_1->SetBranchAddress("time",&evt1.time);
    tmcs_1->SetBranchAddress("realtime",&evt1.realtime);
    
    //pretty sure the PMT pairs are ch1+ch2, ch3+ch4, ch11+ch12, ch13+ch14 (Lara)
    for (long i=0; i<tmcs_0->GetEntries();tmcs_0->GetEntry(i++))
    {
		if (evt0.channel == 1 or evt0.channel == 2) {PMT12.push_back(evt0);}
		else if (evt0.channel == 3 or evt0.channel == 4) {PMT34.push_back(evt0);}

    }

    for (long j=0; j<tmcs_1->GetEntries();tmcs_1->GetEntry(j++))
    {
        if (evt1.channel == 11 or evt1.channel == 12) {PMT1112.push_back(evt1);}
		else if (evt1.channel == 13 or evt1.channel == 14) {PMT1314.push_back(evt1);}
    }

	vector<EventList> result;
	result.push_back(PMT12);
	result.push_back(PMT34);
	result.push_back(PMT1112);
	result.push_back(PMT1314);
	return result;

}

void processfile(string data_folder, string output_folder, string runnum) {

	vector<EventList> result = processfile(data_folder, runnum);
	EventList PMT12, PMT34, PMT1112, PMT1314;

	if (result.size() == 4) {
		PMT12 = result[0];
		PMT34 = result[1];
		PMT1112 = result[2];
		PMT1314 = result[3];
	} else {
		return;
	}

	ofstream output_file;
	string out_part1 = "PECountsRun";
	string out_part3 = ".txt";
	string outfile = output_folder+out_part1+runnum+out_part3;
	cout << "Writing to file: " << outfile << endl;

	output_file.open(outfile, fstream::app);
	auto ti = PMT12.begin();
	while (ti!=PMT12.end()){
		auto PEtime = (*ti).realtime;
		output_file << "12" << ", " << setprecision(15)<< PEtime << "," << (*ti).channel << endl;
		ti++;
	}
	auto ti2 = PMT34.begin();
	while (ti2!=PMT34.end()){
		auto PEtime = (*ti2).realtime;
		output_file << "34" << ", " << setprecision(15)<< PEtime << "," << (*ti2).channel << endl;
		ti2++;
	}
	auto ti3 = PMT1112.begin();
	while (ti3!=PMT1112.end()){
		auto PEtime = (*ti3).realtime;
		output_file << "56" << ", " << setprecision(15)<< PEtime << "," << (*ti3).channel << endl;
		ti3++;
	}
	auto ti4 = PMT1314.begin();
	while (ti4!=PMT1314.end()){
		auto PEtime = (*ti4).realtime;
		output_file << "78" << ", " << setprecision(15)<< PEtime << "," << (*ti4).channel << endl;
		ti4++;
	}
	output_file.close();

    return;
}

void analysis_setup(const vector<EventList> run_data, json params, string output_folder) { // Event format: <time (us), PE #, event #, window width, # of events in window>
	double start = (double)params["fill_time"] + (double)params["hold_time"] + (double)params["clean_time"] + 40;
	double stop = start + 60;
	double bg_start = stop + 50;
	Pulse_Fitting fitter(run_data[0]);
	fitter.setWindow(start * 1e6, stop * 1e6);
	fitter.setBackgroundWindow(bg_start * 1e6);
	fitter.analyze();
	auto signalPulses = fitter.getSignalPulses();
	auto backgroundPulses = fitter.getBackgroundPulses();
	
	string output_file = output_folder + "results/PulseAnalysis_" + to_string(params["run_number"]) + ".csv";
	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "Error opening output file: " << output_file << endl;
		return;
	}

	out << "Time (us), PE, Event\n";

	for (const auto& event : signalPulses) {
		out << get<0>(event)/1e6 << ", "
			<< get<1>(event) << ", "
			<< "1 \n";
	}

	for (const auto& event : backgroundPulses) {
		out << get<0>(event)/1e6 << ", "
			<< get<1>(event) << ", "
			<< "0 \n";
	}

	out.close();
}
	
int main(int argc, char **argv) {

    if (argc < 11) {
        cerr << "Error: Expected at least 9 arguemnts.\n";
        cerr << "Usage: " << argv[0] << " ./data_folder/ ./output_folder/ ./runinfo.json\n";
        cerr << " <startRun #> <endRun #> <threshold>\n";
        cerr << " <coincidence window> <prompt window> <telescoping window> <save_to_txt (true/false)>\n";
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
    int threshold = atoi(argv[6]);
    double wc = atof(argv[7]);
    double wp = atof(argv[8]);
    double wt = atof(argv[9]);
	bool save_to_txt = string(argv[10]) == "true";
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
				analysis_setup(run_data, params[run], output_folder);
			}
		}

	};	

	return 0;
}