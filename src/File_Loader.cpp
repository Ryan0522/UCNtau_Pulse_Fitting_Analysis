#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include <stdexcept>
#include "File_Loader.h"

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
    
    if (tems == NULL || tmcs_0 == NULL || tmcs_1 == NULL) {
		cerr << "Error: Trees not found in file " << filename << endl;
		cerr << "TFile Structure: " << endl;
		fin->ls();
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

Config load_config(int argc, char** argv, const std::string& default_cfg) {
	Config c;

	std::string cfg_path = (argc == 2 ? argv[1] : default_cfg);
	json cfg;
	{
		std::ifstream f(cfg_path);
		if (!f) throw std::runtime_error("Could not open config file: " + cfg_path);
		f >> cfg;
	}

	if (argc >= 8) {
        cfg["data_folder"]   = std::string(argv[1]);
        cfg["output_folder"] = std::string(argv[2]);
        cfg["runinfo"]       = std::string(argv[3]);
        cfg["good_runs"]     = std::string(argv[4]);
        cfg["start_run"]     = std::stoi(argv[5]);
        cfg["end_run"]       = std::stoi(argv[6]);
        cfg["save_to_txt"]   = (string(argv[7]) == "true");
    } else if (argc != 1 && argc != 2) {
        throw std::runtime_error(
            "Expected 7 args in full mode.\n"
            "Usage (config): prog [config.json]\n"
            "Usage (full):   prog ./data/ ./out/ ./runinfo.json ./good_runs.txt <start> <end> <true|false>");
    }

	c.data_folder = cfg.value("data_folder", "../UCNtau_2022_raw_data/");
    c.output_folder = cfg.value("output_folder", "./output/");
    c.runinfo_path = cfg.value("runinfo", "./config/betteroutput2022.json");
    c.good_runs_path = cfg.value("good_runs", "./config/2022runlist.txt");
    c.start_run = cfg.value("start_run", 0);
    c.end_run = cfg.value("end_run", 0);
    c.save_to_txt = cfg.value("save_to_txt", false);

    {
        std::ifstream i(c.runinfo_path);
        if (!i) throw std::runtime_error("Cannot open runinfo JSON: " + c.runinfo_path);
        i >> c.runinfo_json;
    }

	std::set<std::string> good_runs;
	{
		std::ifstream grf(c.good_runs_path);
		std::string line;
		while (std::getline(grf, line)) {
			if (!line.empty()) good_runs.insert(line);
		}
	}
    c.good_runs_set = good_runs;

	return c;
}