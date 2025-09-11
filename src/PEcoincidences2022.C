//to run: make 2022analysisv2, ./2022analysisv2
#include <iostream>
#include <fstream>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <json.hpp>
#include "TCanvas.h"
#include "TLatex.h"
#include <cmath>
#include <TF1.h>
#include <TH1I.h>
#include <numeric>
#include <set>

using json = nlohmann::json;
using namespace std;

typedef struct
{
    Int_t channel;
    Int_t edge;
    Int_t tag;
    Int_t full;
    ULong64_t time;
    Double_t realtime;
} event;


typedef struct
{
    ULong64_t time;
    Double_t realtime;
    Double_t length;
    Int_t edge;
    Int_t tag;
    Int_t full;
    Int_t N;
    Double_t kappa;
    Double_t tau;
    Double_t Npileup;
} coincidence_event;

typedef struct
{
	Double_t unload;
	Double_t background;
	Double_t pilupunload;
	Double_t pileupbackground;
	Double_t deadtimeunload;
	Double_t deadtimebackground;
	Double_t rdeunload;
	Double_t rdebackground;
	Double_t fillUCN;
} unloads_backgrounds;

// some definitions
using EventList = std::list<event>;
using CoincidenceList = std::list<coincidence_event>;
//using UnloadsBackgroundsList = std::list<unloads_backgrounds>;

// a function to take two lists of times, and find coincidences

CoincidenceList find_coincidences(EventList& l1, double coincidence_window, double prompt_window, double tele_window, double Npe, double* ptr_avgku, double* ptr_gammau)
{
    CoincidenceList coincidences;
    auto ti = l1.begin();
    auto testi = ti++;
    //cout << (*ti).realtime << " test " << (*testi).realtime << endl;
    coincidence_event e = {0,-1.0,0.0,0,0,0,0};
	vector<double> karr;
	vector<double> scalarr;
	vector<double> tauarr;
	vector<double> numphotons;
	vector<double> starttimes;
	vector<double> endtimes;
	//TH1F* photons = new TH1F("photons", "photons: counts: time [s]", 100, 0.0, 1e-5);
	//photons->Fill((*ti).realtime-(*ti).realtime);
	double rate = 0;
	//auto tnext = ti;
	//tnext++;
	//double rate = (*tnext).realtime-(*ti).realtime;
	double prevfreephoton = 0;
	double freephotoncount = 0;
	double Npileup;
	srand( time(NULL) );
    while (ti!=l1.end())
    {
	   
	   int nt = 1; //number of photons found within telescoping window
       //int np = 1; //number of photons found within prompt window
       int state = 0; //state 0 = no coincidence, state 1 = coincidence found
       int coincevt = 0; //if coincevt = 1, there's enough PEs for an event
       auto tj = ti;
       tj++; //tj starts as ti+1
	//photons->Fill((*tj).realtime-(*ti).realtime);
       while (tj!=l1.end()) {
       	//if ((*ti).realtime >= 420.08457 and (*ti).realtime <= 420.08459){
       		//cout << setprecision(20)<<(*tj).realtime << ", " <<(*ti).realtime<< endl;
       		//}
		   if (state == 0){
			   if ((*tj).realtime-(*ti).realtime > coincidence_window) { //break if no PEs found within threshold
				   break;
			   }
			   else if ((*tj).realtime-(*ti).realtime < coincidence_window) { 
				   nt = nt + 1;
				   //np = np + 1;
				   //photons->Fill((*tj).realtime-(*ti).realtime);
				   if ((*tj).channel != (*ti).channel){ //state = 1 only if a PE is found in the other PMT within threshold
					   state = 1;
					   //double ku = (*tj).realtime-(*ti).realtime;
					   //karr.push_back(ku);
					   
				   }
			   }
		   }
		   else if (state == 1){
			   auto tj2 = tj;
			   --tj2;
			   if ((*tj).realtime-(*tj2).realtime >= tele_window) { //continue only if the time between current and last event is less than telescoping window
				   if (nt >= Npe){ //normally np
					   coincevt = 1; 
					   break;
				   }
				   else {
					   break;
				   }
			   }
			   else {
				   //if ((*tj).realtime-(*ti).realtime < prompt_window) {
				//	   np = np + 1;
				   //}
				   nt = nt+1;
				   //photons->Fill((*tj).realtime-(*ti).realtime);
			   }
		   }
		   tj++;
	   }
	   if (coincevt == 1){ //if nt reaches threshold, add coincidence event and restart algorithm at current tj
		   auto tj2 = tj;
			--tj2;
			double length = (*tj2).realtime + tele_window - (*ti).realtime;
			karr.push_back(length);
			//TF1 *f1 = new TF1("f1","[0]*(exp([1]*(x)))", 0 , 1e-5);
			//f1->SetParLimits(0,0, 1e10);
			//photons->Fit("f1", "RQ"); //R is use specified range, Q is quiet/minimal printing
			//TF1 *fitresult = photons->GetFunction("f1");
			//double chi2 = fitresult->GetChisquare();
			//double ndf = fitresult->GetNDF();
			double p1 = 1;//f1->GetParameter(1);
			double p0 = 1;//f1->GetParameter(0);
			tauarr.push_back(p1);
			scalarr.push_back(p0);
			numphotons.push_back(nt);
			starttimes.push_back((*ti).realtime);
			endtimes.push_back((*tj2).realtime+tele_window);
			double randNum = (float)rand()/(float)RAND_MAX;
			Npileup = nt - length/rate + randNum;
			e = {(*ti).time, (*ti).realtime, length, (*ti).edge, (*ti).tag, (*ti).full, nt, p0, p1, Npileup};
			coincidences.push_back(e);
				
			ti = tj;
		}
		else if (coincevt == 0) { //if nt doesn't reach threshold, restart algorithm at ti+1
			if (freephotoncount==0){
				rate = (*ti).realtime-prevfreephoton;
			}
			else if (freephotoncount <100){
				rate = ((freephotoncount-1)/freephotoncount)*rate+(1/freephotoncount)*((*ti).realtime-prevfreephoton);
			}
			else {
				rate = (0.99)*rate+(0.01)*((*ti).realtime-prevfreephoton);
			}
			prevfreephoton = (*ti).realtime;
			freephotoncount = freephotoncount+1;
			ti++;
		}
		/*
		TCanvas *cp = new TCanvas("cp");
		cp->cd();
		photons->Draw();
		cp->SaveAs("photons.root");*/

		//photons->Reset();	   
				   
    }
    
    //histogram that shows coincidence lengths
    
    /*
    TH1D* coinclength = new TH1D("coinclength", "chain lengths; time [s]; counts",1600,0.0, 16e-6);
    for(int i = 0; i<=karr.size(); i++){	
		coinclength->Fill(karr[i]);
	}
	*/
	//time between starts of coincidence events
	
	/*
	TH1D* tbe = new TH1D("tbe", "Time between events; time [s]; counts",1000,0.0, 1e-4);
    for(int i = 1; i<=starttimes.size(); i++){	
		double timebtwn = starttimes[i]-endtimes[i-1];
		tbe->Fill(timebtwn);
	}
	*/
	
    float average = accumulate( karr.begin(), karr.end(), 0.0)/karr.size(); 
    float averagescale = accumulate( scalarr.begin(), scalarr.end(), 0.0)/scalarr.size();
    float averagetau = accumulate( tauarr.begin(), tauarr.end(), 0.0)/tauarr.size();       
    float avgphotons = accumulate( numphotons.begin(), numphotons.end(), 0.0)/numphotons.size();
    *ptr_avgku = average;
    *ptr_gammau = avgphotons;        
    //cout << "avg coinc length: " << average << ", avg p0: " << averagescale << ", avg p1: " << averagetau << ", avg photons: " <<avgphotons << endl;   
    return coincidences;
}

void getunloadsbackgrounds(CoincidenceList cl, string runnum, double avgkuval, double gammauval, double run_duration, int holdtime, int filltime, int cleantime, int threshold, TH1I* ch1ch2, int segment){
	
	double start = filltime + holdtime + cleantime + 40;
	double stop  = start + 60; // end of event window
	double bg_start = stop + 50; // start of bg
	double bg_stop = bg_start + 60; // end of bg
	
	TH1D* hcoinc = new TH1D("hcoinc", "coincidences; time [s]; counts", int(run_duration)*10,0.0, run_duration);
	//TFile* fout = TFile::Open("testfile.root","recreate");
    TTree* tcoinc = new TTree("tcoinc","tree for coincidences");
   
    coincidence_event* ce;
    ULong64_t time;
    Double_t realtime, length, kappa, tau;
    Int_t edge, tag, full, N, Npileup;
    tcoinc->Branch("time",&time);
    tcoinc->Branch("realtime",&realtime);
    tcoinc->Branch("length",&length);
    tcoinc->Branch("edge",&edge);
    tcoinc->Branch("tag",&tag);
    tcoinc->Branch("full",&full);
    tcoinc->Branch("N",&N);
    tcoinc->Branch("kappa", &kappa);
    tcoinc->Branch("tau", &tau);
    tcoinc->Branch("Npileup", &Npileup);
   
    for (auto& e : cl)
    {
        time = e.time;
        realtime = e.realtime;
        length = e.length;
        edge = e.edge;
        tag = e.tag;
        full = e.full;
        N = e.N;
        kappa = e.kappa;
        tau = e.tau;
        Npileup = e.Npileup;
        tcoinc->Fill();
        //if (N >= 8){		//tried to make it so it needed N = 8 counts before recording an event
	//		tcoinc->Fill();}
    }
    
    ofstream thisfile;
    string outputfile = "./output/coincidences/CoincRun" + runnum + "_" + to_string(threshold) + "PEthreshold.csv";
    thisfile.open(outputfile, fstream::app);
    tcoinc->GetEntry(0);
    //cout << realtime << "????" << endl;
	
	auto event_type = [start, stop, bg_start, bg_stop](double t) -> int {
		if (t >= start && t <= stop) return 1;
		else if (t >= bg_start && t <= bg_stop) return 0;
		return -1;
	};

    for (int s=1; s<tcoinc->GetEntries(); tcoinc->GetEntry(s++)){
		hcoinc->Fill(realtime);
		int et = event_type(realtime);
		if (et == -1) continue;
		//cout << segment << ", "<<realtime << ", " << N << endl;
		thisfile << segment << ", "<<realtime << ", " << N << ", " << et << endl;
	}
	thisfile.close();
	delete hcoinc;


}


//TH1I* processfile(string runnum, int holdtime, int filltime, int cleantime, json params){//, TCanvas* pdf_canvas, TCanvas* pdf_canvas2){
double processfile(json params, string runnum, int holdtime, int filltime, int cleantime, int threshold, double wc, double wp, double wt, int firstrun, int lastrun){//, TCanvas* pdf_canvas, TCanvas* pdf_canvas2){

	double start = filltime + holdtime + cleantime + 40;
	double stop  = start + 60; // end of event window
	double bg_start = stop + 50; // start of bg
	double bg_stop = bg_start + 60; // end of bg
	cout << "Start, End time for event coincidences: " << start << ", " << stop << endl;
	cout << "Start, End time for bg coincidences: " << bg_start << ", " << bg_stop << endl;

	string part1 = "../UCNtau_2022_raw_data/processed_output_";
	string part3 = ".root";
	string filename = part1+runnum+part3;
	
	//check if file exists and if it is empty
	if (FILE *file = fopen(filename.c_str(), "r")) {
		if (file == NULL){
			cout << filename;
			fclose(file);
			return 0;
		} else {
			fclose(file);
		}
    } else {
        return 0;
    } 
    
    TFile* fin = TFile::Open(filename.c_str());
	
	//TFile* fin = TFile::Open(argv[2]);
    TTree* tmcs_0 = (TTree*)fin->Get("tmcs_0");
    TTree* tmcs_1 = (TTree*)fin->Get("tmcs_1");
    TTree* tmcs_2 = (TTree*)fin->Get("tmcs_2");
    TTree* tems = (TTree*)fin->Get("tems");
    
    //had issue with files "probably not closed", this seems to help?
    if (tems == NULL) {
		return 0;
	}
	
	double run_duration = tems->GetMaximum("time");
	double backup_run_duration = tmcs_0->GetMaximum("realtime");
    //cout << "run time = "<<run_duration << endl;
    //cout << "backup run time = " << backup_run_duration << endl;
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

    EventList PMT_A, PMT_B, PMT12, PMT34, PMT1112, PMT1314;
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
    
    event evt2;
    tmcs_2->SetBranchAddress("channel",&evt2.channel);
    tmcs_2->SetBranchAddress("edge",&evt2.edge);
    tmcs_2->SetBranchAddress("tag",&evt2.tag);
    tmcs_2->SetBranchAddress("full",&evt2.full);
    tmcs_2->SetBranchAddress("time",&evt2.time);
    tmcs_2->SetBranchAddress("realtime",&evt2.realtime);
    
    //cout << "test1" << endl;
    Int_t tag_previous = 0;
    vector<double> arr;
    vector<int> arr2;
    
    for (long j=0;j<tmcs_1->GetEntries();tmcs_1->GetEntry(j++))
    {
		if ((evt1.tag >> 0 & 1) != (tag_previous >> 0 & 1)) {
			double time1 = evt1.realtime;
			int tag1 = evt1.tag >> 0 & 1;
			//cout << time1 << endl;
			for (double t = j; t<tmcs_1->GetEntries(); tmcs_1->GetEntry(t++)){
				//int breakoutflag = 0;
				if (time1+0.2 > evt1.realtime){
					if ((evt1.tag >> 0 & 1) != tag1){
						//breakoutflag = 1;
						break;
					}
				}
				else {
					//j = t;
					arr.push_back(time1);
					arr2.push_back(tag1);
					//cout << "dagger: " <<evt1.realtime<< ", " << time1 << ", " << tag1 << endl;
					break;
				}
			}
									
		}
		tag_previous = evt1.tag;
	}
	int flag = 0;
	float daggerup = 0;
	float daggerdown = 0;
	for (int i=0;i<arr.size();i++){
		if (arr2[i]==1 && flag ==1){
			cout << "holding time: " << arr[i]-300-50  << endl;
			daggerup = arr[i];
			flag = 0;
		}
		if (arr2[i]==0 && flag == 1){
			daggerdown = arr[i];
			flag = 1;
		}
		if (int(arr[i])==300 && arr2[i] ==1){
			flag = 1;
		}
	}
	
	vector<double> arr3;
    vector<int> arr4;
	tag_previous = 0;
	for (long j=0;j<tmcs_1->GetEntries();tmcs_1->GetEntry(j++))
    {
		if ((evt1.tag >> 2 & 1) != (tag_previous >> 2 & 1)) {
			double time1 = evt1.realtime;
			int tag1 = evt1.tag >> 2 & 1;
			for (double t = j; t<tmcs_1->GetEntries(); tmcs_1->GetEntry(t++)){
				if (time1+0.2 > evt1.realtime){
					if ((evt1.tag >> 2 & 1) != tag1){
						break;
					}
				}
				else {
					arr3.push_back(time1);
					arr4.push_back(tag1);
					//cout << "giant cleaner: " << evt1.realtime<< ", " << time1  << ", " << tag1 << endl;
					break;
				}
			}
									
		}
		tag_previous = evt1.tag;
	}
	flag = 0;
	float cleanerup = 0;
	for (int i=0;i<arr3.size();i++){
		if (arr4[i]==0 && flag ==1){
			cleanerup = arr3[i];
			flag = 0;
		}
		if (int(arr3[i])==0 && arr4[i] ==1){
			flag = 1;
		}
	}
	
	float tagholdtime = daggerup-cleanerup;
    
	float m = 0.1;
    
    //create histograms
    string gvtitle = runnum + " GV counts;time [s];counts";
    TH1I* hGV = new TH1I("hGV",                       // name of histogram to be stored in output file
                         gvtitle.c_str(), // title ; x label ; y label
                         int(400*m),        // number of bins
                         0.0,                         // lowest bin value
                         400);               // highest bin value
    
    TH1I* ch4_2 = new TH1I("ch4_2", "ch4_2 counts;time [s];counts",int(400*m),0.0,400);
    TH1I* ch4 = new TH1I("ch4", "ch4 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch3 = new TH1I("ch3", "ch3 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* hGV_2 = new TH1I("hGV_2", "hGV_2 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch5 = new TH1I("ch5", "ch5 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch6 = new TH1I("ch6", "ch6 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch1 = new TH1I("ch1", "ch1 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch2 = new TH1I("ch2", "ch2 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch11 = new TH1I("ch11", "ch11 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch12 = new TH1I("ch12", "ch12 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch13 = new TH1I("ch13", "ch13 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch14 = new TH1I("ch14", "ch14 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch15 = new TH1I("ch15", "ch15 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch16 = new TH1I("ch16", "ch16 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch21 = new TH1I("ch21", "ch21 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch22 = new TH1I("ch22", "ch22 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch23 = new TH1I("ch23", "ch23 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    TH1I* ch24 = new TH1I("ch24", "ch24 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    
    TH1I* ch1ch2 = new TH1I("ch1ch2", "ch1 and ch2 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
	TH1I* ch3ch4 = new TH1I("ch3ch4", "ch3 and ch4 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
	TH1I* ch11ch12 = new TH1I("ch11ch12", "ch11 and ch12 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
	TH1I* ch13ch14 = new TH1I("ch13ch14", "ch13 and ch14 counts;time [s];counts",int(run_duration)*10,0.0,run_duration);
    
    //pretty sure the PMT pairs are ch1+ch2, ch3+ch4, ch11+ch12, ch13+ch14
	//RHDD is probably ch5, GV is probably ch15
    //fill histograms
    for (long i=0; i<tmcs_0->GetEntries();tmcs_0->GetEntry(i++))
    {
		if (evt0.realtime < start || evt0.realtime > bg_stop) continue;
		if (evt0.channel == 1 or evt0.channel == 2) {PMT12.push_back(evt0); ch1ch2->Fill(evt0.realtime);}
		else if (evt0.channel == 3 or evt0.channel == 4) {PMT34.push_back(evt0); ch3ch4->Fill(evt0.realtime);}
        if      (evt0.channel == 1) {PMT_A.push_back(evt0); ch1->Fill(evt0.realtime);} //maybe?
        else if (evt0.channel == 2) {PMT_B.push_back(evt0); ch2->Fill(evt0.realtime);}
        else if (evt0.channel == 3) {ch3->Fill(evt0.realtime);} // the ->Fill() method is used to add a count to the histogram
        else if (evt0.channel == 4) {ch4->Fill(evt0.realtime); ch4_2->Fill(evt0.realtime);}
        else if (evt0.channel == 5) ch5->Fill(evt0.realtime);
        else if (evt0.channel == 6) ch6->Fill(evt0.realtime);
    }
    
    for (long j=0; j<tmcs_1->GetEntries();tmcs_1->GetEntry(j++))
    {
		if (evt1.realtime < start || evt1.realtime > bg_stop) continue;
        if (evt1.channel == 11 or evt1.channel == 12) {PMT1112.push_back(evt1); ch11ch12->Fill(evt1.realtime);}
		else if (evt1.channel == 13 or evt1.channel == 14) {PMT1314.push_back(evt1); ch13ch14->Fill(evt1.realtime);}
		if      (evt1.channel == 11) {ch11->Fill(evt1.realtime);}
        else if (evt1.channel == 12) {ch12->Fill(evt1.realtime);}
        else if (evt1.channel == 13) ch13->Fill(evt1.realtime);
        else if (evt1.channel == 14) ch14->Fill(evt1.realtime);
        else if (evt1.channel == 15) {ch15->Fill(evt1.realtime);hGV->Fill(evt1.realtime);hGV_2->Fill(evt1.realtime);}
        else if (evt1.channel == 16) ch16->Fill(evt1.realtime);
    }
    
    for (long k=0; k<tmcs_2->GetEntries();tmcs_2->GetEntry(k++))
    {
		if (evt2.realtime < start || evt2.realtime > bg_stop) continue;
        if      (evt2.channel == 21) {ch21->Fill(evt2.realtime);}
        else if (evt2.channel == 22) {ch22->Fill(evt2.realtime);}
        else if (evt2.channel == 23) ch23->Fill(evt2.realtime);
        else if (evt2.channel == 24) ch24->Fill(evt2.realtime);
    }

	ofstream somefile;
	string outfile = "./output/coincidences/PECountsRun" + runnum + ".txt";
	std::remove(outfile.c_str());
	somefile.open(outfile, fstream::app);
	auto ti = PMT12.begin();

	while (ti!=PMT12.end()){
		auto PEtime = (*ti).realtime;
		somefile << "12" << ", " << setprecision(15)<< PEtime << ", "<< (*ti).channel << endl;
		ti++;
	}
	auto ti2 = PMT34.begin();
	while (ti2!=PMT34.end()){
		auto PEtime = (*ti2).realtime;
		somefile << "34" << ", " << setprecision(15)<< PEtime << ", "<< (*ti2).channel << endl;
		ti2++;
	}
	auto ti3 = PMT1112.begin();
	while (ti3!=PMT1112.end()){
		auto PEtime = (*ti3).realtime;
		somefile << "56" << ", " << setprecision(15)<< PEtime << ", "<< (*ti3).channel << endl;
		ti3++;
	}
	auto ti4 = PMT1314.begin();
	while (ti4!=PMT1314.end()){
		auto PEtime = (*ti4).realtime;
		somefile << "78" << ", " << setprecision(15)<< PEtime << ", "<< (*ti4).channel << endl;
		ti4++;
	}
	somefile.close();
	
    
    //using pointers to get avg coinc length and number of photons out of the find_coincidences fcn
    double avgku, gammau;
	double avgku12, gammau12;
	double avgku34, gammau34;
	double avgku1112, gammau1112;
	double avgku1314, gammau1314;
    
    //CoincidenceList cl = find_coincidences(PMT_A,PMT_B, 100.0e-9, 1000.0e-9);//20.0e-9,1000.0e-9);//100e-9, 5000e-9);
    //int threshold = 8;
    //find_coincidences(list, coinc window, prompt window, tele window, threshold, avg ku, gammau)
    //was on 25 ns
    CoincidenceList cl12 = find_coincidences(PMT12, wc, wp, wt, threshold, &avgku, &gammau);
	avgku12 = avgku;
	gammau12 = gammau;
	CoincidenceList cl34 = find_coincidences(PMT34, wc, wp, wt, threshold, &avgku, &gammau);
	avgku34 = avgku;
	gammau34 = gammau;
	CoincidenceList cl1112 = find_coincidences(PMT1112, wc, wp, wt, threshold, &avgku, &gammau);
	avgku1112 = avgku;
	gammau1112 = gammau;
	CoincidenceList cl1314 = find_coincidences(PMT1314, wc, wp, wt, threshold, &avgku, &gammau);
	avgku1314 = avgku;
	gammau1314 = gammau;
	
	TH1D* testh = new TH1D("testh", "test; time [s]; counts", int(run_duration)*10,0.0, run_duration);
	TFile* fout = TFile::Open("testfile.root","recreate");
 
	string csv_filename = "./output/coincidences/CoincRun" + runnum + "_" + to_string(threshold) + "PEthreshold.csv";
	std::remove(csv_filename.c_str());

 	//auto [un12, bg12, unp12, bgp12, undt12, bgdt12, unrde12, bgrde12, fillUCN12] = getunloadsbackgrounds(cl12, avgku12, gammau12, run_duration, holdtime, filltime, cleantime, threshold, ch1ch2, 12);
 	//auto [un34, bg34, unp34, bgp34, undt34, bgdt34, unrde34, bgrde34, fillUCN34] = getunloadsbackgrounds(cl34, avgku34, gammau34, run_duration, holdtime, filltime, cleantime, threshold, ch3ch4, 34);
 	//auto [un1112, bg1112, unp1112, bgp1112, undt1112, bgdt1112, unrde1112, bgrde1112, fillUCN1112] = getunloadsbackgrounds(cl1112, avgku1112, gammau1112, run_duration, holdtime, filltime, cleantime, threshold, ch11ch12, 56);
 	//auto [un1314, bg1314, unp1314, bgp1314, undt1314, bgdt1314, unrde1314, bgrde1314, fillUCN1314] = getunloadsbackgrounds(cl1314, avgku1314, gammau1314, run_duration, holdtime, filltime, cleantime, threshold, ch13ch14, 78);
	getunloadsbackgrounds(cl12, runnum, avgku12, gammau12, run_duration, holdtime, filltime, cleantime, threshold, ch1ch2, 12);
	getunloadsbackgrounds(cl34, runnum, avgku34, gammau34, run_duration, holdtime, filltime, cleantime, threshold, ch3ch4, 34);
	getunloadsbackgrounds(cl1112, runnum, avgku1112, gammau1112, run_duration, holdtime, filltime, cleantime, threshold, ch11ch12, 56);
	getunloadsbackgrounds(cl1314, runnum, avgku1314, gammau1314, run_duration, holdtime, filltime, cleantime, threshold, ch13ch14, 78);
	
	delete fin;	
	fout->Close();
	remove("testfile.root");
    return 0;
}
	
int main(int argc, char **argv){
	std::ifstream i("./config/betteroutput2022.json");
	json params;
	i >> params;
   	
   	int startrun = atoi(argv[1]);
    int endrun = atoi(argv[2]);
    int threshold = atoi(argv[3]);
    double wc = atof(argv[4]);
    double wp = atof(argv[5]);
    double wt = atof(argv[6]);
	string run;
	string hold;

	std::set<std::string> good_runs;
	{
		std::ifstream grf("./config/2022runlist.txt");
		std::string line;
		while (std::getline(grf, line)) {
			if (!line.empty()) good_runs.insert(line);
		}
	}

	for (int z =startrun; z<endrun;z++){//was at 23836 to 23936

			string run = std::to_string(z);

			if (good_runs.find(run) == good_runs.end()) {
				cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
				continue;
			}
 
			if (params[run]["run_type"] == "production"){// and params[run]["hold_time"] == 20){
				int holdtime = params[run]["hold_time"];
				int filltime = params[run]["fill_time"];
				int cleantime = params[run]["clean_time"];

				//TH1I* ch4 = ...
				processfile(params[run], run, holdtime, filltime, cleantime, threshold, wc, wp, wt, startrun, endrun);//, pdf_canvas, pdf_canvas2);
				//slopes->Fill(largestvalue);
			}
		//}
	};
	//myfile.close();
	

	return 0;
}
