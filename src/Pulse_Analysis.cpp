#include "Pulse_Analysis.h"
#include "Pulse_Fitting.h"
#include "File_Loader.h"
#include "PDF_Global.h"
#include <json.hpp>
#include <iostream>
#include <fstream>
#include <set>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;


static void plot_fit_window(const std::string& out_png,
                            const std::string& title,
                            const std::vector<int>& hist,
                            const std::vector<double>& x,
                            const std::vector<double>& total,
                            const std::vector<std::vector<double>>& comps)
{
    if (hist.empty() || x.size() != hist.size() || total.size() != hist.size())
        return;
	
    gStyle->SetOptStat(0);
	std::string cname = "c_fit_" + out_png;
    TCanvas* c = new TCanvas(cname.c_str(), "Fit window", 1000, 600);
	c->SetMargin(0.12, 0.05, 0.12, 0.08);
	c->SetGrid(1, 1);
	
	// Histogram (observed) with variable bin edges from x (left edges)
	const int nb = (int)hist.size();
	const double dx = (x.size() > 1 ? (x[1] - x[0]) : 1.0);

    // Histogram (observed)
	std::string hname = "h_obs_" + out_png;
    TH1D* h = new TH1D(hname.c_str(), title.c_str(), nb, 0.0, x.back() + dx);
    h->SetDirectory(nullptr); // don't register to gDirectory (August 20, 2025)
	for (int i = 0; i < nb; ++i) h->SetBinContent((int)i+1, hist[i]);
    h->GetXaxis()->SetTitle("Time offset in window (#mu s)");
    h->GetYaxis()->SetTitle("Counts");
	h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->Draw("HIST");

    // Total model
    TGraph* gtot = new TGraph(nb);
	for (int i = 0; i < nb; ++i) gtot->SetPoint(i, x[i], total[i]);
	gtot->SetLineColor(kBlue+1);
	gtot->SetLineWidth(3);
	gtot->Draw("L SAME");

    // Components
    TLegend* leg = new TLegend(0.62, 0.65, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
    leg->AddEntry(h, "Observed", "l");
    leg->AddEntry(gtot, "Model total", "l");

    int colors[6] = {kRed+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet};
    for (size_t k = 0; k < comps.size(); ++k) {
        TGraph* gk = new TGraph((int)comps[k].size());
        for (int i = 0; i < (int)comps[k].size(); ++i) {
            gk->SetPoint(i, x[i], comps[k][i]);
        }
		gk->SetLineStyle(2);
        gk->SetLineWidth(2);
        gk->SetLineColor(colors[k % 6]);
		gk->Draw("L SAME");
        leg->AddEntry(gk, Form("Pulse %d", (int)k+1), "l");
    }
    leg->Draw();

    c->SaveAs(out_png.c_str());

    delete h;
	delete gtot;
	delete leg;
	delete c;
}

// Set up and run the analysis, output to csv file
void analysis_setup(const vector<EventList> run_data, json params, string output_folder, const Config& cfg) { // Event format: <time (us), PE #, event #, window width, # of events in window>
	
	//define signal and background windows (us) from run parameters
	double start = (double)params["fill_time"] + (double)params["hold_time"] + (double)params["clean_time"] + 40;
	double stop = start + 60;
	double bg_start = stop + 50;
	
	vector<string> segment_labels = {"12", "34", "56", "78"};
	vector<int> seg_ids = {12, 34, 56, 78};
	string output_file = output_folder + "results/PulseAnalysis_" + to_string(params["run_number"]) + ".csv";

	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "Error opening output file: " << output_file << endl;
		return;
	}

	out << "Segment, Time (us), PE, Event, Window Width\n";

	for (size_t seg = 0; seg < run_data.size(); ++seg) {
		// run pulse fitting on each segment independently
		// cout << "Segment: " << segment_labels[seg] << endl;
		Pulse_Fitting fitter(run_data[seg]);
		fitter.initFromConfig(cfg);
		fitter.setSegmentId(seg_ids[seg]);
		fitter.setConfigKnobs(cfg.shift_us, cfg.seed_pe_default, cfg.gradient_threshold, cfg.guard_bin);

		fitter.setWindow(start * 1e6, stop * 1e6);
		fitter.setBackgroundWindow(bg_start * 1e6);
		fitter.analyze();

		// plotting (if enabled and available)
		if (cfg.plot_fits) {
			std::string base = ensureTrailingSlash(output_folder) + "graphs/";
			// make sure the directory exists (cheap way; you might already mkdir -p elsewhere)
			system((std::string("mkdir -p ") + base).c_str());

			const auto& sh = fitter.getSingleHist();
			
			if (!sh.empty()) {
				plot_fit_window(base + "first_single_seg" + segment_labels[seg] + ".png",
								"First single-pulse fit (seg " + segment_labels[seg] + ")",
								sh,
								fitter.getSingleX(),
								fitter.getSingleTotal(),
								fitter.getSingleComps());
			}

			const auto& ph = fitter.getPileHist();
			if (!ph.empty()) {
				plot_fit_window(base + "first_pileup_seg" + segment_labels[seg] + ".png",
								Form("First pileup fit (>= %d pulses) (seg %s)", cfg.pileup_min_pulses, segment_labels[seg].c_str()),
								ph,
								fitter.getPileX(),
								fitter.getPileTotal(),
								fitter.getPileComps());
			}
		}

		auto signalPulses = fitter.getSignalPulses();
		auto backgroundPulses = fitter.getBackgroundPulses();

		// write signal pulsese (Event=1)
		for (const auto& event : signalPulses) {
			out << segment_labels[seg] << ", "
				<< get<0>(event)/1e6 << ", "
				<< get<1>(event) << ", "
				<< get<3>(event) << ", "
				<< "1 \n";
		}

		// write background pulses (Event=0)
		for (const auto& event : backgroundPulses) {
			out << segment_labels[seg] << ", "
				<< get<0>(event)/1e6 << ", "
				<< get<1>(event) << ", "
				<< get<3>(event) << ", "
				<< "0 \n";
		}
	}
	
	out.close();
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
    bool        save_to_txt   = cfg.save_to_txt;
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;
	vector<EventList> run_data;
	
	if (save_to_txt) {
		cout << "** Note: converting data to text, no analysis will be performed **" << endl;
	}

	for (int z =startrun; z<endrun;z++){

		string run = std::to_string(z);
		// skip runs not in good runs list
		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

		if (params.contains(run) && params[run]["run_type"] == "production") {
			if (save_to_txt) {
				// convert ROOT -> txt for this run
				processfile(data_folder, output_folder, run);
			} else {
				// analyze this run and write PulseAnalysis_<run>.csv
				run_data = processfile(data_folder, run);
				if (run_data.empty()) {
					cerr << "No data found for run " << run << ". Skipping analysis." << endl;
					continue;
				}
				analysis_setup(run_data, params[run], output_folder, cfg);
			}
		} else {
			cerr << "Run " << run << " not found or not a production run. Skipping." << endl;
		}

	}	

	return 0;
}