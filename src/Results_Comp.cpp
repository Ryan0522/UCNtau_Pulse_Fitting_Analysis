#include "Results_Comp.h"
#include "File_Loader.h"
#include "PDF_Global.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <map>
#include <limits>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;

static inline std::string trim(const std::string& s) {
    auto a = s.find_first_not_of(" \t\r\n");
    auto b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

static bool try_open(const std::string& p) {
    std::ifstream f(p);
    return (bool)f;
}

vector<PulseRow> load_pulse_results(int run_number, const string& output_folder) {
    vector<PulseRow> rows;
    const string base = ensureTrailingSlash(output_folder);
    const string path = base + "results/PulseAnalysis_" + std::to_string(run_number) + ".csv";

    std::ifstream in(path);
    if (!in) {
        cerr << "[load_pulse_results] Cannot open: " << path << endl;
        return rows;
    }

    string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // skip header
        if(first) {
            first = false;
            std::string header = trim(line);
            if (header.size() && (header[0] == 'S' || header.find("Segment") != std::string::npos)) {
                // header line — skip and continue
                continue;
            }
        }

        std::stringstream ss(line);
        string seg_s, time_s, pe_s;

        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ','))   continue;

        seg_s  = trim(seg_s);
        double t_us = std::stod(time_s);
        double pe   = std::stod(pe_s);
        rows.emplace_back(seg_s, t_us, pe);
    }
    return rows;
}

vector<CoincRow> load_coinc_results(int run_number, const string& output_folder) {
    vector<CoincRow> rows;
    const string base = ensureTrailingSlash(output_folder) + "coincidences/";
    const string path = base + "PECountsPerCoincRun" + std::to_string(run_number) + "_5PEthresholdtest.csv";

    std::ifstream in(path);
    if (!in.is_open()) {
        cerr << "[load_coinc_results] Cannot open: " << path << endl;
        return rows;
    }

    string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        string seg_s, time_s, pe_s;
        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ',')) continue;

        string seg = trim(seg_s);
        double time_val = std::stod(time_s);
        int pe_val = std::stoi(pe_s);
        rows.emplace_back(seg, time_val, pe_val);
    }
    return rows;
}

// Make quick side‑by‑side ROOT histograms:
//  - Pulse PE (double) per segment
//  - Coincidence N (int) per segment
//  - Pulse time (us) per segment (simple distribution)
static void plot_comparisons(const std::vector<PulseRow>& pulse,
                             const std::vector<CoincRow>& coinc,
                             const std::string& out_png_prefix)
{
    gStyle->SetOptStat(0);

    std::vector<std::string> segs = {"12","34","56","78"};

    // --- Pulse PE distributions (0..300 PE) ---
    {
        TCanvas* c = new TCanvas("c_pulse_pe","Pulse PE by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            TH1D* h = new TH1D(("h_pe_"+segs[i]).c_str(),
                               ("Pulse PE (seg "+segs[i]+");PE;Counts").c_str(),
                               300, 0, 300);
            for (const auto& r : pulse) {
                if (std::get<0>(r) == segs[i]) h->Fill(std::get<2>(r));
            }
            h->SetLineWidth(2);
            h->Draw("HIST");
        }
        c->SaveAs((out_png_prefix + "_pulsePE.png").c_str());
        delete c;
    }

    // --- Coincidence PE‑count distributions (N) ---
    {
        TCanvas* c = new TCanvas("c_coinc_N","Coincidence N by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            TH1D* h = new TH1D(("h_N_"+segs[i]).c_str(),
                               ("Coincidence N (seg "+segs[i]+");N;Counts").c_str(),
                               100, 0, 100);
            for (const auto& r : coinc) {
                if (std::get<0>(r) == segs[i]) h->Fill(std::get<2>(r));
            }
            h->SetLineWidth(2);
            h->Draw("HIST");
        }
        c->SaveAs((out_png_prefix + "_coincN.png").c_str());
        delete c;
    }

    // --- Pulse time distributions (µs) ---
    {
        // find min/max time to choose a reasonable range
        double tmin = std::numeric_limits<double>::infinity();
        double tmax = -std::numeric_limits<double>::infinity();
        for (const auto& r : pulse) {
            tmin = std::min(tmin, std::get<1>(r));
            tmax = std::max(tmax, std::get<1>(r));
        }
        if (!std::isfinite(tmin) || !std::isfinite(tmax) || tmin >= tmax) {
            tmin = 0.0; tmax = 60.0 * 1e6; // default 60 s window in µs
        }
        int nbins = 200;
        TCanvas* c = new TCanvas("c_pulse_t","Pulse Time by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            TH1D* h = new TH1D(("h_t_"+segs[i]).c_str(),
                               ("Pulse Time (seg "+segs[i]+");Time (µs);Counts").c_str(),
                               nbins, tmin, tmax);
            for (const auto& r : pulse) {
                if (std::get<0>(r) == segs[i]) h->Fill(std::get<1>(r));
            }
            h->SetLineWidth(2);
            h->Draw("HIST");
        }
        c->SaveAs((out_png_prefix + "_pulseTime.png").c_str());
        delete c;
    }
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

		if (!params.contains(run) && params[run]["run_type"] == "production") {
			cerr << "Run " << run << " not found or not a production run. Skipping." << endl;
            continue;
        } 

        // --- load ---
        auto pulse = load_pulse_results(z, output_folder);

        cout << "Loaded pulse output" << endl;

        auto coinc = load_coinc_results(z, output_folder);

        cout << "Loaded coinc output" << endl;

        if (pulse.empty() && coinc.empty()) {
            cerr << "[Results_Comp] No data to compare for run " << run << endl;
            continue;
        }

        const string out_prefix = output_folder + "graphs/comp_run" + run;
        plot_comparisons(pulse, coinc, out_prefix);
	}	

	return 0;
}