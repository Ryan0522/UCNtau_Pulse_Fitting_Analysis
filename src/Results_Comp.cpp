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
    const string path = base + "testPECountsPerCoincRun" + std::to_string(run_number) + "_5PEthreshold.csv";

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

// Make quick side‑by‑side ROOT histograms
static void plot_comparisons(const std::vector<PulseRow>& pulse,
                             const std::vector<CoincRow>& coinc,
                             const std::string& out_png_prefix,
                             const json& run_params)
{
    gStyle->SetOptStat(0);
    vector<string> segs = {"12","34","56","78"};

    double start = (double)run_params["fill_time"] + (double)run_params["hold_time"] + (double)run_params["clean_time"] + 40;
    double stop  = start + 60;
	cout << "Start, End time for comparison: " <<start << ", " << stop << endl;

    // --- PE distributions Comparison ---
    {
        TCanvas* c = new TCanvas("c_pe_overlay","PE distributions by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();

            auto seg = segs[i];
            auto h_pulse = new TH1D(("h_pulsePE_"+seg).c_str(),
                                    ("PE (seg "+seg+");PE;Counts").c_str(),
                                    160, 0, 160);

            auto h_coinc = new TH1D(("h_coincPE_"+seg).c_str(),
                                    ("PE (seg"+seg+");PE;Counts").c_str(),
                                    160, 0, 160);

            for (const auto& r : pulse) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= start && t_s <= stop) h_pulse->Fill(std::get<2>(r));

                }
            }
            for (const auto& r : coinc) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= start && t_s <= stop) h_coinc->Fill(std::get<2>(r));

                }
            }
            
            h_pulse->SetLineColor(kBlue+1);  h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);   h_coinc->SetLineWidth(2);

            h_pulse->Draw("HIST");
            h_coinc->Draw("HIST SAME");

            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis PE", "l");
            leg->AddEntry(h_coinc, "Coincidence PE",   "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_PE_overlay.png").c_str());
        delete c;
    }

    // --- PE distributions Comparison (LOG-Y) ---
    {
        TCanvas* c = new TCanvas("c_pe_overlay_log","PE distributions (log y) by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1);
            gPad->SetGrid();
            gPad->SetLogy(1);

            auto seg = segs[i];
            auto h_pulse = new TH1D(("h_pulsePE_log_"+seg).c_str(),
                                    ("PE (seg "+seg+");PE;Counts").c_str(),
                                    160, 0, 160);
            auto h_coinc = new TH1D(("h_coincPE_log_"+seg).c_str(),
                                    ("PE (seg "+seg+");PE;Counts").c_str(),
                                    160, 0, 160);

            // Fill with correct gating & units
            for (const auto& r : pulse) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= start && t_s <= stop) h_pulse->Fill(std::get<2>(r));
                }
            }
            for (const auto& r : coinc) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= start && t_s <= stop) h_coinc->Fill((double)std::get<2>(r));
                }
            }

            // Styles + log-safe minima + shared max
            h_pulse->SetLineColor(kBlue+1);  h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);   h_coinc->SetLineWidth(2);
            h_pulse->SetMinimum(0.5);
            h_coinc->SetMinimum(0.5);
            double ymax = std::max(h_pulse->GetMaximum(), h_coinc->GetMaximum());
            h_pulse->SetMaximum(ymax * 1.2);

            h_pulse->Draw("HIST");
            h_coinc->Draw("HIST SAME");

            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis PE", "l");
            leg->AddEntry(h_coinc, "Coincidence PE",   "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_PE_overlay_log.png").c_str());
        delete c;
    }

    // --- Pulse time distributions (s) (active window only) ---
    {
        TCanvas* c = new TCanvas("c_pulse_t","Pulse Time by segment",1200,800);
        c->Divide(2,2);

        int nbins = 60;
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();

            auto seg = segs[i];
            TH1D* h_pulse = new TH1D(("h_pulseT_"+seg).c_str(),
                               ("Pulse Time (seg "+seg+");Time (s);Counts").c_str(),
                               nbins, start, stop);
            TH1D* h_coinc = new TH1D(("h_coincT_"+seg).c_str(),
                               ("Pulse Time (seg "+seg+");Time (s);Counts").c_str(),
                               nbins, start, stop);

            for (const auto& r : pulse) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= start && t_s <= stop) h_pulse->Fill(t_s);
                }
            }
            for (const auto& r : coinc) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r); // coinc times are in seconds → convert
                    if (t_s >= start && t_s <= stop) h_coinc->Fill(t_s);
                }
            }

            h_pulse->SetLineColor(kBlue+1);  h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);   h_coinc->SetLineWidth(2);

            double max_val = std::max(h_pulse->GetMaximum(), h_coinc->GetMaximum());
            h_pulse->SetMaximum(max_val * 1.2);

            h_pulse->Draw("HIST");
            h_coinc->Draw("HIST SAME");

            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis time", "l");
            leg->AddEntry(h_coinc, "Coincidence time",   "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_Time_overlay.png").c_str());
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

        cout << "Loaded pulse output: " << pulse.size() << " Entries" << endl;

        auto coinc = load_coinc_results(z, output_folder);

        cout << "Loaded coinc output: " << coinc.size() << " Entries" << endl;

        if (pulse.empty() && coinc.empty()) {
            cerr << "[Results_Comp] No data to compare for run " << run << endl;
            continue;
        }

        const string out_prefix = output_folder + "graphs/comparisons/" + run;
        plot_comparisons(pulse, coinc, out_prefix, params[run]);
	}	

	return 0;
}