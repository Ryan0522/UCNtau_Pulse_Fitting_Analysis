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
#include <TGaxis.h>
#include <TGraph.h>
#include <TLine.h>

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
        string seg_s, time_s, pe_s, ww_s;

        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ','))   continue;
        if (!std::getline(ss, ww_s, ','))   continue;

        seg_s  = trim(seg_s);
        double t_us = std::stod(time_s);
        double pe = std::stod(pe_s);
        double ww = std::stod(ww_s);
        rows.emplace_back(seg_s, t_us, pe, ww);
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

std::vector<WindowRow> load_window_stats(int run_number, const std::string& output_folder) {
    std::vector<WindowRow> rows;
    const string base = ensureTrailingSlash(output_folder) + "results/";
    const string path = base + "PulseWindowStats_" + std::to_string(run_number) + ".csv";

    std::ifstream in(path);
    if (!in.is_open()) {
        cerr << "[load_coinc_results] Cannot open: " << path << endl;
        return rows;
    }

    string line;
    bool is_header = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (is_header) { is_header = false; continue; }

        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> col;
        while (std::getline(ss, tok, ',')) col.push_back(trim(tok));
        if (col.size() < 8) continue;
        WindowRow r;
        r.segment = col[0];
        r.windowIndex = std::stoi(col[1]);
        r.start = std::stod(col[2]);
        r.binWidth_us = std::stod(col[3]);
        r.nPulses = std::stoi(col[4]);
        r.negLogL = std::stod(col[5]);
        r.N_obs = std::stoi(col[6]);
        r.N_exp = std::stod(col[7]);

        rows.push_back(std::move(r));
    }
    return rows;
}

void plot_neglog_hist(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix, int nbins)
{
    if (ws_all.empty()) return;

    double mn = ws_all.front().negLogL, mx = mn;
    for (const auto& w : ws_all) { mn = std::min(mn, w.negLogL); mx = std::max(mx, w.negLogL); }
    if (mx <= mn) mx = mn + 1.0;

    TCanvas* c = new TCanvas("c_neglogL", "Negative Log-Likelihood", 1000, 650);
    TH1D* h = new TH1D("h_neglog", "Negative Log-Likelihood; -logL; Windows", nbins, mn, mx);
    for (const auto& w : ws_all) h->Fill(w.negLogL);

    c->SetGrid();
    h->SetLineWidth(2);
    h->Draw("HIST");
    c->Update();
    c->SaveAs((out_png_prefix + "_neglogL_dist.png").c_str());
    delete c;

    TCanvas* c_log = new TCanvas("c_neglogL_log", "Negative Log-Likelihood (log_)", 1000, 650);
    TH1D* h_log = new TH1D("h_neglog_log", "Negative Log-Likelihood (log); -logL; Windows", nbins, mn, mx);
    for (const auto& w : ws_all) h_log->Fill(w.negLogL);

    c_log->SetGrid();
    c_log->SetLogy(1);
    h_log->SetLineWidth(2);
    h_log->Draw("HIST");
    c_log->Update();
    c_log->SaveAs((out_png_prefix + "_neglogL_dist_log.png").c_str());
    delete c_log;
}

void plot_obs_exp_corr(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix)
{
    if (ws_all.empty()) return;

    TCanvas* c = new TCanvas("c_obs_pdf", "Observed vs PDF", 1000, 650);
    TGraph* g = new TGraph();
    g->SetTitle("Observed vs PDF PE Counts;N_{obs};N_{PDF}");

    int idx = 0;
    int minObs = ws_all.front().N_obs, maxObs = minObs;
    double minExp = ws_all.front().N_exp, maxExp = minExp;

    for (const auto& w : ws_all) {
        g->SetPoint(idx++, (double)w.N_obs, w.N_exp);
        minObs = std::min(minObs, w.N_obs);
        maxObs = std::max(maxObs, w.N_obs);
        minExp = std::min(minExp, w.N_exp);
        maxExp = std::max(maxExp, w.N_exp);
    }

    c->SetGrid();
    g->SetMarkerStyle(20);
    g->Draw("AP");

    // y=x reference line
    double lo = std::min((double)minObs, minExp);
    double hi = std::max((double)maxObs, maxExp);
    TLine* diag = new TLine(lo, lo, hi, hi);
    diag->SetLineStyle(2);
    diag->Draw("SAME");

    c->Update();
    c->SaveAs((out_png_prefix + "_obs_pdf_correlation.png").c_str());
    delete c;
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
    double bg_start = stop + 50;
    double bg_stop = bg_start + 60;
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

    // --- PE distributions Comparison BG (LOG-Y) ---
    {
        TCanvas* c = new TCanvas("c_pe_overlay_log_bg","PE distributions (log y) by segment (bg)",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1);
            gPad->SetGrid();
            gPad->SetLogy(1);

            auto seg = segs[i];
            auto h_pulse = new TH1D(("h_pulsePE_log_bg_"+seg).c_str(),
                                    ("PE (seg "+seg+");PE;Counts").c_str(),
                                    100, 0, 100);
            auto h_coinc = new TH1D(("h_coincPE_log_bg_"+seg).c_str(),
                                    ("PE (seg "+seg+");PE;Counts").c_str(),
                                    100, 0, 100);

            // Fill with correct gating & units
            for (const auto& r : pulse) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= bg_start && t_s <= bg_stop) h_pulse->Fill(std::get<2>(r));
                }
            }
            for (const auto& r : coinc) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= bg_start && t_s <= bg_stop) h_coinc->Fill((double)std::get<2>(r));
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
        c->SaveAs((out_png_prefix + "_PE_overlay_log_bg.png").c_str());
        delete c;
    }

    // --- Pulse time distributions (s) (active window only) ---
    {
        TCanvas* c = new TCanvas("c_pulse_t","Pulse Time by segment",1200,800);
        c->Divide(2,2);

        int nbins = 60; // 1 bin per second since stop = start + 60 (s)

        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            gPad->SetGrid();
            gPad->SetRightMargin(0.14); // leave space for right axis

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

            double leftMax = std::max(h_pulse->GetMaximum(), h_coinc->GetMaximum());
            if (leftMax <= 0) leftMax = 1.0;
            h_pulse->SetMaximum(leftMax * 1.2);

            h_pulse->Draw("HIST");
            h_coinc->Draw("HIST SAME");

            // NEW --- Ratio: Pulse/Coinc per time bin on right axis ---
            auto h_ratio = (TH1D*)h_coinc->Clone(("h_ratio_"+seg).c_str());
            h_ratio->SetTitle("");
            for (int b=1; b<=h_ratio->GetNbinsX(); ++b) {
                double A = h_coinc->GetBinContent(b);
                double B = h_pulse->GetBinContent(b);
                h_ratio->SetBinContent(b, (A > 0.0) ? (B / A) : 0.0);
            }
            double rightMax = h_ratio->GetMaximum(); if (rightMax <= 0) rightMax = 1.0;

            double scale = (h_pulse->GetMaximum()) / rightMax;
            auto h_ratio_scaled = (TH1D*)h_ratio->Clone(("h_ratio_scaled_"+seg).c_str());
            h_ratio_scaled->Scale(scale);
            h_ratio_scaled->SetLineColor(kGreen+2);
            h_ratio_scaled->SetLineStyle(2);
            h_ratio_scaled->SetLineWidth(2);
            h_ratio_scaled->Draw("HIST SAME");

            gPad->Update();
            double xRight = gPad->GetUxmax();
            double yLow = gPad->GetUymin();
            double yHigh = gPad->GetUymax();
            auto axis = new TGaxis(xRight, yLow, xRight, yHigh, 0.0, rightMax, 510, "+L");
            axis->SetTitle("Pulse/Coinc");
            axis->SetTitleColor(kGreen+2);
            axis->SetLabelColor(kGreen+2);
            axis->SetLineColor(kGreen+2);
            axis->SetLabelSize(0.035);
            axis->SetTitleSize(0.040);
            axis->SetTitleOffset(1.1);
            axis->Draw();

            // --- end: NEW (August 26, 2025) ---

            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis time", "l");
            leg->AddEntry(h_coinc, "Coincidence time",   "l");
            leg->AddEntry(h_ratio_scaled, "Pulse/Coinc Ratio (right axis)", "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_Time_overlay.png").c_str());
        delete c;
    }

    // --- Pulse time distributions (s) (bg window only) ---
    {
        TCanvas* c = new TCanvas("c_pulse_t_bg","Pulse Time by segment (bg)",1200,800);
        c->Divide(2,2);

        int nbins = 60; // 1 bin per second since stop = start + 60 (s)

        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            gPad->SetGrid();
            gPad->SetRightMargin(0.14); // leave space for right axis

            auto seg = segs[i];
            TH1D* h_pulse = new TH1D(("h_pulseT_bg_"+seg).c_str(),
                               ("Pulse Time (seg "+seg+");Time (s);Counts").c_str(),
                               nbins, bg_start, bg_stop);
            TH1D* h_coinc = new TH1D(("h_coincT_bg_"+seg).c_str(),
                               ("Pulse Time (seg "+seg+");Time (s);Counts").c_str(),
                               nbins, bg_start, bg_stop);

            for (const auto& r : pulse) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r);
                    if (t_s >= bg_start && t_s <= bg_stop) h_pulse->Fill(t_s);
                }
            }
            for (const auto& r : coinc) {
                if (std::get<0>(r) == seg) {
                    double t_s = std::get<1>(r); // coinc times are in seconds → convert
                    if (t_s >= bg_start && t_s <= bg_stop) h_coinc->Fill(t_s);
                }
            }

            h_pulse->SetLineColor(kBlue+1);  h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);   h_coinc->SetLineWidth(2);

            double leftMax = std::max(h_pulse->GetMaximum(), h_coinc->GetMaximum());
            if (leftMax <= 0) leftMax = 1.0;
            h_pulse->SetMaximum(leftMax * 1.2);

            h_pulse->Draw("HIST");
            h_coinc->Draw("HIST SAME");

            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis time", "l");
            leg->AddEntry(h_coinc, "Coincidence time",   "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_Time_overlay_bg.png").c_str());
        delete c;
    }

        // --- Window Width (pulse-only) distributions by segment ---
    {
        TCanvas* c = new TCanvas("c_ww_pulse","Pulse Window Width (ww) by segment",1200,800);
        c->Divide(2,2);

        // Collect ww values per segment within the gated [start, stop]
        std::map<std::string, std::vector<double>> ww_by_seg;
        double ww_min = std::numeric_limits<double>::infinity();
        double ww_max = -std::numeric_limits<double>::infinity();

        for (const auto& r : pulse) {
            const std::string& seg = std::get<0>(r);
            double t_s  = std::get<1>(r);
            double ww   = std::get<3>(r); // <- Window Width

            if (t_s >= start && t_s <= stop) {
                ww_by_seg[seg].push_back(ww);
                if (ww < ww_min) ww_min = ww;
                if (ww > ww_max) ww_max = ww;
            }
        }

        // Guard: if no ww data was found, skip plotting
        if (!std::isfinite(ww_min) || !std::isfinite(ww_max)) {
            std::cerr << "[plot_comparisons] No ww entries found in time gate for pulse data.\n";
        } else {
            // Make sure the range is reasonable for ROOT (non-zero width)
            if (ww_max <= ww_min) {
                ww_max = ww_min + 1.0;
            }

            const int nbins = 50; // adjust if you want finer/coarser
            for (size_t i = 0; i < segs.size(); ++i) {
                c->cd((int)i+1);
                gPad->SetGrid();

                const auto& seg = segs[i];
                auto h_ww = new TH1D(("h_ww_"+seg).c_str(),
                                     ("Pulse Window Width (seg "+seg+");ww;Counts").c_str(),
                                      nbins, ww_min, ww_max);

                auto it = ww_by_seg.find(seg);
                if (it != ww_by_seg.end()) {
                    for (double v : it->second) h_ww->Fill(v);
                }

                h_ww->SetLineColor(kBlue+1);
                h_ww->SetLineWidth(2);
                h_ww->Draw("HIST");

                // quick stats in the pad
                double mean = h_ww->GetMean();
                double rms  = h_ww->GetRMS();
                auto leg = new TLegend(0.60,0.75,0.88,0.90);
                leg->AddEntry(h_ww, "Pulse ww", "l");
                leg->AddEntry((TObject*)0, Form("Mean = %.3g", mean), "");
                leg->AddEntry((TObject*)0, Form("RMS  = %.3g", rms),  "");
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->Draw();
            }
            c->SaveAs((out_png_prefix + "_WW_pulse.png").c_str());
        }
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

        auto ws = load_window_stats(z, output_folder);

        cout << "Loaded window stats: " << ws.size() << " Entries" << endl;
        
        if (pulse.empty() && coinc.empty()) {
            cerr << "[Results_Comp] No data to compare for run " << run << endl;
            continue;
        }

        const string out_prefix = output_folder + "graphs/comparisons/" + run;
        plot_comparisons(pulse, coinc, out_prefix, params[run]);
        plot_neglog_hist(ws, out_prefix, 200);
        plot_obs_exp_corr(ws, out_prefix);
	}	

	return 0;
}