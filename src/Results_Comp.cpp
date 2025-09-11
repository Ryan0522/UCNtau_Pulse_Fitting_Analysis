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
#include <TROOT.h>
#include <TCollection.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TLine.h>
#include <THStack.h>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

static inline std::string trim(const std::string& s) {
    auto a = s.find_first_not_of(" \t\r\n");
    auto b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

vector<PulseRow> load_pulse_results(int run_number, std::string epoch, const string& output_folder, int ht) {
    vector<PulseRow> rows;
    const string base = ensureTrailingSlash(output_folder);
    const string path = base + "results/epoch_" + epoch + "/PulseAnalysis_" + std::to_string(run_number) + ".csv";

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
        string seg_s, time_s, pe_s, ww_s, fbw_s, evt_s;

        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ','))   continue;
        if (!std::getline(ss, ww_s, ','))   continue;
        if (!std::getline(ss, fbw_s, ','))   continue;
        if (!std::getline(ss, evt_s, ','))   continue;

        seg_s  = trim(seg_s);
        double t_us = std::stod(time_s);
        double pe = std::stod(pe_s);
        double ww = std::stod(ww_s);
        bool fbw = trim(fbw_s) == "1"; // fbw is 1
        bool evt = trim(evt_s) == "1"; // event is 1
        rows.emplace_back(seg_s, t_us, pe, ww, fbw, evt, ht);
    }
    return rows;
}

vector<CoincRow> load_coinc_results(int run_number, std::string epoch, const string& output_folder, int ht) {
    vector<CoincRow> rows;
    const string base = ensureTrailingSlash(output_folder) + "coincidences/";
    const string path = base + "CoincRun" + std::to_string(run_number) + "_5PEthreshold.csv"; // 5 is fixed here (IMPORTANT!!)

    std::ifstream in(path);
    if (!in.is_open()) {
        cerr << "[load_coinc_results] Cannot open: " << path << endl;
        return rows;
    }

    string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        string seg_s, time_s, pe_s, evt_s;
        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ',')) continue;
        if (!std::getline(ss, evt_s, ',')) continue;

        string seg = trim(seg_s);
        double time_val = std::stod(time_s);
        int pe_val = std::stoi(pe_s);
        bool evt = trim(evt_s) == "1";; // event is 1
        rows.emplace_back(seg, time_val, pe_val, evt, ht);
    }
    return rows;
}

std::vector<WindowRow> load_window_stats(int run_number, std::string epoch, const std::string& output_folder) {
    std::vector<WindowRow> rows;
    const string base = ensureTrailingSlash(output_folder) + "results/epoch_" + epoch + "/";
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
    mx = std::min(mx, 1500.0);

    std::map<int, TH1D*> hmap;
    auto mk = [&](int np) {
        auto* h = new TH1D(Form("h_np_%d", np), ";-logL;Windows", nbins, mn, mx);
        h->SetLineWidth(2);
        return h;
    };
    for (const auto& w : ws_all) {
        if (!hmap.count(w.nPulses)) hmap[w.nPulses] = mk(w.nPulses);
        hmap[w.nPulses]->Fill(w.negLogL);
    }

    std::vector<Color_t> cols = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+2, kOrange+7, kViolet+1, kSpring+1};
    int ci = 0;
    for (auto& kv : hmap) kv.second->SetLineColor(cols[ci++ % cols.size()]);

    THStack hs("hs", "NLL by # Pulses; -logL; Windows");
    TLegend leg(0.65, 0.65, 0.88, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0);
    for (auto& [np, h] : hmap) { hs.Add(h, "HIST"); leg.AddEntry(h, Form("%d pulses (n=%d)", np, (int)h->GetEntries()), "l"); }

     // Linear
    TCanvas c1("c1","",1100,700); c1.SetGrid();
    hs.Draw("NOSTACK HIST"); leg.Draw(); c1.Update();
    c1.SaveAs((out_png_prefix + "_nll_byNp.png").c_str());

    // Log-y
    TCanvas c2("c2","",1100,700); c2.SetGrid(); c2.SetLogy(1);
    hs.Draw("NOSTACK HIST"); leg.Draw(); c2.Update();
    c2.SaveAs((out_png_prefix + "_nll_byNp_log.png").c_str());

    for (auto& kv : hmap) delete kv.second;
}

void plot_obs_exp_corr(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix)
{
    if (ws_all.empty()) return;

    TCanvas* c = new TCanvas("c_obs_pdf", "Observed vs PDF", 1000, 650);
    
    // c->SetLogx(1); c->SetLogy(1);
    
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
    g->SetMarkerStyle(6); // tiny dot
    g->Draw("AP"); // A is axis P is point

    // y=x reference line
    double lo = std::min((double)minObs, minExp);
    double hi = std::max((double)maxObs, maxExp);

    hi = std::min(hi, 300.0); // cut-off

    TLine* diag = new TLine(lo, lo, hi, hi);
    diag->SetLineStyle(1);
    diag->SetLineColor(kRed+1);
    diag->Draw("SAME");

    c->Update();
    c->SaveAs((out_png_prefix + "_obs_pdf_correlation.png").c_str());
    delete c;
}

// Make quick side‑by‑side ROOT histograms
static void plot_comparisons(const std::vector<PulseRow>& pulse,
                             const std::vector<CoincRow>& coinc,
                             const std::string& out_png_prefix,
                             const Config cfg,
                             const bool overall=false)
{
    gStyle->SetOptStat(0);
    vector<string> segs = {"12","34","56","78"};

    // ------------------------------------------
    // PE distributions
    // ------------------------------------------
    {
        TCanvas* c = new TCanvas("c_pe_overlay","PE distributions by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid();
            auto seg = segs[i];
            auto h_pulse = new TH1D(("h_pulsePE_"+seg).c_str(), ("PE (seg "+seg+");PE;Counts").c_str(), 160, 0, 160);
            auto h_coinc = new TH1D(("h_coincPE_"+seg).c_str(), ("PE (seg "+seg+");PE;Counts").c_str(), 160, 0, 160);
            for (const auto& r : pulse) if (std::get<0>(r)==seg && std::get<5>(r)) h_pulse->Fill(std::get<2>(r));
            for (const auto& r : coinc) if (std::get<0>(r)==seg && std::get<3>(r)) h_coinc->Fill(std::get<2>(r));
            h_pulse->SetLineColor(kBlue+1);  h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);   h_coinc->SetLineWidth(2);
            h_pulse->Draw("HIST"); h_coinc->Draw("HIST SAME");
            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis PE", "l");
            leg->AddEntry(h_coinc, "LaraCoincidence PE",   "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_PE_overlay.png").c_str());
        delete c;
    }
    {
        TCanvas* c = new TCanvas("c_pe_overlay_log","PE distributions (log y) by segment",1200,800);
        c->Divide(2,2);
        for (size_t i=0;i<segs.size();++i) {
            c->cd((int)i+1); gPad->SetGrid(); gPad->SetLogy(1);
            auto seg = segs[i];
            auto h_pulse = new TH1D(("h_pulsePE_log_"+seg).c_str(), ("PE (seg "+seg+");PE;Counts").c_str(), 160, 0, 160);
            auto h_coinc = new TH1D(("h_coincPE_log_"+seg).c_str(), ("PE (seg "+seg+");PE;Counts").c_str(), 160, 0, 160);
             for (const auto& r : pulse) if (std::get<0>(r)==seg && std::get<5>(r)) h_pulse->Fill(std::get<2>(r));
            for (const auto& r : coinc) if (std::get<0>(r)==seg && std::get<3>(r)) h_coinc->Fill(std::get<2>(r));
            h_pulse->SetLineColor(kBlue+1); h_pulse->SetLineWidth(2);
            h_coinc->SetLineColor(kRed+1);  h_coinc->SetLineWidth(2);
            h_pulse->SetMinimum(0.5); h_coinc->SetMinimum(0.5);
            double ymax = std::max(h_pulse->GetMaximum(), h_coinc->GetMaximum());
            h_pulse->SetMaximum(ymax * 1.2);
            h_pulse->Draw("HIST"); h_coinc->Draw("HIST SAME");
            auto leg = new TLegend(0.60,0.72,0.88,0.90);
            leg->AddEntry(h_pulse, "PulseAnalysis PE", "l");
            leg->AddEntry(h_coinc, "Lara Coincidence PE", "l");
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();
        }
        c->SaveAs((out_png_prefix + "_PE_overlay_log.png").c_str());
        delete c;
    }

    // ------------------------------------------------------------
    // TIME DISTRIBUTIONS
    // ------------------------------------------------------------

    { 
        // helpers
        const std::map<int, double> HT_ANCHOR = {
            {20, 410.0},
            {50, 440.0},
            {100, 490.0},
            {200, 590.0},
            {1550, 1940.0}
        };
        const std::vector<int> HTS = {20, 50, 100, 200, 1550};

        auto anchor_for_ht = [&](int ht)->double {
            auto it = HT_ANCHOR.find(ht);
            return it->second; // no fallback, assume always works (IMPORTANT)
        };

        auto make_hist = [&](const std::string& name)->TH1D* {
            int nbins = (int)std::ceil(60.0 / cfg.bin_width_us);
            auto* h = new TH1D(name.c_str(), ";Time since window start (s); Counts", nbins, 0.0, 60.0);
            h->SetLineWidth(2);
            return h;
        };

        auto draw_ratio_axis = [&](TH1D* hLeft, TH1D* hPulse, TH1D* hCoinc, const std::string& name_suffix)->TH1D* {
            TH1D* h_ratio = (TH1D*)hCoinc->Clone(("h_ratio_"+name_suffix).c_str());
            h_ratio->SetTitle("");
            for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
                double A = hCoinc->GetBinContent(b); // coinc
                double B = hPulse->GetBinContent(b); // pulse
                h_ratio->SetBinContent(b, (A > 0.0) ? (B / A) : 0.0);
                h_ratio->SetBinError(b, 0.0);
            }

            double rightMax = 3.0;
            double scale = (hLeft->GetMaximum() > 0.0) ? (hLeft->GetMaximum() / rightMax) : 1.0;

            TH1D* h_ratio_scaled = (TH1D*)h_ratio->Clone(("h_ratio_scaled_"+name_suffix).c_str());
            h_ratio_scaled->Scale(scale);
            h_ratio_scaled->SetLineColor(kGreen+2);
            h_ratio_scaled->SetLineStyle(2);
            h_ratio_scaled->SetLineWidth(2);
            h_ratio_scaled->Draw("HIST SAME");

            // Right Y axis
            gPad->Update();
            double xRight = gPad->GetUxmax();
            double yLow   = gPad->GetUymin();
            double yHigh  = gPad->GetUymax();
            auto axis = new TGaxis(xRight, yLow, xRight, yHigh, 0.0, rightMax, 510, "+L");
            axis->SetTitle("Pulse/Coinc");
            axis->SetTitleColor(kGreen+2);
            axis->SetLabelColor(kGreen+2);
            axis->SetLineColor(kGreen+2);
            axis->SetLabelSize(0.035);
            axis->SetTitleSize(0.040);
            axis->SetTitleOffset(1.1);
            axis->Draw();
            return h_ratio_scaled; // return for legend
        };

        auto draw_seg_canvas = [&](const std::string& title,
                                   const std::string& fname_suffix,
                                   const std::map<std::string, TH1D*>& pulse_by_seg,
                                   const std::map<std::string, TH1D*>& coinc_by_seg)
        {
            TCanvas* c = new TCanvas(title.c_str(), title.c_str(), 1200, 800);
            c->Divide(2,2);
            for (size_t i = 0; i < segs.size(); ++i) {
                c->cd((int)i+1); gPad->SetGrid();
                gPad->SetRightMargin(0.14); // space for ratio axis

                const auto& seg = segs[i];
                auto itP = pulse_by_seg.find(seg);
                auto itC = coinc_by_seg.find(seg);

                TH1D* hP = itP==pulse_by_seg.end() ? nullptr : itP->second;
                TH1D* hC = itC==coinc_by_seg.end() ? nullptr : itC->second;
                if (!hP && !hC) continue;

                double ymax = 1.0;
                if (hP) ymax = std::max(ymax, hP->GetMaximum());
                if (hC) ymax = std::max(ymax, hC->GetMaximum());

                if (hP) { hP->SetMaximum(ymax*1.2); hP->SetLineColor(kBlue+1); hP->SetLineWidth(2); }
                if (hC) { hC->SetLineColor(kRed+1);  hC->SetLineWidth(2); }

                if (hP && hC) {
                    hP->Draw("HIST");
                    hC->Draw("HIST SAME");
                } else if (hP) {
                    hP->Draw("HIST");
                } else { // hC only
                    hC->SetMaximum(ymax*1.2);
                    hC->Draw("HIST");
                }

                TH1D* hR = nullptr;
                if (hP && hC) {
                    hR = draw_ratio_axis(hP, hP, hC, seg);
                }

                auto leg = new TLegend(0.56,0.70,0.88,0.90);
                if (hP) leg->AddEntry(hP, "PulseAnalysis", "l");
                if (hC) leg->AddEntry(hC, "Lara Coincidence", "l");
                if (hR) leg->AddEntry(hR, "Pulse/Coinc (right)", "l");
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();
            }

            c->SaveAs((out_png_prefix + fname_suffix).c_str());
            delete c;
        };

        using std::optional;
        using std::nullopt;

        // Return t_rel for EVENT pulses if should keep; else nullopt
        auto mk_event_rel_pulse = [&](optional<int> only_ht) {
            return [&, only_ht](const PulseRow& r)->optional<double> {
                if (!std::get<5>(r)) return nullopt;        // event only
                int ht = (int)std::get<6>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<1>(r) - anchor_for_ht(ht);
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt;
                return t_rel;
            };
        };
        // Return t_rel for EVENT coinc if should keep; else nullopt
        auto mk_event_rel_coinc = [&](optional<int> only_ht) {
            return [&, only_ht](const CoincRow& r)->optional<double> {
                if (!std::get<3>(r)) return nullopt;        // event only
                int ht = (int)std::get<4>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<1>(r) - anchor_for_ht(ht);
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt;
                return t_rel;
            };
        };
        // Return t_rel for BG pulses if should keep; else nullopt
        auto mk_bg_rel_pulse = [&](optional<int> only_ht) {
            return [&, only_ht](const PulseRow& r)->optional<double> {
                if (std::get<5>(r)) return nullopt;         // background only
                int ht = (int)std::get<6>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<1>(r) - anchor_for_ht(ht) - 110.0; // 110.0 is bg_start - start
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt; // BG window
                return t_rel;
            };
        };
        // Return t_rel for BG coinc if should keep; else nullopt
        auto mk_bg_rel_coinc = [&](optional<int> only_ht) {
            return [&, only_ht](const CoincRow& r)->optional<double> {
                if (std::get<3>(r)) return nullopt;         // background only
                int ht = (int)std::get<4>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<1>(r) - anchor_for_ht(ht) - 110.0; // 110.0 is bg_start - start
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt; // BG window
                return t_rel;
            };
        };

        auto fill_pair = [&](std::map<std::string, TH1D*>& pulse_by_seg,
                             std::map<std::string, TH1D*>& coinc_by_seg,
                             auto get_trel_pulse,
                             auto get_trel_coinc)
        {
            for (const auto& seg : segs) {
                pulse_by_seg[seg] = make_hist("hP_"+seg);
                coinc_by_seg[seg] = make_hist("hC_"+seg);
            }
            for (const auto& r : pulse) {
                auto t_rel = get_trel_pulse(r);
                if (!t_rel.has_value()) continue;
                const std::string& seg = std::get<0>(r);
                auto it = pulse_by_seg.find(seg);
                if (it != pulse_by_seg.end()) it->second->Fill(*t_rel);
            }
            for (const auto& r : coinc) {
                auto t_rel = get_trel_coinc(r);
                if (!t_rel.has_value()) continue;
                const std::string& seg = std::get<0>(r);
                auto it = coinc_by_seg.find(seg);
                if (it != coinc_by_seg.end()) it->second->Fill(*t_rel);
            }
        };

        // plotting

        if (overall) {
            // 5 per-HT EVENT plots
            for (int ht : HTS) {
                std::map<std::string, TH1D*> p_by_seg, c_by_seg;
                fill_pair(p_by_seg, c_by_seg,
                          mk_event_rel_pulse(ht),
                          mk_event_rel_coinc(ht));
                draw_seg_canvas(
                    Form("Event Time (relative) — HT %d s", ht),
                    Form("_Time_ev_ht%d_rel.png", ht),
                    p_by_seg, c_by_seg
                );
                for (auto& kv : p_by_seg) delete kv.second;
                for (auto& kv : c_by_seg) delete kv.second;
            }
            // +1 OVERALL EVENT
            {
                std::map<std::string, TH1D*> p_by_seg, c_by_seg;
                fill_pair(p_by_seg, c_by_seg,
                          mk_event_rel_pulse(std::nullopt),
                          mk_event_rel_coinc(std::nullopt));
                draw_seg_canvas("Event Time (relative) — OVERALL",
                                "_Time_ev_summed_rel.png",
                                p_by_seg, c_by_seg);
                for (auto& kv : p_by_seg) delete kv.second;
                for (auto& kv : c_by_seg) delete kv.second;
            }
            // +1 OVERALL BACKGROUND
            {
                std::map<std::string, TH1D*> p_by_seg, c_by_seg;
                fill_pair(p_by_seg, c_by_seg,
                          mk_bg_rel_pulse(std::nullopt),
                          mk_bg_rel_coinc(std::nullopt));
                draw_seg_canvas("Background Time (relative) — OVERALL",
                                "_Time_bg_summed_rel.png",
                                p_by_seg, c_by_seg);
                for (auto& kv : p_by_seg) delete kv.second;
                for (auto& kv : c_by_seg) delete kv.second;
            }
        } else {
            // RUN-overall EVENT + BACKGROUND (sum over whatever HTs exist)
            {
                std::map<std::string, TH1D*> p_by_seg, c_by_seg;
                fill_pair(p_by_seg, c_by_seg,
                          mk_event_rel_pulse(std::nullopt),
                          mk_event_rel_coinc(std::nullopt));
                draw_seg_canvas("Event Time (relative) — RUN",
                                "_Time_ev_run_rel.png",
                                p_by_seg, c_by_seg);
                for (auto& kv : p_by_seg) delete kv.second;
                for (auto& kv : c_by_seg) delete kv.second;
            }
            {
                std::map<std::string, TH1D*> p_by_seg, c_by_seg;
                fill_pair(p_by_seg, c_by_seg,
                          mk_bg_rel_pulse(std::nullopt),
                          mk_bg_rel_coinc(std::nullopt));
                draw_seg_canvas("Background Time (relative) — RUN",
                                "_Time_bg_run_rel.png",
                                p_by_seg, c_by_seg);
                for (auto& kv : p_by_seg) delete kv.second;
                for (auto& kv : c_by_seg) delete kv.second;
            }
        }
    }

    // ------------------------------------------------------------
    // WW (unchanged)
    // ------------------------------------------------------------
    
    {
        TCanvas* c = new TCanvas("c_ww_pulse","Pulse Window Width (ww) by segment",1200,800);
        c->Divide(2,2);
        std::map<std::string, std::vector<double>> ww_by_seg;
        double ww_min = std::numeric_limits<double>::infinity();
        double ww_max = -std::numeric_limits<double>::infinity();
        for (const auto& r : pulse) {
            const std::string& seg = std::get<0>(r);
            double ww   = std::get<3>(r);
            const bool evt = std::get<5>(r);
            if (!evt) continue;
            ww_by_seg[seg].push_back(ww);
            ww_min = std::min(ww_min, ww);
            ww_max = std::max(ww_max, ww);
        }
        if (!std::isfinite(ww_min) || !std::isfinite(ww_max)) {
            std::cerr << "[plot_comparisons] No ww entries found in time gate for pulse data.\n";
        } else {
            if (ww_max <= ww_min) ww_max = ww_min + 1.0;
            const int nbins = 50;
            for (size_t i=0;i<segs.size();++i) {
                c->cd((int)i+1); gPad->SetGrid(); gPad->SetLogy(1);
                const auto& seg = segs[i];
                auto h_ww = new TH1D(("h_ww_"+seg).c_str(),
                                      ("Pulse Window Width (seg "+seg+");ww;Counts").c_str(),
                                      nbins, ww_min, ww_max);
                auto it = ww_by_seg.find(seg);
                if (it != ww_by_seg.end()) for (double v : it->second) h_ww->Fill(v);
                h_ww->SetLineColor(kBlue+1); h_ww->SetLineWidth(2); h_ww->Draw("HIST");
                double mean = h_ww->GetMean(), rms = h_ww->GetRMS();
                auto leg = new TLegend(0.60,0.75,0.88,0.90);
                leg->AddEntry(h_ww, "Pulse ww", "l");
                leg->AddEntry((TObject*)0, Form("Mean = %.3g", mean), "");
                leg->AddEntry((TObject*)0, Form("RMS  = %.3g", rms),  "");
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->Draw();
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
    json epoch_json = cfg.epoch_json;
	std::string epoch = std::to_string(cfg.epoch);
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;
	vector<EventList> run_data;
	
    int epoch_start = startrun;
	int epoch_end = endrun;

	if (epoch_json.contains(epoch)) {
		const auto& e = epoch_json[epoch];
		epoch_start = e["start_run_number"].get<int>();
		epoch_end = e["end_run_number"].get<int>();
	}

    string out_prefix;
    vector<PulseRow> all_pulses;
    vector<CoincRow> all_coincs;
    vector<WindowRow> all_windows;

	for (int z =startrun; z<endrun;z++){

		string run = std::to_string(z);

        if (z < epoch_start || z > epoch_end) {
			std::cerr << "Run " << run << " outside epoch " << epoch
                          << " range [" << epoch_start << ", " << epoch_end << "]. Skipping.\n";
            continue;
		}
        
		// skip runs not in good runs list
		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

		if (!params.contains(run) || params[run]["run_type"] != "production") {
			cerr << "Run " << run << " not found or not a production run. Skipping." << endl;
            continue;
        } 

        // --- load ---
        auto pulse = load_pulse_results(z, epoch, output_folder, (int)params[run]["hold_time"]);
        auto ws = load_window_stats(z, epoch, output_folder);
        auto coinc = load_coinc_results(z, epoch, output_folder, (int)params[run]["hold_time"]);

        cout << "\n[Plotting] plotting run " << run << ":" << endl;
        cout << "Loaded pulse output: " << pulse.size() << " Entries" << endl;
        cout << "Loaded coinc output: " << coinc.size() << " Entries" << endl;
        cout << "Loaded window stats: " << ws.size() << " Entries" << endl;

        if (pulse.size() == 0 || coinc.size() == 0 || ws.size() == 0) {
            cerr << "[Results_Comp] No data to compare for run " << run << endl;
            continue;
        } // don't add missing runs

        all_pulses.insert(all_pulses.end(), pulse.begin(), pulse.end());
        all_coincs.insert(all_coincs.end(), coinc.begin(), coinc.end());
        all_windows.insert(all_windows.end(), ws.begin(), ws.end());

        if (cfg.use_coinc) {
            out_prefix = output_folder + "graphs/epoch_" + epoch + "/comp_coinc/individual_run/run_" + run + "/" + run;
            fs::create_directories(fs::path(output_folder + "graphs/epoch_" + epoch + "/comp_coinc/individual_run/run_" + run));
        } else {
            out_prefix = output_folder + "graphs/epoch_" + epoch + "/comp_no_coinc/individual_run/run_" + run + "/" + run;
            fs::create_directories(fs::path(output_folder + "graphs/epoch_" + epoch + "/comp_no_coinc/individual_run/run_" + run));
        }
        plot_comparisons(pulse, coinc, out_prefix, cfg);
        plot_neglog_hist(ws, out_prefix, 100);
        plot_obs_exp_corr(ws, out_prefix);

        gDirectory->GetList()->Clear();       // clear objects in current directory
        gROOT->GetListOfCanvases()->Clear();  // drop canvases->Reset();
	}	

    if (endrun - startrun == 1) return 0;

    cout << "\n[Plotting] plotting overall:" << endl;
    if (cfg.use_coinc) {
        out_prefix = output_folder + "graphs/epoch_" + epoch + "/comp_coinc/overall";
        fs::create_directories(fs::path(output_folder + "graphs/epoch_" + epoch + "/comp_coinc/"));
    } else {
        out_prefix = output_folder + "graphs/epoch_" + epoch + "/comp_no_coinc/overall";
        fs::create_directories(fs::path(output_folder + "graphs/epoch_" + epoch + "/comp_no_coinc/"));
    }
    plot_comparisons(all_pulses, all_coincs, out_prefix, cfg, true);
    plot_neglog_hist(all_windows, out_prefix, 500);
    plot_obs_exp_corr(all_windows, out_prefix);

	return 0;
}