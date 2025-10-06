#include "Results_Comp.h"
#include "File_Loader.h"
#include "PDF_Global.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cctype>
#include <map>
#include <limits>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TCollection.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TLine.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

std::ostream& operator<<(std::ostream& os, const Config& c) {
	os << "year               = " << c.year << "\n";
    os << "data_folder        = " << c.data_folder << "\n";
    os << "output_folder      = " << c.output_folder << "\n";
    os << "runinfo_path       = " << c.runinfo_path << "\n";
    os << "good_runs_path     = " << c.good_runs_path << "\n";
	os << "coinc_results_path = " << c.coinc_results_path << "\n";
    os << "start_run          = " << c.start_run << "\n";
    os << "end_run            = " << c.end_run << "\n";
    os << "save_to_txt        = " << std::boolalpha << c.save_to_txt << "\n";
	os << "epoch_path         = " << c.epoch_path << "\n";
	os << "epoch              = " << c.epoch << "\n";

    os << "bin_width_us       = " << c.bin_width_us << "\n";
    os << "fine_bin_width_us  = " << c.fine_bin_width_us << "\n";
    os << "min_gap_us         = " << c.min_gap_us << "\n";
    os << "pdf_csv_path       = " << c.pdf_csv_path << "\n";
    os << "good_runs_count    = " << c.good_runs_set.size() << "\n";

	os << "seeding_window_us  = " << c.seeding_window_us << "\n";
    os << "pe_min_thresh      = " << c.pe_min_thresh << "\n";

    os << "shift_us           = " << c.shift_us << "\n";
    os << "seed_pe_default    = " << c.seed_pe_default << "\n";
    os << "gradient_threshold = " << c.gradient_threshold << "\n";
    os << "guard_bin          = " << c.guard_bin << "\n";
	os << "cluster_close_us   = " << c.cluster_close_us << "\n";

    os << "plot_fits          = " << std::boolalpha << c.plot_fits << "\n";
    os << "pileup_min_pulses  = " << c.pileup_min_pulses << "\n";

	os << "plot_outliers      = " << std::boolalpha << c.plot_outliers << "\n";
    os << "outlier_min_obs    = " << c.outlier_min_obs << "\n";
    os << "outlier_ratio_low  = " << c.outlier_ratio_low << "\n";

	os << "use_coinc          = " << std::boolalpha << c.use_coinc << "\n";
	os << "coinc_win_us       = " << c.coinc_win_us << "\n";
	os << "coinc_seed_pe_min  = " << c.coinc_seed_pe_min << "\n";

	os << "debug              = " << c.debug << "\n";
	os << "debug_window_index = " << c.debug_window_index << "\n";
	os << "debug_segment_id   = " << c.debug_segment_id << "\n";

    return os;
}

static inline std::string fmt_r(double r) {
    if (std::isnan(r)) return "n/a";
    char buf[64]; std::snprintf(buf, sizeof(buf), "%.7f", r); // 7 digits after decimal
    return std::string(buf);
}

static inline double pearson_r(const TGraph* g) {
    const int n = g ? g->GetN() : 0;
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    double sx=0, sy=0, x, y;
    for (int i=0;i<n;++i) { g->GetPoint(i, x, y); sx += x; sy += y; }
    const double mx = sx / n, my = sy / n;
    double sxx=0, syy=0, sxy=0;
    for (int i=0;i<n;++i) { 
        g->GetPoint(i, x, y);
        const double dx = x - mx, dy = y - my;
        sxx += dx*dx; syy += dy*dy; sxy += dx*dy;
    }
    if (sxx <= 0 || syy <= 0) return std::numeric_limits<double>::quiet_NaN();
    return sxy / std::sqrt(sxx * syy);
}

static inline double pearson_r_vec(const std::vector<double>& xs, const std::vector<double>& ys) {
    const int n = (int)std::min(xs.size(), ys.size());
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    double sx=0, sy=0;
    for (int i=0;i<n;++i) { sx += xs[i]; sy += ys[i]; }
    const double mx = sx / n, my = sy / n;
    double sxx=0, syy=0, sxy=0;
    for (int i=0;i<n;++i) {
        const double dx = xs[i] - mx, dy = ys[i] - my;
        sxx += dx*dx; syy += dy*dy; sxy += dx*dy;
    }
    if (sxx <= 0 || syy <= 0) return std::numeric_limits<double>::quiet_NaN();
    return sxy / std::sqrt(sxx * syy);
}

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
        string run_s, seg_s, time_s, pe_s, ww_s, fbw_s, evt_s;

        if (!std::getline(ss, run_s, ',')) continue;
        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ','))   continue;
        if (!std::getline(ss, ww_s, ','))   continue;
        if (!std::getline(ss, fbw_s, ','))   continue;
        if (!std::getline(ss, evt_s, ','))   continue;

        int run = std::stoi(run_s);
        seg_s = trim(seg_s);
        double t_us = std::stod(time_s);
        double pe = std::stod(pe_s);
        double ww = std::stod(ww_s);
        bool fbw = trim(fbw_s) == "1"; // fbw is 1
        bool evt = trim(evt_s) == "1"; // event is 1
        rows.emplace_back(run, seg_s, t_us, pe, ww, fbw, evt, ht);
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
        string run_s, seg_s, time_s, pe_s, evt_s;
        if (!std::getline(ss, run_s, ',')) continue;
        if (!std::getline(ss, seg_s, ',')) continue;
        if (!std::getline(ss, time_s, ',')) continue;
        if (!std::getline(ss, pe_s, ',')) continue;
        if (!std::getline(ss, evt_s, ',')) continue;

        int run = std::stoi(run_s);
        string seg = trim(seg_s);
        double time_val = std::stod(time_s);
        int pe_val = std::stoi(pe_s);
        bool evt = trim(evt_s) == "1";; // event is 1
        rows.emplace_back(run, seg, time_val, pe_val, evt, ht);
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
        r.run = run_number;
        r.segment = col[0];
        r.windowIndex = std::stoi(col[1]);
        r.start = std::stod(col[2]);
        r.end = std::stod(col[3]);
        r.binWidth_us = std::stod(col[4]);
        r.nPulses = std::stoi(col[5]);
        r.negLogL = std::stod(col[6]);
        r.N_obs = std::stoi(col[7]);
        r.N_exp = std::stod(col[8]);

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
        auto* h = new TH1D(Form("h_np_%d", np), ";-logL;Event #", nbins, mn, mx);
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

    THStack hs("hs", "NLL by # Pulses; -logL; Event #");
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

void plot_window_sep(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix) 
{
    if (ws_all.size() < 2) return;

    std::map<std::string, std::vector<WindowRow>> by_seg;
    for (const auto& w : ws_all) by_seg[w.segment].push_back(w);
    for (auto& kv : by_seg) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end(), [](const WindowRow& a, const WindowRow& b){ return a.start < b.start; });
    }

    static int uid = 0; ++uid;
    TCanvas* c = new TCanvas(Form("c_winsep_seg_%d", uid), "Window separation by segment", 1200, 900);
    c->Divide(2,2);

    std::vector<std::string> segs = {"12","34","56","78"};
    int pad = 1;

    for (const auto& seg : segs) {
        c->cd(pad++); 
        gPad->SetGrid();
        gPad->SetLogx();

        const auto it = by_seg.find(seg);
        if (it == by_seg.end()) {
            continue;
        }
        const auto& W = it->second;
        if (W.size() <= 1) continue;

        std::vector<double> seps; seps.reserve(W.size()-1);
        double mn, mx;
        {
            double sep0 = W[1].start - W[0].end;
            mn = mx = sep0;
            for (size_t i = 0; i + 1 < W.size(); ++i) {
                if (W[i+1].run != W[i].run) continue;
                double sep = W[i+1].start - W[i].end;
                if (sep < 0) std::cout << "[winsep] Overlap in seg " << seg
                                       << " at windowIndex " << W[i].windowIndex << "\n";
                seps.push_back(sep);
                mn = std::max(std::min(mn, sep), 1e-9);
                mx = std::max(mx, sep);
            }
        }

        // build log-spaced bin edges
        int nbins = 100; // or tune
        double logmin = std::log10(mn * 0.9);
        double logmax = std::log10(mx * 1.1);
        std::vector<double> edges(nbins+1);
        for (int i=0; i<=nbins; ++i) {
            edges[i] = std::pow(10, logmin + (logmax-logmin)*i/nbins);
        }

        TH1D* h = new TH1D(Form("h_winsep_%s_%d", seg.c_str(), uid),
                           Form("seg %s (min #Delta t = %.3g s);#Delta t = start_{i+1} - end_{i} (s);#windows (log)", 
                                seg.c_str(), mn),
                           nbins, edges.data());
        h->SetLineWidth(2);
        h->GetXaxis()->SetMoreLogLabels();
        h->GetXaxis()->SetNoExponent();
        for (double v : seps) h->Fill(v);
        h->Draw("HIST");
    }

    c->Modified(); c->Update();
    c->SaveAs((out_png_prefix + "_window_separation_bySeg.png").c_str());
    delete c;
}

void plot_time_ev_byNp(const std::vector<PulseRow>& pulses,
                                      const std::vector<WindowRow>& ws_all,
                                      const std::string& out_png_prefix,
                                      const Config& cfg)
{
    using std::string; using std::vector; using std::map; using std::pair;

    const std::map<int, double> HT_ANCHOR = {
        {20, 410.0}, {50, 440.0}, {100, 490.0}, {200, 590.0}, {1550, 1940.0}
    };
    const std::vector<int> HTS = {20, 50, 100, 200, 1550};
    auto anchor_for_ht = [&](int ht)->double { return HT_ANCHOR.at(ht); };

    struct W { double start_s, end_s; int np; };
    std::map<std::pair<std::string, int>, std::vector<W>> windows;

    for (const auto& w : ws_all) {
        for (int ht : HTS) {
            double a = anchor_for_ht(ht);
            if (w.start >= a && w.start < a + 60.0) {
                windows[{w.segment, ht}].push_back({w.start, 0.0, w.nPulses});
                break;
            }
        }
    }

    for (auto& kv : windows) {
        auto& vec = kv.second;
        std::sort(vec.begin(), vec.end(), [](const W& a, const W& b){ return a.start_s < b.start_s; });
        int ht = kv.first.second;
        double cap = anchor_for_ht(ht) + 60.0;
        for (size_t i=0;i<vec.size();++i) {
            double end = (i+1<vec.size()? vec[i+1].start_s : cap);
            if (end <= vec[i].start_s)
                end = std::nextafter(vec[i].start_s, std::numeric_limits<double>::infinity());
            if (end > cap) end = cap;
            vec[i].end_s = end;
        }
    }

    const std::vector<std::string> segs = {"12","34","56","78"};
    auto mk_time_hist = [&](const std::string& name)->TH1D* {
        int nbins = (int)std::ceil(60.0 / cfg.bin_width_us);
        auto* h = new TH1D(name.c_str(), ";Time since window start (s);Counts", nbins, 0.0, 60.0);
        h->SetLineWidth(2);
        return h;
    };

    // Prepare histograms (1,2,3,4,5+) per segment
    std::map<std::string, TH1D*> h_np1, h_np2, h_np3, h_np4, h_np5p;
    for (const auto& seg : segs) {
        h_np1[seg]  = mk_time_hist("h_np1_"+seg);
        h_np2[seg]  = mk_time_hist("h_np2_"+seg);
        h_np3[seg]  = mk_time_hist("h_np3_"+seg);
        h_np4[seg]  = mk_time_hist("h_np4_"+seg);
        h_np5p[seg] = mk_time_hist("h_np5p_"+seg);
    }

    // Fill from EVENT pulses, aligned to their HT anchor
    for (const auto& r : pulses) {
        if (!std::get<6>(r)) continue;  // Event only
        const std::string& seg = std::get<1>(r);
        int ht = (int)std::get<7>(r);
        auto itA = HT_ANCHOR.find(ht);
        if (itA == HT_ANCHOR.end()) continue;

        double t_abs = std::get<2>(r);
        double t_rel = t_abs - itA->second;
        if (t_rel < 0.0 || t_rel >= 60.0) continue;

        // find window multiplicity
        int np = 0;
        auto key = std::make_pair(seg, ht);
        auto itW = windows.find(key);
        if (itW != windows.end()) {
            const auto& vec = itW->second;
            int L=0, R=(int)vec.size();
            while (L < R) {
                int M=(L+R)/2;
                if (vec[M].start_s <= t_abs) L=M+1; else R=M;
            }
            int idx = L-1;
            if (idx >= 0) {
                const auto& w = vec[(size_t)idx];
                if (t_abs >= w.start_s && t_abs < w.end_s) np = w.np;
            }
        }

        if (np <= 1)      h_np1[seg]->Fill(t_rel);
        else if (np == 2) h_np2[seg]->Fill(t_rel);
        else if (np == 3) h_np3[seg]->Fill(t_rel);
        else if (np == 4) h_np4[seg]->Fill(t_rel);
        else              h_np5p[seg]->Fill(t_rel);
    }

    // Draw
    TCanvas* c = new TCanvas("c_time_byNp_SUM", "Event time by #pulses (SUMMED)", 1200, 800);
    c->Divide(2,2);
    for (size_t i=0;i<segs.size();++i) {
        c->cd((int)i+1); gPad->SetGrid();
        auto seg = segs[i];
        auto* h1 = h_np1[seg];
        auto* h2 = h_np2[seg];
        auto* h3 = h_np3[seg];
        auto* h4 = h_np4[seg];
        auto* h5 = h_np5p[seg];

        double ymax = 1.0;
        ymax = std::max(ymax, h1->GetMaximum());
        ymax = std::max(ymax, h2->GetMaximum());
        ymax = std::max(ymax, h3->GetMaximum());
        ymax = std::max(ymax, h4->GetMaximum());
        ymax = std::max(ymax, h5->GetMaximum());
        h1->SetMaximum(ymax*1.2);

        h1->SetLineColor(kBlue+1);
        h2->SetLineColor(kRed+1);
        h3->SetLineColor(kGreen+2);
        h4->SetLineColor(kMagenta+1);
        h5->SetLineColor(kCyan+2);

        h1->Draw("HIST");
        h2->Draw("HIST SAME");
        h3->Draw("HIST SAME");
        h4->Draw("HIST SAME");
        h5->Draw("HIST SAME");

        auto leg = new TLegend(0.60,0.72,0.88,0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(h1, "1 pulse", "l");
        leg->AddEntry(h2, "2 pulses", "l");
        leg->AddEntry(h3, "3 pulses", "l");
        leg->AddEntry(h4, "4 pulses", "l");
        leg->AddEntry(h5, "5+ pulses", "l");
        leg->Draw();
    }
    c->SaveAs((out_png_prefix + "_Time_ev_byNp_SUM_rel.png").c_str());
    delete c;

    for (auto& kv : h_np1)  delete kv.second;
    for (auto& kv : h_np2)  delete kv.second;
    for (auto& kv : h_np3)  delete kv.second;
    for (auto& kv : h_np4)  delete kv.second;
    for (auto& kv : h_np5p) delete kv.second;
}

void plot_obs_exp_corr(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix)
{
    if (ws_all.empty()) return;

    std::vector<std::string> segs = {"12", "34", "56", "78"};
    std::map<std::string, TGraph*> g_by_seg;
    for (auto& s : segs) g_by_seg[s] = new TGraph();

    int minObs = ws_all.front().N_obs;
    double minExp = ws_all.front().N_exp;
    const double maxObs = 220.0;
    const double maxExp = 220.0;
    
    for (const auto& w : ws_all) {
        auto it = g_by_seg.find(w.segment);
        if (it == g_by_seg.end()) continue;
        int p = it->second->GetN();
        it->second->SetPoint(p, (double)w.N_obs, w.N_exp);
        minObs = std::min(minObs, w.N_obs);
        // maxObs = std::max(maxObs, w.N_obs);
        minExp = std::min(minExp, w.N_exp);
        // maxExp = std::max(maxExp, w.N_exp);
    }

    double lo = std::min((double)minObs, minExp);
    if (lo > 0.0) lo = 0.0;
    double hi = std::max((double)maxObs, maxExp);
    if (hi <= lo) hi = lo + 1.0;

    TCanvas* c = new TCanvas("c_obs_pdf_byseg", "Observed vs PDF by segment", 1200, 800);
    c->Divide(2,2);
    int colors[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};

    for (size_t i = 0; i < segs.size(); ++i) {
        c->cd((int)i+1); gPad->SetGrid();
        gPad->DrawFrame(
            lo, lo, hi, hi,
            Form("Observed vs PDF - seg %s;N_{obs};N_{PDF}", segs[i].c_str())
        );

        auto* g = g_by_seg[segs[i]];
        g->SetMarkerStyle(6);
        g->SetMarkerSize(0.7);
        g->SetMarkerColor(colors[i]);
        g->Draw("P SAME");

        auto* diag = new TLine(lo, lo, hi, hi);
        diag->SetLineColor(kGray+2);
        diag->SetLineStyle(2);
        diag->SetLineWidth(2);
        diag->Draw("SAME");

        c->Modified(); c->Update();
    }
    c->SaveAs((out_png_prefix + "_obs_pdf_correlation.png").c_str());
    delete c;

    // Overlay
    TCanvas* c2 = new TCanvas("c_obs_pdf_overlay_byseg", "Observed vs PDF (overlay by segment)", 1000, 700);
    c2->SetGrid();
    TH2D frameAll("fr_all", "Observed vs PDF — all segments;N_{obs};N_{PDF}", 1, lo, hi, 1, lo, hi);
    frameAll.SetStats(0); frameAll.Draw();

    TLegend leg(0.62, 0.16, 0.88, 0.34);
    leg.SetBorderSize(0); leg.SetFillStyle(0);

    for (size_t i = 0; i < segs.size(); ++i) {
        auto* g = g_by_seg[segs[i]];
        g->SetMarkerStyle(6);
        g->SetMarkerSize(0.6);
        g->SetMarkerColor(colors[i]);
        g->Draw("P SAME");
        leg.AddEntry(g, Form("seg %s (n=%d)", segs[i].c_str(), g->GetN()), "p");
    }

    TLine diag2(lo, lo, hi, hi);
    diag2.SetLineColor(kGray+2);
    diag2.SetLineStyle(2);
    diag2.Draw("SAME");
    leg.Draw();

    c2->SaveAs((out_png_prefix + "_obs_pdf_correlation_overlay.png").c_str());
    delete c2;

    for (auto& kv : g_by_seg) delete kv.second;
}

// Make quick side‑by‑side ROOT histograms
void plot_comparisons(const std::vector<PulseRow>& pulse,
                             const std::vector<CoincRow>& coinc,
                             const std::string& out_png_prefix,
                             const Config cfg,
                             const bool overall)
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
            for (const auto& r : pulse) if (std::get<1>(r)==seg && std::get<6>(r)) h_pulse->Fill(std::get<3>(r));
            for (const auto& r : coinc) if (std::get<1>(r)==seg && std::get<4>(r)) h_coinc->Fill(std::get<3>(r));
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
             for (const auto& r : pulse) if (std::get<1>(r)==seg && std::get<6>(r)) h_pulse->Fill(std::get<3>(r));
            for (const auto& r : coinc) if (std::get<1>(r)==seg && std::get<4>(r)) h_coinc->Fill(std::get<3>(r));
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
            double rightMin = 0.5, rightMax = 1.3;
            hPulse->Sumw2(); hCoinc->Sumw2();

            auto* h_ratio = (TH1D*)hPulse->Clone(("h_ratio_"+name_suffix).c_str());
            h_ratio->Divide(hCoinc);

            gPad->Update();
            double yLow  = gPad->GetUymin();
            double yHigh = gPad->GetUymax();
            double m = (yHigh - yLow) / (rightMax - rightMin);   // slope
            double b = yLow - m * rightMin;                      // intercept

            auto* h_ratio_scaled = (TH1D*)h_ratio->Clone(("h_ratio_scaled_"+name_suffix).c_str());
            for (int bni=1; bni<=h_ratio_scaled->GetNbinsX(); ++bni) {
                double r = h_ratio->GetBinContent(bni);         // true ratio value
                h_ratio_scaled->SetBinContent(bni, m*r + b);    // draw in left-axis coords
                h_ratio_scaled->SetBinError(bni, 0.0);
            }

            h_ratio_scaled->SetLineColor(kGreen+2);
            h_ratio_scaled->SetLineStyle(2);
            h_ratio_scaled->SetLineWidth(2);
            h_ratio_scaled->Draw("HIST SAME");

            // Right Y axis (now consistent)
            auto axis = new TGaxis(gPad->GetUxmax(), yLow, gPad->GetUxmax(), yHigh,
                                rightMin, rightMax, 510, "+L");
            axis->SetTitle("Pulse/Coinc");
            axis->SetTitleColor(kGreen+2);
            axis->SetLabelColor(kGreen+2);
            axis->SetLineColor(kGreen+2);
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
                if (!std::get<6>(r)) return nullopt;        // event only
                int ht = (int)std::get<7>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<2>(r) - anchor_for_ht(ht);
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt;
                return t_rel;
            };
        };
        // Return t_rel for EVENT coinc if should keep; else nullopt
        auto mk_event_rel_coinc = [&](optional<int> only_ht) {
            return [&, only_ht](const CoincRow& r)->optional<double> {
                if (!std::get<4>(r)) return nullopt;        // event only
                int ht = (int)std::get<5>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<2>(r) - anchor_for_ht(ht);
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt;
                return t_rel;
            };
        };
        // Return t_rel for BG pulses if should keep; else nullopt
        auto mk_bg_rel_pulse = [&](optional<int> only_ht) {
            return [&, only_ht](const PulseRow& r)->optional<double> {
                if (std::get<6>(r)) return nullopt;         // background only
                int ht = (int)std::get<7>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<2>(r) - anchor_for_ht(ht) - 110.0; // 110.0 is bg_start - start
                if (t_rel < 0.0 || t_rel >= 60.0) return nullopt; // BG window
                return t_rel;
            };
        };
        // Return t_rel for BG coinc if should keep; else nullopt
        auto mk_bg_rel_coinc = [&](optional<int> only_ht) {
            return [&, only_ht](const CoincRow& r)->optional<double> {
                if (std::get<4>(r)) return nullopt;         // background only
                int ht = (int)std::get<5>(r);
                if (only_ht && ht != *only_ht) return nullopt;
                double t_rel = std::get<2>(r) - anchor_for_ht(ht) - 110.0; // 110.0 is bg_start - start
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
                const std::string& seg = std::get<1>(r);
                auto it = pulse_by_seg.find(seg);
                if (it != pulse_by_seg.end()) it->second->Fill(*t_rel);
            }
            for (const auto& r : coinc) {
                auto t_rel = get_trel_coinc(r);
                if (!t_rel.has_value()) continue;
                const std::string& seg = std::get<1>(r);
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
            const std::string& seg = std::get<1>(r);
            double ww   = std::get<4>(r);
            const bool evt = std::get<6>(r);
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

static inline bool in_window(double t_s, double start_s, double end_s) {
    return (t_s >= start_s - 1e-7) && (t_s < end_s); // start time within 5 us of my start_s
}

void plot_pulses_vs_coinc_by_window(const std::vector<PulseRow>& pulses,
                                    const std::vector<CoincRow>& coincs,
                                    const std::vector<WindowRow>& windows,
                                    const std::string& out_png_prefix)
{
    if (windows.empty()) return;
    
    using RunSeg = std::pair<int, std::string>; // {run_id, "12"/"34"/...}

    std::map<RunSeg, std::vector<WindowRow>> win_by_runseg;
    for (const auto& w : windows) {
        int run_id = w.run;
        win_by_runseg[{run_id, w.segment}].push_back(w);
    }
    for (auto& kv : win_by_runseg) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end(), 
                    [](const WindowRow& a, const WindowRow& b){ return a.start < b.start; });
    }

    struct PT { int nPulse=0, nCoinc=0; };
    std::map<RunSeg, std::vector<PT>> per_runseg_pts;

    for (auto& kv : win_by_runseg) {
        const RunSeg& rs = kv.first;
        const auto& W = kv.second;
        per_runseg_pts[rs].assign(W.size(), PT{});
    }

    // Pulses
    for (const auto& r : pulses){
        if (!std::get<6>(r)) continue; // Event only
        int run_id = std::get<0>(r);
        const std::string& s = std::get<1>(r);
        double t = std::get<2>(r);

        RunSeg rs{run_id, s};
        auto itW = win_by_runseg.find(rs);
        if (itW == win_by_runseg.end()) continue;
        const auto& W = itW->second;

        int L=0, R=(int)W.size();
        while (L<R) { int M=(L+R)/2; if (W[M].start <= t) L=M+1; else R=M; }
        int i = L-1;
        if (i>=0 && i<(int)W.size() && in_window(t, W[i].start, W[i].end))
            per_runseg_pts[rs][i].nPulse++;
    }

    // Coinc
    for (const auto& r : coincs) {
        if (!std::get<4>(r)) continue; // Event only
        int run_id = std::get<0>(r);
        const std::string& s = std::get<1>(r);
        double t = std::get<2>(r);

        RunSeg rs{run_id, s};
        auto itW = win_by_runseg.find(rs);
        if (itW == win_by_runseg.end()) continue;
        const auto& W = itW->second;

        int L=0, R=(int)W.size();
        while (L<R) { int M=(L+R)/2; if (W[M].start <= t) L=M+1; else R=M; }
        int i = L-1;

        auto in_win = [&](int k)->bool {
            if (k < 0 || k >= (int)W.size()) return false;
            return in_window(t, W[k].start, W[k].end);
        };

        if      (in_win(i))   per_runseg_pts[rs][i].nCoinc++;
        else if (in_win(i+1)) per_runseg_pts[rs][i+1].nCoinc++;
    }

    std::map<std::string, std::vector<PT>> per_seg_pts;
    for (const auto& kv : per_runseg_pts) {
        const std::string& seg = kv.first.second;
        const auto& vec = kv.second;
        per_seg_pts[seg].insert(per_seg_pts[seg].end(), vec.begin(), vec.end());
    }

    // by segment
    {
        TCanvas* c = new TCanvas("c_pulses_vs_coinc_bywin_hists", "Pulses-per-window hists (per segment)", 1200, 800);
        c->Divide(2,2);

        std::vector<std::string> segs = {"12","34","56","78"};
        int colMine = kBlue+1;
        int colCoin = kRed+1;

        for (size_t i=0;i<segs.size();++i) {
            const auto& seg = segs[i];
            auto it = per_seg_pts.find(seg);
            if (it == per_seg_pts.end()) continue;

            // find max bin
            int maxCnt = 1;
            for (const auto& p : it->second) {
                maxCnt = std::max(maxCnt, std::max(p.nPulse, p.nCoinc));
            }

            // build histograms of "count per window"
            int nbins = maxCnt + 1; // bins for 0..maxCnt
            TH1I* hMine = new TH1I(Form("hMine_%s", seg.c_str()),
                                Form("seg %s;#pulses per window;#windows", seg.c_str()),
                                nbins, -0.5, maxCnt + 0.5);
            TH1I* hCoin = new TH1I(Form("hCoin_%s", seg.c_str()),
                                Form("seg %s;#pulses per window;#windows", seg.c_str()),
                                nbins, -0.5, maxCnt + 0.5);

            for (const auto& p : it->second) {
                hMine->Fill(p.nPulse);
                hCoin->Fill(p.nCoinc);
            }

            hMine->SetLineColor(colMine);
            hMine->SetLineWidth(2);
            hCoin->SetLineColor(colCoin);
            hCoin->SetLineWidth(2);

            c->cd((int)i+1); gPad->SetGrid();
            // draw with y-range that fits both
            int ymax = std::max(hMine->GetMaximum(), hCoin->GetMaximum());
            hMine->SetMaximum(ymax * 1.15);
            hMine->Draw("HIST");
            hCoin->Draw("HIST SAME");

            TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(hMine, "my pulses", "l");
            leg->AddEntry(hCoin, "coinc pulses", "l");
            leg->Draw();
        }

        c->SaveAs((out_png_prefix + "_byWindow_myVsCoinc.png").c_str());
        delete c;
    }

    // overall combined (all segments merged) hist overlay: my vs coinc
    {
        TCanvas* c = new TCanvas("c_pulses_vs_coinc_bywin_overlay_hists", "Pulses-per-window hists (overall)", 1000, 700);
        c->SetGrid();

        std::vector<int> allMine, allCoinc;
        int maxCnt = 1;
        for (const auto& kv : per_seg_pts) {
            for (const auto& p : kv.second) {
                allMine.push_back(p.nPulse);
                allCoinc.push_back(p.nCoinc);
                maxCnt = std::max(maxCnt, std::max(p.nPulse, p.nCoinc));
            }
        }

        int nbins = maxCnt + 1;
        TH1I* hMine = new TH1I("hMine_all", "All segments;#pulses per window;#windows",
                            nbins, -0.5, maxCnt + 0.5);
        TH1I* hCoin = new TH1I("hCoin_all", "All segments;#pulses per window;#windows",
                            nbins, -0.5, maxCnt + 0.5);

        for (int v : allMine) hMine->Fill(v);
        for (int v : allCoinc) hCoin->Fill(v);

        int colMine = kBlue+1, colCoin = kRed+1;
        hMine->SetLineColor(colMine); hMine->SetLineWidth(3);
        hCoin->SetLineColor(colCoin); hCoin->SetLineWidth(3);

        int ymax = std::max(hMine->GetMaximum(), hCoin->GetMaximum());
        hMine->SetMaximum(ymax * 1.20);

        hMine->Draw("HIST");
        hCoin->Draw("HIST SAME");

        TLegend leg(0.62, 0.72, 0.88, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(hMine, "my pulses (all segs)", "l");
        leg.AddEntry(hCoin, "coinc pulses (all segs)", "l");
        leg.Draw();

        c->SaveAs((out_png_prefix + "_byWindow_myVsCoinc_overlay.png").c_str());
        delete c;
    }
}

static bool load_RDE_table(const std::string& path,
                           std::map<int, std::map<std::string, double>>& rde_by_run // run -> { "12":val, ... }
                           )
{
    std::ifstream in(path);
    if (!in) { std::cerr << "[RDE] cannot open: " << path << "\n"; return false; }

    std::string line;
    if (!std::getline(in, line)) { std::cerr << "[RDE] empty file: " << path << "\n"; return false; }

    std::vector<std::string> hdr;
    {
        std::stringstream ss(line);
        for (std::string tok; std::getline(ss, tok, ','); ) {
            auto a = tok.find_first_not_of(" \t\r\n");
            auto b = tok.find_last_not_of(" \t\r\n");
            hdr.push_back(a==std::string::npos ? "" : tok.substr(a, b-a+1));
        }
    }

    auto idx_of = [&](const std::string& name)->int {
        for (size_t i=0;i<hdr.size();++i) if (hdr[i]==name) return (int)i;
        return -1;
    };

    const int col_run      = idx_of("Run Number");
    const int col_un12rde  = idx_of("Un12 RDE");
    const int col_un34rde  = idx_of("Un34 RDE");
    const int col_un1112rde= idx_of("Un1112 RDE");
    const int col_un1314rde= idx_of("Un1314 RDE");
    const int col_RHDD     = idx_of("RHDD");

    if (col_run < 0) { std::cerr << "[RDE] 'Run Number' column not found\n"; return false; }

    // read rows
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> cols;
        for (std::string tok; std::getline(ss, tok, ','); ) {
            auto a = tok.find_first_not_of(" \t\r\n");
            auto b = tok.find_last_not_of(" \t\r\n");
            cols.push_back(a==std::string::npos ? "" : tok.substr(a, b-a+1));
        }
        if ((int)cols.size() <= col_run) continue;

        int run;
        try { run = std::stoi(cols[col_run]); } catch(...) { continue; }
        auto& m = rde_by_run[run];

        auto fetch = [&](int idx)->double{
            if (idx < 0 || (int)cols.size() <= idx || cols[idx].empty()) return 0.0;
            try { return std::stod(cols[idx]); } catch(...) { return 0.0; }
        };

        if (col_un12rde   >= 0) m["12"]   = fetch(col_un12rde);
        if (col_un34rde   >= 0) m["34"]   = fetch(col_un34rde);
        if (col_un1112rde >= 0) m["56"] = fetch(col_un1112rde);
        if (col_un1314rde >= 0) m["78"] = fetch(col_un1314rde);
        if (col_RHDD      >= 0) m["RHDD"] = fetch(col_RHDD);
    }
    return true;
}

void plot_event_count_corr_vs_RDE(const std::vector<PulseRow>& all_pulses,
                                  const std::string& runinfo_csv_path,
                                  const std::string& out_png_prefix)
{
    std::map<int, std::map<std::string,double>> RDE;
    if (!load_RDE_table(runinfo_csv_path, RDE)) return;

    std::vector<std::string> segs = {"12","34","56","78"};
    int colors[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};

    // ------------------------------------------------
    // Count my events: run -> seg -> count
    // ------------------------------------------------
    std::map<int,std::map<std::string,int>> my_counts;
    for (const auto& r : all_pulses) {
        if (!std::get<6>(r)) continue;           // Event flag
        int run = std::get<0>(r);        // run number
        std::string seg = std::get<1>(r);        // segment
        my_counts[run][seg] += 1;
    }

    // -----------------------
    // 1) Per-segment scatter
    // -----------------------
    TCanvas* c1 = new TCanvas("c_my_vs_rde_byseg", "My vs RDE (per segment)", 1000, 700);
    c1->SetGrid();

    TLegend leg1(0.62, 0.16, 0.88, 0.38); leg1.SetBorderSize(0); leg1.SetFillStyle(0);

    std::map<std::string,TGraph*> g_by_seg;
    for (size_t i=0;i<segs.size();++i) {
        g_by_seg[segs[i]] = new TGraph();
        g_by_seg[segs[i]]->SetMarkerStyle(20);
        g_by_seg[segs[i]]->SetMarkerSize(0.8);
        g_by_seg[segs[i]]->SetMarkerColor(colors[i]);
    }

    std::vector<double> xs_all, ys_all;

    for (const auto& kv : my_counts) {
        int run = kv.first;
        auto itR = RDE.find(run);
        if (itR == RDE.end()) continue;

        for (const auto& s : segs) {
            double x = kv.second.count(s) ? kv.second.at(s) : 0;
            double y = (itR->second.count(s) ? itR->second.at(s) : 0.0);
            int idx = g_by_seg[s]->GetN();
            g_by_seg[s]->SetPoint(idx, x, y);

            xs_all.push_back(x);
            ys_all.push_back(y);
        }
    }

    double xmin=0, ymin=0, xmax=1, ymax=1;
    for (auto& s : segs) {
        if (g_by_seg[s]->GetN()==0) continue;
        double gx0,gy0,gx1,gy1;
        g_by_seg[s]->ComputeRange(gx0,gy0,gx1,gy1);
        xmin = std::min(xmin,gx0);
        ymin = std::min(ymin,gy0);
        xmax = std::max(xmax,gx1);
        ymax = std::max(ymax,gy1);
    }
    if (xmax<=0) xmax=1;
    if (ymax<=0) ymax=1;

    const double r_all = pearson_r_vec(xs_all, ys_all);

    TH2D* fr1 = new TH2D("fr_my_vs_rde_byseg",
                        Form("Per-segment: My EVENTs vs RDE (r_all=%s);#my events (seg);#RDE events (seg)",
                            fmt_r(r_all).c_str()),
                        100, 0, xmax*1.1+1,
                        100, 0, ymax*1.1+1);
    fr1->SetStats(0); fr1->Draw();

    for (size_t i=0;i<segs.size();++i) {
        auto* g = g_by_seg[segs[i]];
        g->Draw("P SAME");
        const double r_seg = pearson_r(g);
        leg1.AddEntry(g, Form("seg %s (r=%s)", segs[i].c_str(), fmt_r(r_seg).c_str()), "p");
    }

    double dmax = std::max(xmax,ymax);
    TLine diag1(0,0,dmax,dmax);
    diag1.SetLineStyle(2); diag1.SetLineColor(kGray+2); diag1.Draw("SAME");

    leg1.Draw();
    c1->SaveAs((out_png_prefix + "_myEvent_vs_RDE_bySeg.png").c_str());
    delete c1;

    // -----------------------
    // 2) Overall scatter
    // -----------------------
    TGraph* g_tot = new TGraph();
    g_tot->SetMarkerStyle(21);
    g_tot->SetMarkerSize(1.0);
    g_tot->SetMarkerColor(kBlue+1);

    for (const auto& kv : my_counts) {
        int run = kv.first;
        auto itR = RDE.find(run);
        if (itR == RDE.end()) continue;

        int my_sum = 0; for (auto& kv2 : kv.second) my_sum += kv2.second;
        double rde_sum = 0.0; 
        for (auto& kv2 : itR->second) {
            if (kv2.first != "RHDD") rde_sum += kv2.second;
        }

        int idx = g_tot->GetN();
        g_tot->SetPoint(idx, my_sum, rde_sum);
    }

    double x0,y0,x1,y1;
    g_tot->ComputeRange(x0,y0,x1,y1);

    const double r_tot = pearson_r(g_tot);

    TCanvas* c2 = new TCanvas("c_my_vs_rde_total", "My vs RDE (total)", 900, 700);
    c2->SetGrid();
    TH2D* fr2 = new TH2D("fr_my_vs_rde_total",
                        Form("Total: My EVENTs vs RDE (r=%s);#my events (total);#RDE events (total)",
                            fmt_r(r_tot).c_str()),
                        1, 0, x1*1.1+1, 1, 0, y1*1.1+1);
    fr2->SetStats(0); fr2->Draw();


    g_tot->Draw("P SAME");
    TLine diag2(0,0,std::max(x1,y1),std::max(x1,y1));
    diag2.SetLineStyle(2); diag2.SetLineColor(kGray+2); diag2.Draw("SAME");

    c2->SaveAs((out_png_prefix + "_myEvent_vs_RDE_total.png").c_str());
    delete c2;
}

static std::vector<double> get_pe_times_abs(const std::vector<EventList>& run_data, std::string seg, double t0, double t1) {
    int idx = -1;
    if      (seg == "12") idx = 0;
    else if (seg == "34") idx = 1;
    else if (seg == "56") idx = 2;
    else if (seg == "78") idx = 3;

    std::vector<double> peTimes;
    if (idx < 0 || idx >= (int)run_data.size()) return peTimes;

    auto& events = run_data[idx];
    for (const auto& e : events) {
        if (e.realtime >= t0 && e.realtime < t1) 
            peTimes.push_back(e.realtime);
    }
    return peTimes;
}

void plot_disagreeing_coinc_pulses(const std::vector<PulseRow>& pulses,
                                   const std::vector<CoincRow>& coincs,
                                   const std::string& out_png_prefix,
                                   const std::vector<EventList>& run_data,
                                   double tol_us, int N_show, double span_us, double bin_us)
{
    const double tol = tol_us * 1e-6; // seconds

    using RunSeg = std::pair<int, std::string>;

    std::map<RunSeg, std::vector<double>> my_times_by_rs;
    my_times_by_rs.clear();
    for (const auto& r : pulses) {
        if (!std::get<6>(r)) continue; // Event only
        int run_id = std::get<0>(r);
        const std::string& s = std::get<1>(r);
        double t = std::get<2>(r);
        my_times_by_rs[{run_id, s}].push_back(t);
    }
    for (auto& kv : my_times_by_rs) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end());
    }

    struct Disagree { int run; std::string seg; double t; double abs_dt; };
    std::vector<Disagree> disagree_all;
    std::map<std::string, std::vector<double>> disagree_t_by_seg;

    auto nearest_abs_dt = [](const std::vector<double>& v, double t)->double {
        if (v.empty()) return std::numeric_limits<double>::infinity();
        auto it = std::lower_bound(v.begin(), v.end(), t);
        double dt = std::numeric_limits<double>::infinity();
        if (it != v.end())  dt = std::min(dt, std::abs(*it - t));
        if (it != v.begin()) dt = std::min(dt, std::abs(*(it-1) - t));
        return dt;
    };

    for (const auto& r : coincs) {
        if (!std::get<4>(r)) continue; // Event only
        int run_id = std::get<0>(r);
        const std::string& s = std::get<1>(r);
        double t = std::get<2>(r);

        RunSeg rs{run_id, s};
        auto it = my_times_by_rs.find(rs);
        double dt = (it == my_times_by_rs.end()) ? std::numeric_limits<double>::infinity()
                                                 : nearest_abs_dt(it->second, t);
        if (dt > tol) {
            disagree_all.push_back({run_id, s, t, dt});
            disagree_t_by_seg[s].push_back(t);
        }
    }

    if (disagree_all.empty()) {
        std::cerr << "[DISAGREE/PE] no disagreeing pulses.\n";
        return;
    }

    int M = std::min<int>((int)disagree_all.size(), N_show);

    // panel layout: up to 10 (2x5)
    int nx = 2, ny = (M <= 2 ? 1 : (M <= 4 ? 2 : (M <= 6 ? 3 : (M <= 8 ? 4 : 5))));
    if (M <= 1) nx = 1;

    const double span_s = span_us * 1e-6;
    const int nbins = std::max(1, (int)std::ceil(span_us / bin_us));

    // --- multi-pad canvas with up to N_show tiles ---------------------------
    {
        TCanvas* c = new TCanvas("c_disagree_firstN_pe", "Disagreeing coinc: PE locations", 4800, 400*ny);
        c->Divide(nx, ny);

        for (int i = 0; i < M; ++i) {
            const auto& d = disagree_all[i];
            double t0 = d.t - 50e-6; // anchor at the coinc time - 50us
            double t1 = t0 + span_s; // plot forward window; change if you want symmetric

            auto pe_abs = get_pe_times_abs(run_data, d.seg, t0, t1);
            // std::vector<double> pe_rel_us; pe_rel_us.reserve(pe_abs.size());
            // for (double t : pe_abs) pe_rel_us.push_back(t);

            c->cd(i+1); gPad->SetGrid();

            TH1I* h = new TH1I(Form("h_disagree_pe_%d", i),
                               Form("Run %d seg %s  t_{0}=%.7f s  |#Delta t|=%.3g s; t_{rel} [#mu s]; Counts per bin",
                                    d.run, d.seg.c_str(), d.t, d.abs_dt),
                               nbins, t0, t1);
            for (double u : pe_abs) {
                if (u >= t0 && u < t1) h->Fill(u);
            }
            h->SetLineColor(kBlack);
            h->SetLineWidth(2);
            h->Draw("HIST");

            double ymax = h->GetMaximum() * 1.1;
            TLine* l20 = new TLine((d.t), 0, (d.t), ymax);
            l20->SetLineColor(kRed);
            l20->SetLineStyle(2); // dashed
            l20->SetLineWidth(2);
            l20->Draw("SAME");
        }

        c->SaveAs((out_png_prefix + "_disagree_firstN_pe.png").c_str());
        delete c;
    }
}

static std::map<int, std::map<std::string,int>>
count_pulses_with_threshold(const std::vector<PulseRow>& pulses,
                            double pe_threshold)
{
    std::map<int, std::map<std::string,int>> counts;

    bool event = false;
    for (const auto& r : pulses) {
        event = std::get<6>(r); // Event only
        const double pe = std::get<3>(r);
        if (pe < pe_threshold) continue;

        const int run = std::get<0>(r);
        const std::string& seg = std::get<1>(r);
        if (event) ++counts[run][seg];
        else --counts[run][seg]; // only valid because counting window = bg window
    }
    return counts;
}

static std::map<int, std::map<std::string,double>>
normalize_by_rhdd(const std::map<int, std::map<std::string,int>>& counts,
                  const std::map<int, std::map<std::string,double>>& RDE)
{
    std::map<int, std::map<std::string,double>> norm;

    for (const auto& [run, seg_counts] : counts) {
        double rhdd = 0.0;

        if (auto itRun = RDE.find(run); itRun != RDE.end()) {
            if (auto itRHDD = itRun->second.find("RHDD"); itRHDD != itRun->second.end()) {
                rhdd = itRHDD->second;
            }
        }

        if (rhdd <= 0.0) {
            std::cerr << "[plot_lifetime] Warning: RHDD missing or <= 0 for run "
                      << run << ". Normalized counts set to 0 (skipped would also be reasonable).\n";
        }

        for (const auto& [seg, cnt] : seg_counts) {
            norm[run][seg] = (rhdd > 0.0) ? static_cast<double>(cnt) / rhdd : 0.0;
        }
    }
    return norm;
}

void plot_lifetime(const std::vector<PulseRow>& all_pulses,
                                  const std::string& runinfo_csv_path,
                                  const std::string& out_png_prefix,
                                  const json& params)
{
    std::map<int, std::map<std::string,double>> RDE;
    if (!load_RDE_table(runinfo_csv_path, RDE)) return;

    std::vector<std::string> segs = {"12","34","56","78"};
    int colors[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};
    const int color_overall = kBlack;

    const std::string csv_counts_path = out_png_prefix + "_norm_counts_by_run_seg.csv";
    std::ofstream csv_counts(csv_counts_path);
    csv_counts << "PE,run,seg,count,RHDD,norm_count,hold_time\n";

    std::map<int, std::vector<std::pair<double, double>>> tau_all; // pe -> (tau, dtau) pairs (12, 34, 56, 78, ovr)
    for(int pe_thres = 5; pe_thres < 50; ++pe_thres) {
        auto my_counts = count_pulses_with_threshold(all_pulses, (double)pe_thres);
        auto my_counts_norm = normalize_by_rhdd(my_counts, RDE);
        
        for (const auto& [run, segmap] : my_counts_norm) {
            double ht = std::numeric_limits<double>::quiet_NaN();
            lifefit::get_hold_time(params, run, ht); // ok if missing; stays NaN

            double rdhh = std::numeric_limits<double>::quiet_NaN();
            if (auto rit = RDE.find(run); rit != RDE.end()) {
                auto it = rit->second.find("RHDD");
                rdhh = it->second;
            }

            auto run_raw_it = my_counts.find(run);
            const auto* segmap_raw = (run_raw_it != my_counts.end()) ? &run_raw_it->second : nullptr;

            for (const auto& seg : segs) {
                auto it = segmap.find(seg);
                if (it == segmap.end()) continue;
                const double y = it->second; // RHDD-normalized net counts (event-bg)
                if (!std::isfinite(y)) continue;
                double count_raw = std::numeric_limits<double>::quiet_NaN();
                if (segmap_raw) {
                    if (auto itr = segmap_raw->find(seg); itr != segmap_raw->end()) count_raw = itr->second;
                }
                csv_counts << pe_thres << ','
                           << run      << ','
                           << seg      << ','
                           << count_raw<< ','
                           << rdhh     << ','
                           << y        << ','
                           << ht       << '\n';
            }
        }

        auto fit_12 = lifefit::fit_tau(my_counts_norm, params, lifefit::DaggerMode::Segment, "12");
        auto fit_34 = lifefit::fit_tau(my_counts_norm, params, lifefit::DaggerMode::Segment, "34");
        auto fit_56 = lifefit::fit_tau(my_counts_norm, params, lifefit::DaggerMode::Segment, "56");
        auto fit_78 = lifefit::fit_tau(my_counts_norm, params, lifefit::DaggerMode::Segment, "78");
        auto fit_all = lifefit::fit_tau(my_counts_norm, params, lifefit::DaggerMode::Overall);

        std::vector<std::pair<double,double>> row;
        row.reserve(5);
        row.emplace_back(fit_12.tau,  fit_12.dtau);
        row.emplace_back(fit_34.tau,  fit_34.dtau);
        row.emplace_back(fit_56.tau,  fit_56.dtau);
        row.emplace_back(fit_78.tau,  fit_78.dtau);
        row.emplace_back(fit_all.tau, fit_all.dtau);
        
        tau_all[pe_thres] = std::move(row);
    }

    csv_counts.close();
    std::cout << "[plot_lifetime] Wrote normalized counts CSV: " << csv_counts_path << "\n";

    // ---- prep data for graphs & CSV ----
    std::vector<double> x12, y12, e12;
    std::vector<double> x34, y34, e34;
    std::vector<double> x56, y56, e56;
    std::vector<double> x78, y78, e78;
    std::vector<double> xAll, yAll, eAll;

    for (const auto& [pe, row] : tau_all) {
        // row indices: 0->12, 1->34, 2->56, 3->78, 4->Overall
        if (row.size() != 5) continue;

        const auto& [t12, dt12] = row[0];
        const auto& [t34, dt34] = row[1];
        const auto& [t56, dt56] = row[2];
        const auto& [t78, dt78] = row[3];
        const auto& [ta,  dta ] = row[4];

        if (std::isfinite(t12) && std::isfinite(dt12)) { x12.push_back(pe); y12.push_back(t12); e12.push_back(dt12); }
        if (std::isfinite(t34) && std::isfinite(dt34)) { x34.push_back(pe); y34.push_back(t34); e34.push_back(dt34); }
        if (std::isfinite(t56) && std::isfinite(dt56)) { x56.push_back(pe); y56.push_back(t56); e56.push_back(dt56); }
        if (std::isfinite(t78) && std::isfinite(dt78)) { x78.push_back(pe); y78.push_back(t78); e78.push_back(dt78); }
        if (std::isfinite(ta ) && std::isfinite(dta )) { xAll.push_back(pe); yAll.push_back(ta ); eAll.push_back(dta ); }
    }

    // ---- write CSV ----
    {
        const std::string csv_path = out_png_prefix + "_tau_vs_pe.csv";
        std::ofstream csv(csv_path);
        csv << "PE,"
               "tau12,dtau12,"
               "tau34,dtau34,"
               "tau56,dtau56,"
               "tau78,dtau78,"
               "tauAll,dtauAll\n";

        for (const auto& [pe, row] : tau_all) {
            const auto val = [&](size_t i)->std::pair<double,double>{
                if (i < row.size()) return row[i];
                return {std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN()};
            };
            auto [t12, dt12] = val(0);
            auto [t34, dt34] = val(1);
            auto [t56, dt56] = val(2);
            auto [t78, dt78] = val(3);
            auto [ta,  dta ] = val(4);

            csv << pe << ","
                << t12 << "," << dt12 << ","
                << t34 << "," << dt34 << ","
                << t56 << "," << dt56 << ","
                << t78 << "," << dt78 << ","
                << ta  << "," << dta  << "\n";
        }
        csv.close();
        std::cout << "[plot_lifetime] Wrote CSV: " << csv_path << "\n";
    }

    // ---- build ROOT graphs ----
    auto make_graph = [](const std::vector<double>& xs,
                         const std::vector<double>& ys,
                         const std::vector<double>& es,
                         Color_t color,
                         Style_t mstyle,
                         const char* title)->TGraphErrors*
    {
        const int n = static_cast<int>(xs.size());
        if (n == 0) return nullptr;

        // TGraphErrors wants raw pointers; vectors are contiguous
        auto* g = new TGraphErrors(n,
                                   const_cast<double*>(xs.data()),
                                   const_cast<double*>(ys.data()),
                                   nullptr, // ex = nullptr -> treated as 0
                                   const_cast<double*>(es.data()));
        g->SetTitle(title);
        g->SetMarkerStyle(mstyle);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(2);
        return g;
    };

    TGraphErrors* g12  = make_graph(x12,  y12,  e12,  colors[0], kFullCircle,  "Seg 12");
    TGraphErrors* g34  = make_graph(x34,  y34,  e34,  colors[1], kFullSquare,  "Seg 34");
    TGraphErrors* g56  = make_graph(x56,  y56,  e56,  colors[2], kFullTriangleUp, "Seg 56");
    TGraphErrors* g78  = make_graph(x78,  y78,  e78,  colors[3], kFullTriangleDown, "Seg 78");
    TGraphErrors* gAll = make_graph(xAll, yAll, eAll, color_overall, kOpenCircle, "Overall");

    // ---- draw ----
    TCanvas* c = new TCanvas("c_tau_vs_pe", "Lifetime vs PE threshold", 1000, 700);
    c->SetGrid();

    // Frame via first available graph
    TGraphErrors* gframe = gAll ? gAll : (g12 ? g12 : (g34 ? g34 : (g56 ? g56 : g78)));
    if (!gframe) {
        std::cerr << "[plot_lifetime] No valid points to draw.\n";
        return;
    }

    // Establish reasonable axis ranges
    double xmin =  1e9, xmax = -1e9, ymin =  1e9, ymax = -1e9;
    auto update_range = [&](const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& es) {
        for (size_t i=0;i<xs.size();++i) {
            xmin = std::min(xmin, xs[i]);
            xmax = std::max(xmax, xs[i]);
            ymin = std::min(ymin, ys[i] - es[i]);
            ymax = std::max(ymax, ys[i] + es[i]);
        }
    };
    if (g12)  update_range(x12,  y12,  e12);
    if (g34)  update_range(x34,  y34,  e34);
    if (g56)  update_range(x56,  y56,  e56);
    if (g78)  update_range(x78,  y78,  e78);
    if (gAll) update_range(xAll, yAll, eAll);

    // add margins
    const double dx = (xmax > xmin) ? 0.05*(xmax-xmin) : 1.0;
    const double dy = (ymax > ymin) ? 0.10*(ymax-ymin) : 1.0;

    // Draw an empty frame using a dummy graph
    auto* ax = new TH1F("ax",";PE threshold;Lifetime #tau", 100, xmin - dx, xmax + dx);
    ax->SetMinimum(ymin - dy);
    ax->SetMaximum(ymax + dy);
    ax->Draw();

    if (g12)  g12->Draw("P SAME");
    if (g34)  g34->Draw("P SAME");
    if (g56)  g56->Draw("P SAME");
    if (g78)  g78->Draw("P SAME");
    if (gAll) gAll->Draw("P SAME");

    auto* leg = new TLegend(0.15, 0.70, 0.45, 0.90);
    if (g12)  leg->AddEntry(g12,  "Seg 12", "lep");
    if (g34)  leg->AddEntry(g34,  "Seg 34", "lep");
    if (g56)  leg->AddEntry(g56,  "Seg 56", "lep");
    if (g78)  leg->AddEntry(g78,  "Seg 78", "lep");
    if (gAll) leg->AddEntry(gAll, "Overall", "lep");
    leg->SetBorderSize(0);
    leg->Draw();

    const std::string png_path = out_png_prefix + std::string("_tau_vs_pe.png");
    c->SaveAs(png_path.c_str());
    std::cout << "[plot_lifetime] Saved plots:\n  " << png_path << "\n";
}

int main(int argc, char **argv) {
	Config cfg;
	try {
		cfg = load_config(argc, argv); // parse CLI/config, load runinfo + good runs
		std::cout << "\n========== Loaded Config ===========\n" << cfg << std::endl;
		std::cout << "====================================" << std::endl;
	} catch (const std::exception& e) {
		cerr << "Error starting program: " << e.what() << endl;
		return 1;
	}

	init_global_pdf(cfg);

    std::string data_folder = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
    json epoch_json = cfg.epoch_json;
	std::string epoch = std::to_string(cfg.epoch);
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;
	vector<EventList> run_data;
	
    int year = cfg.year;
    int epoch_start = startrun;
	int epoch_end = endrun;

	if (epoch_json.contains(std::to_string(year))) {
		const auto& epoch_year = epoch_json[std::to_string(year)];
		if (epoch_year.contains(epoch)) {
			const auto& e = epoch_year[epoch];
			epoch_start = e["start_run_number"].get<int>();
			epoch_end = e["end_run_number"].get<int>();
		}
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

        run_data = processfile(data_folder, run);

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
        plot_window_sep(ws, out_prefix);
        plot_pulses_vs_coinc_by_window(pulse, coinc, ws, out_prefix);
        plot_disagreeing_coinc_pulses(pulse, coinc, out_prefix, run_data);

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
    plot_window_sep(all_windows, out_prefix);
    plot_time_ev_byNp(all_pulses, all_windows, out_prefix, cfg);
    plot_pulses_vs_coinc_by_window(all_pulses, all_coincs, all_windows, out_prefix);
    plot_event_count_corr_vs_RDE(all_pulses, cfg.coinc_results_path, out_prefix);
    plot_lifetime(all_pulses, cfg.coinc_results_path, out_prefix, params);

	return 0;
}