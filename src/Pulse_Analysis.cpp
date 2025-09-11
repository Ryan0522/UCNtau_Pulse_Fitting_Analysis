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
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

std::ostream& operator<<(std::ostream& os, const Config& c) {
    os << "data_folder        = " << c.data_folder << "\n";
    os << "output_folder      = " << c.output_folder << "\n";
    os << "runinfo_path       = " << c.runinfo_path << "\n";
    os << "good_runs_path     = " << c.good_runs_path << "\n";
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
	string cname = "c_fit_" + out_png;
    TCanvas* c = new TCanvas(cname.c_str(), "Fit window", 1000, 600);
	c->SetMargin(0.12, 0.05, 0.12, 0.08);
	c->SetGrid(1, 1);
	
	// Histogram (observed) with variable bin edges from x (left edges)
	const int nb = (int)hist.size();
	const double dx = (x.size() > 1 ? (x[1] - x[0]) : 1.0);

    // Histogram (observed)
	string hname = "h_obs_" + out_png;
    TH1D* h = new TH1D(hname.c_str(), title.c_str(), nb, 0.0, x.back() + dx);
    h->SetDirectory(nullptr); // don't register to gDirectory (August 20, 2025)
	for (int i = 0; i < nb; ++i) h->SetBinContent((int)i+1, hist[i]);
    h->GetXaxis()->SetTitle("Time offset in window (#mu s)");
    h->GetYaxis()->SetTitle("Counts");
	h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->Draw("HIST");

	auto make_model_hist = [&](const string& name,
							   const vector<double>& yvals,
							   Color_t lineColor,
							   Style_t lineStyle = 1,
							   Width_t lineWidth = 2) {
		TH1D* hm = new TH1D(name.c_str(), "", nb, 0.0, x.back() + dx);
        hm->SetDirectory(nullptr);
        for (int i = 0; i < nb; ++i) hm->SetBinContent(i+1, yvals[i]);
        hm->SetFillStyle(0);         // no fill — outline only (bar-like top)
        hm->SetLineColor(lineColor);
        hm->SetLineStyle(lineStyle); // 1=solid, 2=dashed
        hm->SetLineWidth(lineWidth);
        hm->Draw("HIST SAME");       // flat-top bars
        return hm;
	};

    // Make observed histogram explicitly outline-only too
    h->SetFillStyle(0);              // no fill for the observed bars
    h->SetMarkerStyle(0);
    h->Draw("HIST");                 // already created above

    // Total model as bar-like (flat-top) curve
    std::string htot_name = "h_tot_" + out_png;
    TH1D* htot = make_model_hist(htot_name, total, kBlue+1, 1, 3);

    // Components as bar-like curves (dashed)
    std::vector<TH1D*> hcomps;
    hcomps.reserve(comps.size());
    int colors[6] = {kRed+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet};
    for (size_t k = 0; k < comps.size(); ++k) {
        std::string hk_name = "h_comp_" + std::to_string(k) + "_" + out_png;
        TH1D* hk = make_model_hist(hk_name, comps[k], colors[k % 6], 2, 2);
        hcomps.push_back(hk);
    }

    // Legend (same entries, now referring to hist-style lines)
    TLegend* leg = new TLegend(0.62, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h,    "Observed",     "l");
    leg->AddEntry(htot, "Model total",  "l");
    for (size_t k = 0; k < hcomps.size(); ++k) {
        leg->AddEntry(hcomps[k], Form("Pulse %d", (int)k+1), "l");
    }
    leg->Draw();

    c->SaveAs(out_png.c_str());

    delete htot;
	for (auto* hk : hcomps) delete hk;
	delete h;
	delete leg;
	delete c;
}

static void plot_outlier_overlay(int run_number, const OutlierRecord& r, const string& out_dir, const string seg) {
	const int N = (int)r.hist.size();
	const double bw = r.binWidthUs;

	TH1D hObs("hObs", "", N, 0.0, N*bw);
	TH1D hExp("hExp", "", N, 0.0, N*bw);

	for (int b=0; b<N; ++b) {
		hObs.SetBinContent(b+1, r.hist[b]);
		hExp.SetBinContent(b+1, r.totalExpected[b]);
	}

	hObs.SetLineColor(kBlack);
    hObs.SetMarkerStyle(20);
    hObs.SetMarkerSize(0.7);

    hExp.SetLineColor(kRed+1);
    hExp.SetLineStyle(2);
    hExp.SetLineWidth(2);
    hExp.SetFillStyle(0);

	hObs.SetStats(0);
    hExp.SetStats(0);

    TCanvas c("c","c",800,500);
	c.SetGrid();
	
    hObs.GetXaxis()->SetTitle("t_{rel} [#mu s]");
    hObs.GetYaxis()->SetTitle("Counts per bin");
    hObs.SetTitle(Form("Outlier window %d  (obs=%d PDF=%.1f ratio=%.2f start=%.1f us width=%.1f us)",
                       r.windowIndex, r.nObserved, r.nExpected, r.ratioExpOverObs, r.startTimeUs));
    
	double ymax = std::max(hObs.GetMaximum(), hExp.GetMaximum());
	hObs.SetMaximum(std::max(1.0, 1.15*ymax));
	
	hObs.Draw("HIST");
    hExp.Draw("HIST SAME");

    TLegend L(0.6,0.75,0.88,0.88);
    L.SetFillStyle(0);
    L.SetBorderSize(0);
    L.AddEntry(&hObs,"Observed","l");
    L.AddEntry(&hExp,"PDF","l");
    L.Draw();

    std::string base = out_dir + std::to_string(run_number) + "_segment_" + seg + "_outlier_" + std::to_string(r.windowIndex);
    c.SaveAs((base+".png").c_str());
}

// Set up and run the analysis, output to csv file
void analysis_setup(const vector<EventList> run_data, json params, string output_folder, const Config& cfg) { // Event format: <time (us), PE #, event #, window width, # of events in window>
	
	//define signal and background windows (us) from run parameters
	double start = (double)params["fill_time"] + (double)params["hold_time"] + (double)params["clean_time"] + 40;
	double stop = start + 60;
	double bg_start = stop + 50;
	
	vector<string> segment_labels = {"12", "34", "56", "78"};
	vector<int> seg_ids = {12, 34, 56, 78};
	string pulses_file = output_folder + "results/epoch_" + to_string(cfg.epoch) + "/PulseAnalysis_" + to_string(params["run_number"]) + ".csv";
	string windows_file = output_folder + "results/epoch_" + to_string(cfg.epoch) + "/PulseWindowStats_" + to_string(params["run_number"]) + ".csv";

	fs::create_directories(fs::path(pulses_file).parent_path());
	fs::create_directories(fs::path(windows_file).parent_path());

	ofstream out(pulses_file);
	ofstream ws(windows_file);
	if (!out.is_open() || !ws.is_open()) {
		cerr << "Error opening output file: " << pulses_file << ", or: " << windows_file << endl;
		return;
	}

	out.setf(std::ios::fixed, std::ios::floatfield);
	ws.setf(std::ios::fixed, std::ios::floatfield);

	out << std::setprecision(9);
	ws << std::setprecision(9);

	out << "Segment, Time (s), PE, Window Width, isFineBinWidth, Event\n";
	ws << "Segment,Window,Start,BinWidth_us,nPulses,neglogL,N_obs,N_model\n";

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
			std::string base = ensureTrailingSlash(output_folder) + "graphs/epoch_" + to_string(cfg.epoch) + "/PE_Fitting/";
			fs::create_directories(base);
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

		if (cfg.plot_outliers) {
			std::string base = ensureTrailingSlash(output_folder) + "graphs/epoch_" + to_string(cfg.epoch) + "/debug/";
			fs::create_directories(base);

			const auto& bad = fitter.getOutliers();
			for (const auto& r : bad) {
				plot_outlier_overlay((int)params["run_number"], r, base, segment_labels[seg]);
			}
		}

		auto signalPulses = fitter.getSignalPulses();
		auto backgroundPulses = fitter.getBackgroundPulses();
		const auto& stats = fitter.getWindowStats();
		
		// write signal pulsese (Event=1)
		for (const auto& event : signalPulses) {
			out << segment_labels[seg] << ", " // Segment
				<< get<0>(event)/1e6 << ", " // event time (us)
				<< get<1>(event) << ", " // PE count
				<< get<3>(event) << ", " // WindowWidth
				<< get<5>(event) << ", " // binWidth (1) or fineBinWidth (0)
				<< "1 \n";
		}

		// write background pulses (Event=0)
		for (const auto& event : backgroundPulses) {
			out << segment_labels[seg] << ", "
				<< get<0>(event)/1e6 << ", "
				<< get<1>(event) << ", "
				<< get<3>(event) << ", "
				<< get<5>(event) << ", "
				<< "0 \n";
		}

		// write window statistics (Event only)
		for (const auto& s : stats) {
			if (s.startTimeUs / 1e6 > stop) continue;
			ws << segment_labels[seg] << ","
			<< s.windowIndex << ","
			<< s.startTimeUs/1e6 << ","
			<< s.binWidthUs << ","
			<< s.nPulsesChosen << ","
			<< -s.logL << ","
			<< s.nObserved << ","
			<< s.nExpected << "\n";
		}
	}
	
	out.close();
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

	std::string data_folder   = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
    bool        save_to_txt   = cfg.save_to_txt;
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

	if (save_to_txt) {
		cout << "** Note: converting data to text, no analysis will be performed **" << endl;
	}

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