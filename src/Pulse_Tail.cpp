#include "Pulse_Tail.h"
#include "File_Loader.h"
#include "Pulse_Fitting.h"
#include "PDF_Global.h"
#include <json.hpp>
#include <fstream>
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <filesystem>

using json = nlohmann::json;
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

// Accumulate the histogram for the tail response of PMTs
std::vector<double> accumulateTailHistogram(
    const std::vector<std::tuple<double, double, int, double, bool, bool>>& pulses,
    const EventList& run_data,
    double binWidth, 
    double maxTime)
{
    int nBins = static_cast<int>(std::ceil(maxTime / binWidth));
    std::vector<double> xCenters(nBins), hist(nBins, 0.0);

    for (int i = 0; i < nBins; ++i) {
        xCenters[i] = i * binWidth;
    }

    for (const auto& pulse : pulses) {
        if (std::get<4>(pulse)) continue; // skip windows with pileup

        double pulse_time = std::get<0>(pulse); // pulse time (us)

        for (const auto& e : run_data) {
            double t = e.realtime * 1e6; // PE time (us)
            double dt = t - (pulse_time - 10.0); // allow dt >= -5us by shifting origin
            if (dt >= 0 && dt < maxTime) {
                int bin = static_cast<int>(dt / binWidth);
                hist[bin] += 1.0;
            }
        }
    }

    return hist;
}

// Save the tail as a CSV file
void PlotTail(const std::vector<std::vector<double>>& tails, const std::vector<std::string>& segment_labels, const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << output_path << std::endl;
        return;
    }

    out << "Time(us)";
    for (const auto& label : segment_labels) {
        out << ",Segment_" << label;
    }
    out << "\n";

    int nBins = tails[0].size();
    double binWidth = 0.01; // bin width (us) must match accumulateTailHistogram call

    for (int i = 0; i < nBins; ++i) {
        out << i * binWidth;
        for (const auto& segment : tails) {
            out << "," << segment[i];
        }
        out << "\n";
    }

    out.close();
}

using namespace std;

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

    // disable plotting to save time
    cfg.plot_fits = false;
    cfg.plot_outliers = false;

	std::string data_folder   = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
    bool        save_to_txt   = cfg.save_to_txt;
    json epoch_json = cfg.epoch_json;
    std::string epoch = std::to_string(cfg.epoch);
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;

    std::vector<std::vector<double>> pulse_tails(4, vector<double>(10000, 0.0)); // 100us @ 0.01us/bin
    int is_valid = 0;
    std::vector<EventList> run_data;

    int year = cfg.year;
    int epoch_start = startrun;
	int epoch_end = endrun;

    std::vector<std::string> segment_labels;
    if (year == 2022) {
        segment_labels = {"12", "34", "56", "78"};
    } else {
        segment_labels = {"12"};
    }

	if (epoch_json.contains(std::to_string(year))) {
		const auto& epoch_year = epoch_json[std::to_string(year)];
		if (epoch_year.contains(epoch)) {
			const auto& e = epoch_year[epoch];
			epoch_start = e["start_run_number"].get<int>();
			epoch_end = e["end_run_number"].get<int>();
		}
	}

    for (int z = startrun; z < endrun; z++) {

        string run = std::to_string(z);

        if (z < epoch_start || z > epoch_end) {
			std::cerr << "Run " << run << " outside epoch " << epoch
                          << " range [" << epoch_start << ", " << epoch_end << "]. Skipping.\n";
            continue;
		}

		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

        if (params.contains(run) && params[run]["run_type"] == "production") {
            std::vector<std::vector<double>> pulse_tails_single(4, vector<double>(10000, 0.0)); // per-run accumulation
            run_data = processfile(data_folder, run);
            if (run_data.empty()) {
                cerr << "No data found for run " << run << ". Skipping." << endl;
                continue;
            }

            double start = (double)params[run]["fill_time"] + (double)params[run]["hold_time"] + (double)params[run]["clean_time"] + 70;
            double stop = start + 60;
            double bg_start = stop + 50;

            for (size_t seg = 0; seg < run_data.size(); ++seg) {
                // fit pulses on this segment to identify neutron events
                Pulse_Fitting fitter(run_data[seg], start);
                fitter.initFromConfig(cfg);
                fitter.setSegmentId(stoi(segment_labels[seg]));
                fitter.setConfigKnobs(cfg.shift_us, cfg.seed_pe_default, cfg.gradient_threshold, cfg.guard_bin);

                fitter.setWindow(start * 1e6, stop * 1e6);
                fitter.setBackgroundWindow(bg_start * 1e6);
                fitter.analyze();

                auto signalPulses = fitter.getSignalPulses();
                auto tail = accumulateTailHistogram(signalPulses, run_data[seg], 0.01, 100.0); // 0.01us bins, 100us range
                for (size_t b = 0; b < tail.size(); ++b) {
                    pulse_tails_single[seg][b] += tail[b];
                    pulse_tails[seg][b] += tail[b];
                }
            }
            fs::create_directories(output_folder + "tail/epoch_" + epoch);
            PlotTail(pulse_tails_single, segment_labels, output_folder + "tail/epoch_" + epoch + "/summed_tail_response_" + run + ".csv");
            // write per-run CSV of cumulative tails (all segments)
            run_data.clear();
            pulse_tails_single.clear();
            is_valid++;
        }
    }    

    if (is_valid == 0) return 0;

    TCanvas* c1 = new TCanvas("c1", "Tail Histograms", 1200, 600);
    c1->Divide(2, 1);

    gStyle->SetOptStat(0);
    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta};
    double binWidth = 0.01;
    int nBins = pulse_tails[0].size();

    std::vector<TH1D*> hists;
    hists.reserve(pulse_tails.size());
    for (size_t seg = 0; seg < pulse_tails.size(); ++seg) {
        std::string name  = "h" + segment_labels[seg];
        std::string title = "Segment " + segment_labels[seg];
        TH1D* h = new TH1D(name.c_str(), title.c_str(), nBins, 0, nBins * binWidth);
        for (int i = 0; i < nBins; ++i) h->SetBinContent(i+1, pulse_tails[seg][i]);
        h->SetLineColor(colors[seg % colors.size()]);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitle("Time after pulse (#mu s)");
        h->GetYaxis()->SetTitle("Counts");
        hists.push_back(h);
    }

    // ---- Pad 1: linear ----
    c1->cd(1);
    gPad->SetGrid();
    TLegend* leg1 = new TLegend(0.65, 0.70, 0.88, 0.88);
    for (size_t i = 0; i < hists.size(); ++i) {
        if (i == 0) hists[i]->SetTitle("Summed Tail Response (linear)");
        hists[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        leg1->AddEntry(hists[i], ("Segment " + segment_labels[i]).c_str(), "l");
    }
    leg1->Draw();

    // ---- Pad 2: log ----
    c1->cd(2);
    gPad->SetGrid();
    gPad->SetLogy();  // log-scale y-axis
    TLegend* leg2 = new TLegend(0.65, 0.70, 0.88, 0.88);
    for (size_t i = 0; i < hists.size(); ++i) {
        if (i == 0) hists[i]->SetTitle("Summed Tail Response (log y)");
        hists[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        leg2->AddEntry(hists[i], ("Segment " + segment_labels[i]).c_str(), "l");
    }
    leg2->Draw();

    fs::create_directories(output_folder + "graphs/epoch_" + epoch + "/PE_Response");
    c1->SaveAs((output_folder + "/graphs/epoch_" + epoch + "/PE_Response/summed_tail_response" +
                std::to_string(startrun) + "_" + std::to_string(endrun) + ".png").c_str());
    return 0;
}