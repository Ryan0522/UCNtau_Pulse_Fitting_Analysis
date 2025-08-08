#include "Pulse_Tail.h"
#include "File_Loader.h"
#include "Pulse_Fitting.h"
#include <json.hpp>
#include <fstream>
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

using json = nlohmann::json;

// Accumulate the histogram for the tail response of PMTs
std::vector<double> accumulateTailHistogram(
    const std::vector<std::tuple<double, double, int, double, bool>>& pulses,
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
            double dt = t - (pulse_time - 5.0); // allow dt >= -5us by shifting origin
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
    double binWidth = 0.1; // bin width (us) must match accumulateTailHistogram call

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
		cfg = load_config(argc, argv);
		std::cout << "====================================" << std::endl;
		std::cout << "Data folder: "   << cfg.data_folder   << "\n";
        std::cout << "Output folder: " << cfg.output_folder << "\n";
		std::cout << "Runinfo path: "  << cfg.runinfo_path  << "\n";
		std::cout << "Good runs path: "<< cfg.good_runs_path<< "\n";
        std::cout << "Start run: "     << cfg.start_run     << "\n";
        std::cout << "End run: "       << cfg.end_run       << "\n";
        std::cout << "Good runs loaded: " << cfg.good_runs_set.size() << " entries\n";
		std::cout << "====================================" << std::endl;
	} catch (const std::exception& e) {
		cerr << "Error starting program: " << e.what() << endl;
		return 1;
	}

	std::string data_folder   = ensureTrailingSlash(cfg.data_folder);
    std::string output_folder = ensureTrailingSlash(cfg.output_folder);
    int         startrun      = cfg.start_run;
    int         endrun        = cfg.end_run;
    bool        save_to_txt   = cfg.save_to_txt;
	json params = cfg.runinfo_json;
	const std::set<std::string>& good_runs = cfg.good_runs_set;

    std::vector<std::string> segment_labels = {"12", "34", "56", "78"};
    std::vector<std::vector<double>> pulse_tails(4, vector<double>(750, 0.0)); // 75us @ 0.1us/bin
    int is_valid = 0;

    for (int z = startrun; z < endrun; z++) {

        string run = std::to_string(z);
		if (good_runs.find(run) == good_runs.end()) {
			cerr << "Run " << run << " not found in good runs list. Skipping." << endl;
			continue;
		}

        if (params.contains(run) && params[run]["run_type"] == "production") {
            std::vector<std::vector<double>> pulse_tails_single(4, vector<double>(750, 0.0)); // per-run accumulation
            vector<EventList> run_data = processfile(data_folder, run);
            if (run_data.empty()) {
                cerr << "No data found for run " << run << ". Skipping." << endl;
                continue;
            }

            double start = (double)params[run]["fill_time"] + (double)params[run]["hold_time"] + (double)params[run]["clean_time"] + 70;
            double stop = start + 60;
            double bg_start = stop + 50;

            for (size_t seg = 0; seg < run_data.size(); ++seg) {
                // fit pulses on this segment to identify neutron events
                Pulse_Fitting fitter(run_data[seg]);
                fitter.setWindow(start * 1e6, stop * 1e6);
                fitter.setBackgroundWindow(bg_start * 1e6);
                fitter.analyze();

                auto signalPulses = fitter.getSignalPulses();
                auto tail = accumulateTailHistogram(signalPulses, run_data[seg], 0.1, 75.0); // 0.1us bins, 75us range
                for (size_t b = 0; b < tail.size(); ++b) {
                    pulse_tails_single[seg][b] += tail[b];
                    pulse_tails[seg][b] += tail[b];
                }
            }
            PlotTail(pulse_tails_single, segment_labels, output_folder + "/tail/summed_tail_response_" + run + ".csv");
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
    double binWidth = 0.1;
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

    c1->SaveAs((output_folder + "/graphs/summed_tail_response" +
                std::to_string(startrun) + "_" + std::to_string(endrun) + ".png").c_str());
    return 0;
}