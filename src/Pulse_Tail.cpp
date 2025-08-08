#include "Pulse_Tail.h"
#include "File_Loader.h"
#include "Pulse_Fitting.h"
#include <fstream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include <json.hpp>

using json = nlohmann::json;

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
        if (std::get<4>(pulse)) continue;

        double pulse_time = std::get<0>(pulse); // in us

        for (const auto& e : run_data) {
            double t = e.realtime * 1e6;
            double dt = t - pulse_time;
            if (dt >= 0 && dt < maxTime) {
                int bin = static_cast<int>(dt / binWidth);
                hist[bin] += 1.0;
            }
        }
    }

    return hist;
}

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
    double binWidth = 0.1; // Assuming 0.1 us bin width

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
    if (argc < 6) {
        std::cerr << "Error: Expected at least 5 arguments.\n";
        std::cerr << "Usage: " << argv[0] << " ./data_folder/ ./output_folder/ ./runinfo.json <startRun #> <endRun #>\n";
        return 1;
    }

    std::string data_folder = ensureTrailingSlash(argv[1]);
    std::string output_folder = ensureTrailingSlash(argv[2]);
    std::string json_filename = argv[3];

    std::ifstream i(json_filename);
    json params;
    i >> params;

    int startrun = atoi(argv[4]);
    int endrun = atoi(argv[5]);

    std::vector<std::string> segment_labels = {"12", "34", "56", "78"};
    std::vector<std::vector<double>> pulse_tails(4, vector<double>(750, 0.0));

    for (int z = startrun; z < endrun; z++) {
        string run = to_string(z);
        if (params.contains(run) && params[run]["run_type"] == "production") {
            std::vector<std::vector<double>> pulse_tails_single(4, vector<double>(750, 0.0));
            vector<EventList> run_data = processfile(data_folder, run);
            if (run_data.empty()) {
                cerr << "No data found for run " << run << ". Skipping." << endl;
                continue;
            }

            double start = (double)params[run]["fill_time"] + (double)params[run]["hold_time"] + (double)params[run]["clean_time"] + 70;
            double stop = start + 60;
            double bg_start = stop + 50;

            for (size_t seg = 0; seg < run_data.size(); ++seg) {
                Pulse_Fitting fitter(run_data[seg]);
                fitter.setWindow(start * 1e6, stop * 1e6);
                fitter.setBackgroundWindow(bg_start * 1e6);
                fitter.analyze();

                auto signalPulses = fitter.getSignalPulses();
                auto tail = accumulateTailHistogram(signalPulses, run_data[seg], 0.1, 75.0); // 1us bins, 50us total
                for (size_t b = 0; b < tail.size(); ++b) {
                    pulse_tails_single[seg][b] += tail[b];
                    pulse_tails[seg][b] += tail[b];
                }
            }
            PlotTail(pulse_tails, segment_labels, output_folder + "/tail/summed_tail_response_" + run + ".csv");
            run_data.clear();
            pulse_tails_single.clear();
        }
    }    
    TCanvas* c1 = new TCanvas("c1", "Tail Histograms", 1000, 700);
    c1->SetGrid();
    gStyle->SetOptStat(0);

    std::vector<TH1D*> histograms;
    std::vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta};
    TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
    double binWidth = 0.1;
    int nBins = pulse_tails[0].size();

    for (size_t seg = 0; seg < pulse_tails.size(); ++seg) {
        std::string name = "h" + segment_labels[seg];
        std::string title = "Tail Segment " + segment_labels[seg];
        TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, 0, nBins * binWidth);

        for (int i = 0; i < nBins; ++i) {
            hist->SetBinContent(i + 1, pulse_tails[seg][i]);
        }

        hist->SetLineColor(colors[seg % colors.size()]);
        hist->SetLineWidth(2);
        hist->SetTitle("Summed Tail Response;Time after pulse (#mus);Counts");
        histograms.push_back(hist);
        legend->AddEntry(hist, ("Segment " + segment_labels[seg]).c_str(), "l");

        if (seg == 0)
            hist->Draw("HIST");
        else
            hist->Draw("HIST SAME");
    }

    legend->Draw();
    c1->SaveAs((output_folder + "/tail_png/summed_tail_response" + to_string(startrun) + "_" + to_string(endrun) + ".png").c_str());

    return 0;
}