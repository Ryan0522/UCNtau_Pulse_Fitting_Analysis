#ifndef RESULTS_COMP_H
#define RESULTS_COMP_H

#include <string>
#include <json.hpp>
#include <tuple>
#include <vector>
#include "File_Loader.h" // for EventList

using json = nlohmann::json;

// PulseAnalysis_<run>.csv rows: Segment(string), Time(us)(double), PE(double)
using PulseRow = std::tuple<std::string, double, double, double>; // (seg, time_us, pe, window_width)
// Coincidence rows: segment(string), time(s)(double), N(int)
using CoincRow = std::tuple<std::string, double, int>;    // (seg, time_s, pe)

struct WindowRow {
    std::string segment;
    int windowIndex;
    double start;
    double binWidth_us;
    int nPulses;
    double negLogL;
    int N_obs;
    double N_exp;
};

std::vector<PulseRow> load_pulse_results(int run_number, const std::string& output_folder);
std::vector<CoincRow> load_coinc_results(int run_number, const std::string& output_folder);
std::vector<WindowRow> load_window_stats(int run_number, const std::string& output_folder);
void plot_neglog_hist(const std::vector<WindowRow>& ws_all ,const std::string& out_png_prefix, int nbins = 100);
void plot_obs_exp_corr(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix);

#endif // RESULTS_COMP_H