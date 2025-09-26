#ifndef RESULTS_COMP_H
#define RESULTS_COMP_H

#include <string>
#include <json.hpp>
#include <tuple>
#include <vector>
#include "File_Loader.h" // for EventList

using json = nlohmann::json;

// PulseAnalysis_<run>.csv rows: Segment(string), Time(us)(double), PE(double)
using PulseRow = std::tuple<int, std::string, double, double, double, bool, bool, int>; // (run, seg, time_us, pe, window_width, is_fineBin, is_Event, holding time)
// Coincidence rows: segment(string), time(s)(double), N(int)
using CoincRow = std::tuple<int, std::string, double, int, bool, int>;    // (run, seg, time_s, pe, is_event, holding time)

struct WindowRow {
    int run;
    std::string segment;
    int windowIndex;
    double start;
    double end;
    double binWidth_us;
    int nPulses;
    double negLogL;
    int N_obs;
    double N_exp;
};

std::vector<PulseRow> load_pulse_results(int run_number, std::string epoch, const std::string& output_folder, const int hold_t);
std::vector<CoincRow> load_coinc_results(int run_number, std::string epoch, const std::string& output_folder, const int hold_t);
std::vector<WindowRow> load_window_stats(int run_number, std::string epoch, const std::string& output_folder);
void plot_neglog_hist(const std::vector<WindowRow>& ws_all ,const std::string& out_png_prefix, int nbins = 100);
void plot_obs_exp_corr(const std::vector<WindowRow>& ws_all, const std::string& out_png_prefix);

void plot_pulses_vs_coinc_by_window(const std::vector<PulseRow>& pulses,
                                    const std::vector<CoincRow>& coincs,
                                    const std::vector<WindowRow>& windows,
                                    const std::string& out_png_prefix);

// void plot_event_count_corr_vs_RDE(const std::vector<PulseRow>& all_pulses,
//                                   const std::string& runinfo_txt_path,
//                                   const std::string& out_png_prefix);

#endif // RESULTS_COMP_H