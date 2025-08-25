#ifndef RESULTS_COMP_H
#define RESULTS_COMP_H

#include <string>
#include <json.hpp>
#include <tuple>
#include <vector>
#include "File_Loader.h" // for EventList

using json = nlohmann::json;

// PulseAnalysis_<run>.csv rows: Segment(string), Time(us)(double), PE(double)
using PulseRow = std::tuple<std::string, double, double>; // (seg, time_us, pe)
// Coincidence rows: segment(string), time(s)(double), N(int)
using CoincRow = std::tuple<std::string, double, int>;    // (seg, time_s, pe)

std::vector<PulseRow> load_pulse_results(int run_number, const std::string& output_folder);
std::vector<CoincRow> load_coinc_results(int run_number, const std::string& output_folder);

#endif // RESULTS_COMP_H