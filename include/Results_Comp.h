#ifndef RESULTS_COMP_H
#define RESULTS_COMP_H

#include <string>
#include <json.hpp>
#include <tuple>
#include <vector>
#include "File_Loader.h" // for EventList

using json = nlohmann::json;

using PulseRow = std::tuple<std::string, double, double>; // PE can be fractional
using CoincRow = std::tuple<std::string, double, int>; // PE count is integer

std::vector<PulseRow> load_pulse_results(int run_number, const )
std::vector<PulseRow> load_coinc_results(int run_number, const )

#endif // RESULTS_COMP_H