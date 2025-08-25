#ifndef RESULTS_COMP_H
#define RESULTS_COMP_H

#include <string>
#include <json.hpp>
#include <tuple>
#include <vector>
#include "File_Loader.h" // for EventList

using json = nlohmann::json;

std::vector<std::tuple<std::string, double, int>> load_pe_counts_csv(const std::string& filepath);

#endif // RESULTS_COMP_H