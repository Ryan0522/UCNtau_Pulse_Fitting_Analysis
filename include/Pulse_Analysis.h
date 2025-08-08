#ifndef PULSE_ANALYSIS_H
#define PULSE_ANALYSIS_H

#include <string>
#include <json.hpp>
#include <vector>
#include "File_Loader.h" // For EventList

using json = nlohmann::json;

void analysis_setup(const std::vector<EventList> run_data, json params, std::string output_folder);

#endif // PULSE_ANALYSIS_H