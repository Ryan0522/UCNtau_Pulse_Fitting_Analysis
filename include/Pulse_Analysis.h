#ifndef PULSE_ANALYSIS_H
#define PULSE_ANALYSIS_H

#include "TObject.h"
#include "File_Loader.h"
#include <string>
#include <list>
#include <json.hpp>
#include <vector>

using json = nlohmann::json;

void analysis_setup(const std::vector<EventList> run_data, json params, std::string output_folder);

#endif // PULSE_ANALYSIS_H