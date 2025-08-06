#ifndef PULSE_ANALYSIS_H
#define PULSE_ANALYSIS_H

#include "TObject.h"
#include <string>
#include <list>
#include <json.hpp>
#include <vector>

using json = nlohmann::json;

typedef struct
{
    Int_t channel;
    Int_t edge;
    Int_t tag;
    Int_t full;
    ULong64_t time;
    Double_t realtime;
} event;

using EventList = std::list<event>;

std::string ensureTrailingSlash(const std::string& folder);

std::vector<EventList> processfile( // Not writing to txt
    std::string data_folder,
    std::string runnum
);

void processfile( // Writing to txt
    std::string data_folder,
    std::string output_folder,
    std::string runnum
);

void analysis_setup(const std::vector<EventList> run_data, json params, std::string output_folder);

#endif // PULSE_ANALYSIS_H