#ifndef FILE_LOADER_H
#define FILE_LOADER_H

#include <list>
#include <string>
#include <set>
#include <json.hpp>
#include <Rtypes.h>

using json = nlohmann::json;

typedef struct
{
    Int_t channel;
    Int_t edge;
    Int_t tag;
    Int_t full;
    ULong64_t time; // in ticks
    Double_t realtime; // in seconds
} event;

typedef struct 
{
    std::string data_folder;
    std::string output_folder;
    std::string runinfo_path;
    std::string good_runs_path;
    int start_run;
    int end_run;
    bool save_to_txt;

    json runinfo_json;
    std::set<std::string> good_runs_set;
} Config;

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

Config load_config(int argc, char** argv, const std::string& default_cfg = "./config/default_config.json");

#endif // FILE_LOADER_H