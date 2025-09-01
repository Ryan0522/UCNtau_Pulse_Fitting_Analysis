#ifndef FILE_LOADER_H
#define FILE_LOADER_H

#include <vector>
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
    
    double bin_width_us; // default 1.0
    double fine_bin_width_us; // default 0.25
    int min_gap_us; // default 10 
    std::string pdf_csv_path; // default "./config/all_tail_response.csv"

    int seeding_window;
    double pe_min_thresh;

    json runinfo_json;
    std::set<std::string> good_runs_set;

    bool plot_fits;
    int pileup_min_pulses;

    bool plot_outliers;
    int outlier_min_obs;
    double outlier_ratio_low;

    // --- NEW knobs to change in config.json
    double shift_us; // default 5.0 (time zero offset to add at output due to PDF generation)
    double seed_pe_default; // default 20.0 (default PE for the first pulse in each valid window)
    double gradient_threshold; // default 2.0 (findGradientPeaks threshold factor)
    int guard_bin; // default 1 (exclude extra peaks near the first pulse)

    bool use_coinc;
    double coinc_win_us;
} Config;

using EventList = std::vector<event>;

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