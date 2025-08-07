#ifndef FILE_LOADER_H
#define FILE_LOADER_H

#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <iomanip>
#include <list>

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

#endif