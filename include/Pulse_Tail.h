#ifndef PULSE_TAIL_H
#define PULSE_TAIL_H

#include <string>
#include <vector>
#include <tuple>
#include "File_Loader.h"

std::vector<double> accumulateTailHistogram(
    const std::vector<std::tuple<double, double, int, double, bool>>& pulses,
    const EventList& run_data,
    double binWidth, 
    double maxTime);

void PlotTail(const std::vector<std::vector<double>>& tails, const std::vector<std::string>& segment_labels, const std::string& output_path);

#endif // PULSE_TAIL_H
