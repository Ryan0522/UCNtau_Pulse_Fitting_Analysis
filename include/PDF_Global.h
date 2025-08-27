#ifndef PDF_GLOBAL_H
#define PDF_GLOBAL_H

#include <memory>
#include "File_Loader.h"
#include "PDF_Lookup.h"

extern std::unique_ptr<PDF_Lookup> g_pdf_lookup;

void init_global_pdf(const Config& cfg);
std::vector<double> shifted_pdf(int seg_id, double binWidth, int shift_bins, int nBins);
std::vector<double> get_full_pdf(int seg_id, double binWidth);

#endif
