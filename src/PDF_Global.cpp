#include "PDF_Global.h"
#include <cmath>

std::unique_ptr<PDF_Lookup> g_pdf_lookup;

void init_global_pdf(const Config& cfg) {
    g_pdf_lookup = std::make_unique<PDF_Lookup>(
        cfg.bin_width_us, cfg.fine_bin_width_us, cfg.pdf_csv_path
    );
}

static inline bool same_bw(double a, double b) {
    return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(a));
}

std::vector<double> shifted_pdf(int seg_id, double binWidth, int shift_bins, int nBins) {
    std::vector<double> dst(nBins, 0.0);
    if (!g_pdf_lookup) return dst;

    const auto* base = g_pdf_lookup->get(seg_id);
    if (!base) return dst;

    const std::vector<double>* src = nullptr;
    if (same_bw(binWidth, base->binWidth)) src = &base->p;
    else if (same_bw(binWidth, base->fineBinWidth)) src = &base->p_f;
    else return dst;

    for (int i = 0; i < nBins; ++i) {
        int j = i - shift_bins;
        if (j >= 0 && j < (int)src->size()) dst[i] = (*src)[j];
    }
    return dst;
}