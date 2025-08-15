#include "PDF_Lookup.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <iostream>

using std::string;
using std::vector;

PDF_Lookup::PDF_Lookup(double coarse_bw, double fine_bw, const std::string& csv_path)
    : coarseBinWidth_(coarse_bw), fineBinWidth_(fine_bw), csv_path_(csv_path)
{
    if (!csv_path.empty()) {
        if (!load_csv(csv_path)) {
            throw std::runtime_error("PDF_Lookup: Failed to load CSV: " + csv_path);
        }
    }
}

bool PDF_Lookup::load_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;

    string line;
    if (!std::getline(in, line)) return false;
    std::stringstream hs(line);
    std::vector<string> headers;
    for (std::string h; std::getline(hs, h, ',');) headers.push_back(h);
    if (headers.empty() || headers[0] != "Time(us)") return false;

    const size_t C = headers.size();
    if (C < 2) return false;

    std::vector<std::vector<double>> cols(C);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        for (size_t c = 0; c < C; ++c) {
            std::string tok;
            if (!std::getline(ss, tok, ',')) return false;
            cols[c].push_back(std::stod(tok));
        }
    }
    if (cols[0].size() < 100) return false;

    const vector<double>& t = cols[0];
    double dt001 = (t.size() > 1) ? (t[1] - t[0]) : 0.01;
    if (std::abs(dt001 - 0.01) > 1e-9) throw std::runtime_error("CSV time step not 0.01 us.");

    for (size_t c = 1; c < C; ++c) {
        int seg_id = parse_segment_name(headers[c]);
        vector<double> p001 = cols[c];

        BasePDF b;
        b.binWidth = coarseBinWidth_;
        b.fineBinWidth = fineBinWidth_;

        downsample(p001, dt001, b.p_f, b.fineBinWidth);
        normalize_series(b.p_f, b.fineBinWidth);

        downsample(b.p_f, b.fineBinWidth, b.p, b.binWidth);
        normalize_series(b.p, b.binWidth);

        seg_[seg_id] = std::move(b);
    }
    return true;
}

int PDF_Lookup::parse_segment_name(const std::string& header) {
    // Accept "Segment_12"
    const std::string prefix = "Segment_";
    if (header.rfind(prefix, 0) != 0) throw std::runtime_error("Invalid segment header: " + header);
    std::string num_str = header.substr(prefix.size());
    if (num_str.empty() || num_str.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error("Invalid segment number in header: " + header);
    }
    return std::stoi(num_str);
}

void PDF_Lookup::normalize_series(std::vector<double>& v, double dt) {
    double area = 0.0; for (double x : v) area += x * dt;
    if (area > 0.0) {
        double inv = 1.0 / area;
        for (double& x : v) x *= inv;
    }
}

void PDF_Lookup::downsample(const std::vector<double>& src, double dt_src, std::vector<double>& dst, double dt_dst) {
    if (!PDF_Lookup::approx_int_multiple(dt_dst, dt_src)) {
        throw std::runtime_error("rebin_sum: dt_dst must be an integer multiple of dt_src.");
    }
    int R = (int)std::llround(dt_dst / dt_src);
    int N = (int)(src.size() / R);
    dst.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        int a = i * R;
        for (int k = 0; k < R; ++k) s += src[a + k];  // <- SUM, not average
        dst[i] = s;
    }
}

bool PDF_Lookup::approx_int_multiple(double value, double step, double rel_eps) {
    double q = value / step;
    double rq = std::llround(q);
    return std::abs(q - rq) <= rel_eps * std::max(1.0, std::abs(q));
}
