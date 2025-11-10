#include "Pulse_Fitting.h"
#include "PDF_Global.h"
#include "FitCore.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <deque>

using namespace std;

static inline int round_to_index(double x, int max_index_inclusive) {
    long long idx = std::llround(x);
    if (idx < 0) idx = 0;
    if (idx > (long long)max_index_inclusive) idx = (long long)max_index_inclusive;
    return static_cast<int>(idx);
}

Pulse_Fitting::Pulse_Fitting(const EventList& events, double start) : events_(events) { extractTimes(events, start); } // copy event.realtime to peTimes_us_ (us)

void Pulse_Fitting::setWindow(double start_us, double stop_us) {
    // set signal window in absolute microseconds
    startAfterUs_ = start_us;
    stopAfterUs_ = stop_us;
}

void Pulse_Fitting::setBackgroundWindow(double start_us) {
    backgroundAfterUs_ = start_us;
}

void Pulse_Fitting::analyze() { // Assume 60s is the length for both the counting and the background windows
    cout << "Segment: " << segmentId_ << endl;
    cout << "Event size: " << global_peTimes_us_.size() << endl; // total PE hits loaded

    vector<double> signalTimes = applyTimeWindow(global_peTimes_us_, startAfterUs_, stopAfterUs_);
    cout << "SignalTime PE Event size: " << signalTimes.size() << "  |  ";
    
    vector<double> backgroundTimes;
    if (backgroundAfterUs_ > 0)
        backgroundTimes = applyTimeWindow(global_peTimes_us_, backgroundAfterUs_, backgroundAfterUs_ + 60e6);
    cout << "Background PE Event size: " << backgroundTimes.size() << endl;
    
    const double bg_duration_sec = 60.0;
    peBackgroundRate_ = std::max(0.0, backgroundTimes.size() / bg_duration_sec);
    fitting_bg_ = false;

    const int MAX_ITERS = 10;
    const double tolRel = 1e-2;
    bool temp_debug = debug_;
    debug_ = false;

    for (int it = 0; it < MAX_ITERS; ++it) {
        backgroundPulses_.clear();
        fitRegion(backgroundTimes, backgroundPulses_);

        double sumPE = 0.0;
        for (const auto& pr : backgroundPulses_) sumPE += std::get<1>(pr);
        
        double rate_hat = std::max(0.0, (backgroundTimes.size() - sumPE) / bg_duration_sec);
        double delta = std::fabs(rate_hat - peBackgroundRate_);
        double tol = tolRel * std::max(1.0, peBackgroundRate_);

        peBackgroundRate_ = rate_hat;
        fitting_bg_ = true;
        if (delta < tol) break; // converged
    }
    
    debug_ = temp_debug;
    signalTimes = applyTimeWindow(global_peTimes_us_, startAfterUs_, stopAfterUs_);
    fitRegion(signalTimes, signalPulses_); // parse windows, fit pulses

    cout << "SignalTime Neutron Event count: " << signalPulses_.size() << "  |  ";
    cout << "Background Neutron Event count: " << backgroundPulses_.size() << "\n" << endl;

    eventBackgroundRate_ = backgroundPulses_.size() / 60.0;
}

vector<double> Pulse_Fitting::buildCarryExpected(double winStartUs, double binWidth_us, const vector<double>& xCenters) const
{
    const int L = (int)xCenters.size();
    vector<double> carry(L, 0.0);
    if (L == 0) return carry;

    const vector<double> base = get_full_pdf(segmentId_, binWidth_us);
    if (base.empty()) return carry;
    const int fullBins = (int)base.size();
    const double pdf_len_us = fullBins * binWidth_us;

    for (const auto& pr : carryPulses_) {
        const double t0_us = pr.first; // absolute time (us)
        const double PE = pr.second;
        const double age = winStartUs - t0_us; // >= 0 means the peak of the pdf is to the left of current window

        if (age < 0.0) continue; // skip pulses in current window (if any) (this is for safety)
        if (age >= pdf_len_us) continue; // the entire PDF is to the left of current window

        const int dx = (int)std::floor(age / binWidth_us);
        const int copy_len = std::min(L, fullBins - dx);
        for (int j = 0; j < copy_len; ++j) {
            carry[j] += PE * base[j + dx]; // density -> mass distribution (Sept 1st, 2025)
        }
    }
    return carry;
}

void Pulse_Fitting::extractTimes(const EventList& events, double start) {
    // linearize event.realtime (s) -> vector of times (us)
    global_peTimes_us_.clear();
    global_peChans_.clear();
    global_peTimes_us_.reserve(events.size());
    global_peChans_.reserve(events.size());
    
    for (const auto& e : events) {
        if (e.realtime < start) continue;
        global_peTimes_us_.push_back(e.realtime * 1e6);
        global_peChans_.push_back(e.channel);
    }
}

vector<double> Pulse_Fitting::applyTimeWindow(const vector<double>& times, double start, double end) {
    peTimes_us_.clear();
    peChans_.clear();

    struct Hit {
        double t;
        int ch;
    };

    vector<Hit> hits;
    hits.reserve(times.size());
    
    for (size_t i = 0; i < times.size(); ++i) {
        double t = times[i];
        int ch = global_peChans_[i];
        if (t >= start && t < end) {
            hits.push_back({t, ch});
        }
    }

    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        return a.t < b.t;
    });

    vector<double> filtered_times;
    filtered_times.reserve(hits.size());

    for (const auto& h : hits) {
        filtered_times.push_back(h.t);
        peTimes_us_.push_back(h.t);
        peChans_.push_back(h.ch);
    }

    return filtered_times;
}

tuple<double, int, double, double> Pulse_Fitting::movingWindow(const vector<double>& times, int startIdx) {
    // grow a window starting at 'startIdx' until an inter-hit gap > minGap_us_
    int N = static_cast<int>(times.size());
    if (startIdx >= N - 1) {
        return {0.0, startIdx + 1, 0.0, 0.0};
    }

    double start = times[startIdx];
    int j = startIdx + 1;

    // debug only //

    // constexpr double TARGET_US = 0 * 423.1150133 * 1e6; // <- your focus time (abs µs on the same clock as `times`)
    constexpr double TARGET_US = 0 * 411.4623209208 * 1e6; // <- your focus time (abs µs on the same clock as `times`)
    constexpr double EPS_US = 1; // 1 us window
    
    const bool dbg_hit = std::fabs(start - TARGET_US) < EPS_US;

    auto DBG = [&](auto&&... args) {
        if (!dbg_hit) return;
        std::cerr << std::fixed << std::setprecision(13);
        (std::cerr << ... << args) << '\n';
    };

    DBG("[MW] startIdx=" , startIdx
        , "  start_us=" , start
        , "  N=" , N
        , "  coinc=" , (use_coinc_ ? "on" : "off")
        , "  coinc_win_us=" , coinc_win_us_
        , "  minGap_us=" , minGap_us_);

    bool armed = false;
    int last = startIdx;

    if (use_coinc_) {
        while (j < N) {
            double dt_seed = abs(times[j] - times[startIdx]);
            if (dt_seed > coinc_win_us_) {
                DBG("  [ARM] break: dt_seed=", dt_seed, " > ", coinc_win_us_, " at j=", j);
                break;
            }

            const int ch_i = peChans_[startIdx];
            const int ch_j = peChans_[j];

            DBG("  [ARM] j=", j, "  t[j]=", times[j],
                "  dt_seed=", dt_seed, "  ch[i]/ch[j]=", ch_i, "/", ch_j);


            if (peChans_[j] != peChans_[startIdx]) {
                armed = true;
                ++j;
                DBG("  [ARM] ARMED at j=", j-1, "  (advance to j=", j, ")");
                break;
            }
            ++j;
        }
        if (!armed) {
            const double t = times[startIdx];
            DBG("  [ARM] NOT ARMED — return empty window at t=", t);
            return {0.0, startIdx + 1, t, t};
        }
        int last = j -  1;
        while (j < N) {
            const double dt_step = times[j] - times[last];
            DBG("  [GROW] j=", j, " last=", last,
                "  t[last]=", times[last], "  t[j]=", times[j],
                "  dt_step=", dt_step);

            if (dt_step >= minGap_us_) {
                DBG("  [GROW] break: dt_step=", dt_step, " >= minGap_us_=", minGap_us_);
                break;
            }
            ++j;
        }

    } else {
        while (j < N) {
            const double gap = times[j] - times[j - 1];
            DBG("  [GROW-nc] j=", j, "  gap=", gap);
            if (gap > minGap_us_) break;
            ++j;
        }
    }

    int i0 = startIdx;
    while (i0 > 0) {
        double front_gap = times[i0] - times[i0-1];
        if (front_gap > std::max(5.0, minGap_us_ / 2)) break;
        i0--;
    }
    start = times[i0];

    // previously: double end = times[j - 1]; // August 6th 2025
    // pad by minGap_us_ so the histogram includes explicit zeros after the last hit // August 10th 2025
    const double t_last   = times[j - 1];
    const double nextHit  = (j < N ? times[j] : std::numeric_limits<double>::infinity());
    const double end_pad  = t_last + minGap_us_;
    const double cap_by_next = std::isfinite(nextHit)
        ? std::nextafter(nextHit, -std::numeric_limits<double>::infinity())
        : end_pad;
    const double end = std::min(end_pad, cap_by_next);
    const double width = end - start;
    
    DBG("  [END] start_us=", start,
        "  last_idx=", last,
        "  t_last=", t_last,
        "  nextHit=", (std::isfinite(nextHit) ? nextHit : -1.0),
        "  end_pad=", end_pad,
        "  cap_by_next=", cap_by_next,
        "  end=", end,
        "  width=", width,
        "  j(return_nextIdx)=", j);

    return {width, j, start, end};
}

bool Pulse_Fitting::makeHistogram(const vector<double>& times, int i, double binWidth_us,
                                  double& windowWidth, int& j, double& startTime, double& endTime,
                                  vector<int>& hist, vector<double>& xCenters) 
{
    // compute [startTime, endTime] window and bin hits into 'hist' with given binWidth_us
    tie(windowWidth, j, startTime, endTime) = movingWindow(times, i);
    if (windowWidth < binWidth_us) {
        // cout << "Case 1 ";
        return false;
    }

    int nBins = static_cast<int>(ceil(windowWidth / binWidth_us));
    if (nBins < 1) {
        // cout << "Case 2 "; 
        return false;
    }

    xCenters.resize(nBins);
    hist.assign(nBins, 0);
    for (int b = 0; b < nBins; ++b) {
        xCenters[b] = b * binWidth_us;
    }

    for (int k = i; k < j; ++k) { // fill counts per bin relative to startTime
        double t = times[k] - startTime;
        int bin = static_cast<int>(t / binWidth_us);
        if (bin >= 0 && bin < nBins) {
            hist[bin]++;
        }
    }

    // --- NEW: prepend 5 us worth of empty bins so fitter origin is at -5 us (due to PDF generate algorithm in Pulse_Tail.cpp) ---
    int pre_bins = (binWidth_us == binWidth_us_) ? preBins_ : finePreBins_;
    if (pre_bins > 0) {
        hist.insert(hist.begin(), pre_bins, 0);

        // rebuild xCenters to stay aligned (0, 1*binWidth_us, ...)
        vector<double>newX(hist.size());
        for (size_t b = 0; b < newX.size(); ++b) newX[b] = b * binWidth_us;
        xCenters.swap(newX);
    }
    // --- end NEW (August 20, 2025) ---

    return true;
}

void Pulse_Fitting::fitRegion(const vector<double>& data_us,
                              vector<tuple<double, double, int, double, bool, bool>>& output) 
{
    // slide over data, window by window, fit pulses per window
    int i = 0;
    int N = static_cast<int>(data_us.size());
    int windowIndex = 0;
    
    while (i < N) {
        vector<int> hist;
        vector<double> xCenters;
        double windowWidth, startTime, endTime;
        int j;
        double binWidth_us = binWidth_us_;

        if (!makeHistogram(data_us, i, binWidth_us, windowWidth, j, startTime, endTime, hist, xCenters)) {
            // cout << "unable to make histogram at: " << setprecision(13) << startTime / 1e6 << " to " << endTime / 1e6 << endl;
            // if (debug_ && startTime / 1e6 >= 415) throw std::runtime_error("No....");
            i = j;
            continue;
        }

        if (xCenters.size() < 2) {
            binWidth_us = fineBinWidth_us_;
            if (!makeHistogram(data_us, i, binWidth_us, windowWidth, j, startTime, endTime, hist, xCenters)) {
                i = j;
                // cout << "unable to make histogram at: " << setprecision(13) << startTime / 1e6 << " to " << endTime / 1e6 << endl;
                continue;
            }
        }

        const vector<double> base = get_full_pdf(segmentId_, binWidth_us);
        const double pdf_len_us = base.empty() ? 0.0 : base.size() * binWidth_us;

        if (pdf_len_us > 0.0) {
            carryPulses_.erase(
                std::remove_if(carryPulses_.begin(), carryPulses_.end(),
                                [&](const auto& pr){ return (startTime - pr.first) >= pdf_len_us; }),
                carryPulses_.end()
            );
        }

        std::vector<double> carry = buildCarryExpected(startTime, binWidth_us, xCenters); // FIX
        vector<vector<double>> pdfLookup = generatePDFLookup(xCenters); // shifted PDFs cache
        vector<double> fittedPEs, fittedDTs;

        curWindowIndex_ = windowIndex;
        curWindowStartUs_ = startTime;
        curWindowEndUs_ = endTime;
        curWindowBinWidth_us_ = binWidth_us;

        bool success = fitPulses(hist, xCenters, pdfLookup, fittedPEs, fittedDTs, binWidth_us, &carry);
        if (!success) {
            i = j;
            continue;
        }

        const bool need_single = capturePlots_ && !haveFirstSingle_ && (fittedPEs.size() == 1);
        const bool need_pile = capturePlots_ && !haveFirstPileup_ && ((int)fittedPEs.size() >= pileupMinPulses_);
        
        if (need_single || need_pile) {

            auto build_components = [&](const vector<double>& PEs,
                                        const vector<double>& DTs,
                                        const vector<vector<double>>& look) {
                vector<vector<double>> comps;
                const double bg = peBackgroundRate_ * binWidth_us * 1e-6; // pe bg rate is in unit of s
                comps.reserve(PEs.size());
                vector<double> total(hist.size(), 0.0);
                for (size_t m = 0; m < PEs.size(); ++m) {
                    int dt_idx = (int)std::llround(DTs[m]);
                    dt_idx = std::max(0, std::min(dt_idx, (int)look.size() - 1));
                    vector<double> comp(hist.size(), 0.0);
                    const auto& row = look[dt_idx];
                    for (size_t j = 0; j < hist.size(); ++j) {
                        double val = PEs[m] * row[j]; // density -> mass distribution (Sept 1, 2025)
                        comp[j] = val;
                        total[j] += val;
                        if (m == 0) total[j] += bg;
                    }
                    comps.push_back(std::move(comp));
                }
                return std::make_pair(std::move(total), std::move(comps));
            };

            auto [totalExp, comps] = build_components(fittedPEs, fittedDTs, pdfLookup);

            if (need_single) {
                single_hist_ = hist;
                single_x_ = xCenters;
                single_total_ = totalExp;
                single_components_ = comps;
                haveFirstSingle_ = true;
            }
            if (need_pile) {
                pile_hist_ = hist;
                pile_x_ = xCenters;
                pile_total_ = totalExp;
                pile_components_ = comps;
                haveFirstPileup_ = true;
            }
        }

        // --- NEW: Added pulse subtraction functionalities ---
        for (size_t k = 0; k < fittedPEs.size(); ++k) {
            int dt_bins = round_to_index(fittedDTs[k], (int)xCenters.size() - 1);
            dt_bins = std::max(0, std::min(dt_bins, (int)xCenters.size() - 1)); // FIX
            double pulse_time_us = startTime + dt_bins * binWidth_us; // +5us correction
            output.emplace_back(pulse_time_us, fittedPEs[k], windowIndex, windowWidth, fittedPEs.size() > 1, binWidth_us == fineBinWidth_us_); // store result
            // cout << (double)j/(double)N << ", " << pulse_time_us / 1e6 << ", " << fittedPEs[k] << ", " << endl;
            carryPulses_.emplace_back(pulse_time_us, fittedPEs[k]); // FIX
        }
        // --- end NEW (August 26, 2025) ---

        windowIndex++;
        i = j;
    }
}

vector<vector<double>> Pulse_Fitting::generatePDFLookup(const vector<double>& xCenters) {
    if (xCenters.size() < 2) return {};
    const int length = static_cast<int>(xCenters.size());
    const double binWidth_us = xCenters[1] - xCenters[0];

    auto key = std::make_pair(length, std::round(binWidth_us * 100.0) / 100.0); // PDF is 0.01us in resolution
    auto it = pdfCache_.find(key);
    if (it != pdfCache_.end()) return it->second;

    vector<vector<double>> lookup(length, vector<double>(length, 0.0));
    for (int dx = 0; dx < length; ++dx) {
        lookup[dx] = shifted_pdf(segmentId_, binWidth_us, dx, length);
    }

    pdfCache_[key] = lookup;
    return lookup;
}

vector<double> Pulse_Fitting::findGradientPeaks(const vector<int>& hist, double thresholdFactor) {
    // simple gradient-based seed find; thresholdFactor in units of grad "std"
    if ((int)hist.size() <= 2) {
        return {};
    }

    vector<double> grad(hist.size());
    for (size_t i = 1; i + 1 < hist.size(); ++i) {
        grad[i] = static_cast<double>(hist[i + 1] - hist[i - 1]) / 2.0;
    }

    double sumSq = 0.0;
    for (double g : grad) sumSq += g * g;
    double stdGrad = sqrt(sumSq / grad.size());
    double threshold = thresholdFactor * stdGrad;

    vector<double> peaks;
    for (size_t i = 1; i + 1 < grad.size(); ++i) {
        if (grad[i] > threshold) {
            peaks.push_back((double)i);
        }
    }
    return peaks;
}

vector<double> Pulse_Fitting::findCoincidenceSeeds(double startUs, double endUs, double binWidthUs, int pre_bins) const
{
    std::vector<double> seeds;
    const int N = (int)peTimes_us_.size();
    if (N < 2) return seeds;

    auto it0 = std::lower_bound(peTimes_us_.begin(), peTimes_us_.end(), startUs);
    int i0 = (int)std::distance(peTimes_us_.begin(), it0);
    double prev_time_us = 0.0; // us

    for (int i = i0; i < N; ++i) {
        const double t0 = peTimes_us_[i];
        if (t0 < prev_time_us) continue; // skip hits that are too close to previous seed
        if (t0 >= endUs) break;
        const int ch0 = peChans_[i];

        bool armed = false;
        int total = 1;
        for (int j = i + 1; j < N; ++j) {
            const double tj = peTimes_us_[j];
            if (tj >= endUs) break;
            const double dt = abs(tj - t0);
            if (dt > coinc_win_us_) break;
            if (!armed && peChans_[j] != ch0) armed = true;
            ++total;
        }

        const bool accepted = (armed && total >= coinc_seed_pe_min_);
        if (armed && total >= coinc_seed_pe_min_) {
            const double dt_bins = pre_bins + (t0 - startUs) / binWidthUs; // fractional
            seeds.push_back(dt_bins);
            prev_time_us = t0 + seeding_window_us_; // us
        }
    }

    std::sort(seeds.begin(), seeds.end());
    const double eps = 1e-4; // bins (min separation)
    seeds.erase(std::unique(seeds.begin(), seeds.end(),
                            [eps](double a, double b){ return std::fabs(a - b) < eps; }),
                seeds.end());

    return seeds;
}

// struct Hit { double t; int ch; };

// vector<double> Pulse_Fitting::findCoincidenceSeeds(double startUs, double endUs, double binWidthUs, int pre_bins) const {
//     std::vector<Hit> hits;
//     hits.reserve(peTimes_us_.size());
//     for (size_t k = 0; k < peTimes_us_.size(); ++k) {
//         double t = peTimes_us_[k];
//         if (t >= startUs && t < endUs) hits.push_back({t, peChans_[k]});
//     }
//     std::sort(hits.begin(), hits.end(), [](auto& a, auto& b){ return a.t < b.t; });

//     std::deque<Hit> w;
//     std::array<int, 256> chCount{};
//     int distinctCh = 0;
    
//     auto add = [&](const Hit& h) {
//         if (chCount[h.ch]++ == 0) ++ distinctCh;
//         w.push_back(h);
//     };
//     auto pop_front = [&](){
//         const Hit& h = w.front();
//         if (--chCount[h.ch] == 0) --distinctCh;
//         w.pop_front();
//     };

//     std::vector<double> seeds;
//     double next_ok_time = -1e99;

//     for (const Hit& h : hits) {
//         if (h.t < next_ok_time) continue;
//         add(h);
//         while (!w.empty() && (w.back().t - w.front().t) > coinc_win_us_) pop_front();

//         const int count = (int)w.size();
//         const bool cross = (distinctCh >= 2);
//         if (count >= coinc_seed_pe_min_ && cross) {
//             const double t_seed = w.front().t;
//             const double dt_bins = pre_bins + (t_seed - startUs) / binWidthUs;
//             seeds.push_back(dt_bins);

//             next_ok_time = t_seed + seeding_window_us_;
//             while (!w.empty() && w.front().t < next_ok_time) pop_front();
//         }
//     }

//     return seeds;
// }

bool Pulse_Fitting::fitPulses(const vector<int>& hist, const vector<double>& xCenters,
                              const vector<vector<double>>& pdfLookup,
                              vector<double>& fittedPEs, vector<double>& fittedDTs,
                              const double binWidth_us,
                              const vector<double>* fixedExpected) 
{
    const int pre_bins = static_cast<int>(std::llround(shiftUs_ / binWidth_us)); // Since PDF is shifted by +5 us
    bool dbg = (curWindowIndex_ == debug_window_index_ && 
                    segmentId_ == debug_segment_id_ && debug_);

    constexpr double TARGET_US = 1 * 415.062156 * 1e6; // <- your focus time (abs µs on the same clock as `times`)
    constexpr double EPS_US = 1; // 1 us window
    dbg = (debug_ && std::fabs(curWindowStartUs_ - TARGET_US) < EPS_US);

    double startUs = curWindowStartUs_;
    double endUs = curWindowEndUs_;
    
    std::vector<double> seeds_us = findCoincidenceSeeds(startUs, endUs, binWidth_us, pre_bins);
    if (dbg) {
        std::cerr << "[DBG] fitPulses() window#" << curWindowIndex_
                    << " coincidences found: " << seeds_us.size() << "\n";
        if (!seeds_us.empty()) {
            std::cerr << "  seeds (us rel start): ";
            for (double s : seeds_us) std::cerr << std::fixed << std::setprecision(2) << s << " ";
            std::cerr << "\n";
        }
    }

    if (seeds_us.empty()) {
        if (dbg) std::cerr << "  [DBG] no coincidences in window → skipping\n";
        return false;
    }

    std::vector<double> binEdges_us;
    binEdges_us.reserve(xCenters.size() + 1);
    for (double xc : xCenters)
        binEdges_us.push_back(xc - 0.5 * binWidth_us);
    binEdges_us.push_back(xCenters.back() + 0.5 * binWidth_us);

    std::vector<double> pmf_unit(pdfLookup.size());
    for (size_t i = 0; i < pdfLookup.size(); ++i)
        pmf_unit[i] = pdfLookup[i][0];

    FitResult fit = runGreedyLRT(
        std::vector<double>(hist.begin(), hist.end()),
        binEdges_us, pmf_unit, seeds_us, binWidth_us, /*max_offset_us=*/5.0,
        /*delta_nll_cut=*/10, /*d_min_us=*/2.0, /*amp_min_cut_pe=*/5.0,
        fixedExpected, peBackgroundRate_ * binWidth_us * 1e-6, dbg
    );

    if (fit.pulses.empty()) {
        if (dbg) std::cerr << "  [DBG] LRT found no pulses\n";
        return false;
    }

    fittedPEs.clear();
    fittedDTs.clear();
    for (auto& p : fit.pulses) {
        fittedDTs.push_back((p.t_us - 5.0) / binWidth_us); // shift +5us correction
        fittedPEs.push_back(p.amp_pe);
    }

    // --- Build expected histogram ---
    std::vector<double> expected(hist.size(), 0.0);

    // Fixed expected
    if (fixedExpected) {
        const auto& fix = *fixedExpected;
        for (size_t j = 0; j < expected.size(); ++j)
            expected[j] = (j < fix.size() ? std::max(0.0, fix[j]) : 0.0);
    }

    // Background add
    if (fitting_bg_ && peBackgroundRate_ > 0.0) {
        const double b = peBackgroundRate_ * binWidth_us * 1e-6;
        if (std::isfinite(b) && b > 0.0)
            for (double& e : expected) e += b;
    }

    // Add each pulse’s contribution
    for (size_t m = 0; m < fittedPEs.size(); ++m) {
        int dt_idx = std::clamp((int)std::llround(fittedDTs[m]), 0, (int)pdfLookup.size() - 1);
        const auto& row = pdfLookup[dt_idx];
        for (size_t j = 0; j < expected.size(); ++j)
            expected[j] += fittedPEs[m] * row[j];
    }

    // --- Evaluate log-likelihood ---
    auto poissonLogL = [&](const std::vector<int>& y, const std::vector<double>& mu){
        double s = 0.0;
        for (size_t j = 0; j < y.size(); ++j) {
            const int    k = y[j];
            const double lam = std::max(1e-12, mu[j]);
            s += k * std::log(lam) - lam - std::lgamma(k + 1.0);
        }
        return s;
    };
    const double logL = poissonLogL(hist, expected);

    int nObs = std::accumulate(hist.begin(), hist.end(), 0);
    double nExp = std::accumulate(expected.begin(), expected.end(), 0.0);

    if (dbg) {
        std::cerr << std::fixed
                  << "[DBG] win#" << curWindowIndex_
                  << " seg=" << segmentId_
                  << " pulses=" << fit.pulses.size()
                  << " logL=" << logL
                  << " nObs=" << nObs
                  << " nExp=" << nExp << "\n";
    }

    // --- Save WindowStat ---
    WindowStat ws;
    ws.windowIndex   = curWindowIndex_;
    ws.startTimeUs   = curWindowStartUs_;
    ws.endTimeUs     = curWindowEndUs_;
    ws.binWidthUs    = curWindowBinWidth_us_;
    ws.nPulsesChosen = (int)fittedPEs.size();
    ws.logL          = logL;
    ws.nObserved     = nObs;
    ws.nExpected     = nExp;
    windowStats_.push_back(ws);

    // --- Capture outliers ---
    if (captureOutliers_) {
        const double ratio = (nObs > 0 ? nExp / std::max(1.0, (double)nObs) : 0.0);
        if (nObs >= outlierMinObs_ && ratio <= outlierRatioLow_) {
            OutlierRecord rec;
            rec.windowIndex = curWindowIndex_;
            rec.startTimeUs = curWindowStartUs_;
            rec.binWidthUs  = curWindowBinWidth_us_;
            rec.nObserved   = nObs;
            rec.nExpected   = nExp;
            rec.ratioExpOverObs = ratio;
            rec.hist        = hist;
            rec.totalExpected = expected;
            outliers_.push_back(std::move(rec));
        }
    }

    return !fittedPEs.empty();
}