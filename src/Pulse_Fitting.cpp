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

using namespace std;

static inline int round_to_index(double x, int max_index_inclusive) {
    long long idx = std::llround(x);
    if (idx < 0) idx = 0;
    if (idx > (long long)max_index_inclusive) idx = (long long)max_index_inclusive;
    return static_cast<int>(idx);
}

Pulse_Fitting::Pulse_Fitting(const EventList& events) {extractTimes(events);} // copy event.realtime to peTimes_ (us)

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
    cout << "Event size: " << peTimes_.size() << endl; // total PE hits loaded

    vector<double> signalTimes = applyTimeWindow(peTimes_, startAfterUs_, stopAfterUs_);
    vector<double> backgroundTimes;

    if (backgroundAfterUs_ > 0)
        backgroundTimes = applyTimeWindow(peTimes_, backgroundAfterUs_, backgroundAfterUs_ + 60e6);

    cout << "SignalTime PE Event size: " << signalTimes.size() << "  |  ";
    cout << "Background PE Event size: " << backgroundTimes.size() << endl;
    
    fitRegion(signalTimes, signalPulses_); // parse windows, fit pulses
    fitRegion(backgroundTimes, backgroundPulses_); // ditto for background

    cout << "SignalTime Neutron Event count: " << signalPulses_.size() << "  |  ";
    cout << "Background Neutron Event count: " << backgroundPulses_.size() << "\n" << endl;

    peBackgroundRate_ = backgroundTimes.size() / 60.0;
    eventBackgroundRate_ = backgroundPulses_.size() / 60.0;
}

vector<double> Pulse_Fitting::buildCarryExpected(double winStartUs, double binWidth, const vector<double>& xCenters) const
{
    const int L = (int)xCenters.size();
    vector<double> carry(L, 0.0);
    if (L == 0) return carry;

    const vector<double> base = get_full_pdf(segmentId_, binWidth);
    if (base.empty()) return carry;
    const int fullBins = (int)base.size();
    const double pdf_len_us = fullBins * binWidth;

    for (const auto& pr : carryPulses_) {
        const double t0_us = pr.first; // absolute time (us)
        const double PE = pr.second;
        const double age = winStartUs - t0_us; // >= 0 means the peak of the pdf is to the left of current window

        if (age < 0.0) continue; // skip pulses in current window (if any) (this is for safety)
        if (age >= pdf_len_us) continue; // the entire PDF is to the left of current window

        const int dx = (int)std::floor(age / binWidth);
        const int copy_len = std::min(L, fullBins - dx);
        for (int j = 0; j < copy_len; ++j) {
            carry[j] += PE * base[j + dx];
        }
    }
    return carry;
}

void Pulse_Fitting::extractTimes(const EventList& events) {
    // linearize event.realtime (s) -> vector of times (us)
    vector<double> times;
    times.reserve(events.size());
    for (const auto& e : events) {
        times.push_back(e.realtime * 1e6);
    }
    peTimes_ = move(times);
}

vector<double> Pulse_Fitting::applyTimeWindow(const vector<double>& times, double start, double end) {
    vector<double> filtered_times;
    for (double t : times) {
        if (t >= start && t < end) {
            filtered_times.push_back(t);
        }
    }
    return filtered_times;
}

tuple<double, int, double, double> Pulse_Fitting::movingWindow(const vector<double>& times, int startIdx) {
    // grow a window starting at 'startIdx' until an inter-hit gap > minGap_
    int N = static_cast<int>(times.size());
    double start = times[startIdx];
    int j = startIdx + 1;

    while (j < N && (times[j] - times[j - 1]) <= minGap_) {
        ++j;
    }

    // previously: double end = times[j - 1]; // August 6th 2025
    // pad by minGap_ so the histogram includes explicit zeros after the last hit // August 10th 2025
    const double t_last = times[j - 1];
    const double nextHit = (j < N ? times[j] : numeric_limits<double>::infinity());
    const double end_pad = t_last + minGap_;
    const double cap_by_next = std::isfinite(nextHit) ? nextafter(nextHit, -numeric_limits<double>::infinity()) : end_pad;
    const double end = std::min(end_pad, cap_by_next);

    const double windowWidth = end - start;
    return make_tuple(windowWidth, j, start, end);
}

bool Pulse_Fitting::makeHistogram(const vector<double>& times, int i, double binWidth,
                                  double& windowWidth, int& j, double& startTime, double& endTime,
                                  vector<int>& hist, vector<double>& xCenters) 
{
    // compute [startTime, endTime] window and bin hits into 'hist' with given binWidth
    tie(windowWidth, j, startTime, endTime) = movingWindow(times, i);
    if (windowWidth < binWidth) return false;

    int nBins = static_cast<int>(ceil(windowWidth / binWidth));
    if (nBins < 1) return false;

    xCenters.resize(nBins);
    hist.assign(nBins, 0);
    for (int b = 0; b < nBins; ++b) {
        xCenters[b] = b * binWidth;
    }

    for (int k = i; k < j; ++k) { // fill counts per bin relative to startTime
        double t = times[k] - startTime;
        int bin = static_cast<int>(t / binWidth);
        if (bin >= 0 && bin < nBins) {
            hist[bin]++;
        }
    }

    // --- NEW: prepend 5 us worth of empty bins so fitter origin is at -5 us (due to PDF generate algorithm in Pulse_Tail.cpp) ---
    int pre_bins = (binWidth == binWidth_) ? preBins_ : finePreBins_;
    if (pre_bins > 0) {
        hist.insert(hist.begin(), pre_bins, 0);

        // rebuild xCenters to stay aligned (0, 1*binWidth, ...)
        vector<double>newX(hist.size());
        for (size_t b = 0; b < newX.size(); ++b) newX[b] = b * binWidth;
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
        double binWidth = binWidth_;

        if (!makeHistogram(data_us, i, binWidth, windowWidth, j, startTime, endTime, hist, xCenters)) {
            i = j;
            continue;
        }

        if (xCenters.size() < 2) {
            binWidth = fineBinWidth_;
            if (!makeHistogram(data_us, i, binWidth, windowWidth, j, startTime, endTime, hist, xCenters)) {
                i = j;
                continue;
            }
        }

        const vector<double> base = get_full_pdf(segmentId_, binWidth);
        const double pdf_len_us = base.empty() ? 0.0 : base.size() * binWidth;

        if (pdf_len_us > 0.0) {
            carryPulses_.erase(
                std::remove_if(carryPulses_.begin(), carryPulses_.end(),
                                [&](const auto& pr){ return (startTime - pr.first) >= pdf_len_us; }),
                carryPulses_.end()
            );
        }

        std::vector<double> carry = buildCarryExpected(startTime, binWidth, xCenters); // FIX
        vector<vector<double>> pdfLookup = generatePDFLookup(xCenters); // shifted PDFs cache
        vector<double> fittedPEs, fittedDTs;

        curWindowIndex_ = windowIndex;
        curWindowStartUs_ = startTime;
        curWindowBinWidth_ = binWidth;

        bool success = fitPulses(hist, xCenters, pdfLookup, fittedPEs, fittedDTs, binWidth, &carry);
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
                comps.reserve(PEs.size());
                vector<double> total(hist.size(), 0.0);
                for (size_t m = 0; m < PEs.size(); ++m) {
                    int dt_idx = (int)std::llround(DTs[m]);
                    dt_idx = std::max(0, std::min(dt_idx, (int)look.size() - 1));
                    vector<double> comp(hist.size(), 0.0);
                    const auto& row = look[dt_idx];
                    for (size_t j = 0; j < hist.size(); ++j) {
                        double val = PEs[m] * row[j];
                        comp[j] = val;
                        total[j] += val;
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
            double pulse_time_us = startTime + dt_bins * binWidth + shiftUs_; // +5us correction
            output.emplace_back(pulse_time_us, fittedPEs[k], windowIndex, windowWidth, fittedPEs.size() > 1, binWidth == fineBinWidth_); // store result
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
    const double binWidth = xCenters[1] - xCenters[0];

    auto key = std::make_pair(length, std::round(binWidth * 1e6) / 1e6);
    auto it = pdfCache_.find(key);
    if (it != pdfCache_.end()) return it->second;

    vector<vector<double>> lookup(length, vector<double>(length, 0.0));
    for (int dx = 0; dx < length; ++dx) {
        lookup[dx] = shifted_pdf(segmentId_, binWidth, dx, length);
    }

    pdfCache_[key] = lookup;
    return lookup;
}

vector<int> Pulse_Fitting::findGradientPeaks(const vector<int>& hist, double thresholdFactor) {
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

    vector<int> peaks;
    for (size_t i = 1; i + 1 < grad.size(); ++i) {
        if (grad[i] > threshold) {
            peaks.push_back(static_cast<int>(i));
        }
    }
    return peaks;
}

bool Pulse_Fitting::fitPulses(const vector<int>& hist, const vector<double>& xCenters,
                              const vector<vector<double>>& pdfLookup,
                              vector<double>& fittedPEs, vector<double>& fittedDTs,
                              const double binWidth,
                              const vector<double>* fixedExpected) 
{
    const int pre_bins = static_cast<int>(std::llround(shiftUs_ / binWidth));
    vector<int> peaks = findGradientPeaks(hist, gradThr_);

    vector<double> peGuess = { seedPE_ };
    vector<double> dtGuess = { 0.0 };
    const int window = 5, minPE = 5;

    for (int p : peaks) {
        if (p < pre_bins + guardBin_) continue;
        int idx = p - pre_bins;
        if (idx < 0 || idx > (int)hist.size()) continue;

        int end = min(idx + window, static_cast<int>(hist.size()));
        int sum = accumulate(hist.begin() + idx, hist.begin() + end, 0);

        if (sum >= minPE) {
            peGuess.push_back(sum);
            dtGuess.push_back(idx);
        }
    }

    vector<double> newPE, newDT;
    for (size_t i = 0; i < peGuess.size(); ++i) {
        if (peGuess[i] >= 5) {
            newPE.push_back(peGuess[i]);
            newDT.push_back(dtGuess[i]);
        }
    }

    if (newPE.empty()) {
        return false;
    }

    FitProblem prob;
    prob.observed = &hist;
    prob.pdfLookup = &pdfLookup;
    prob.fixedExpected = fixedExpected;
    prob.nTime = (int)hist.size();

    double peMin = 1.0, peMax = 300.0;
    double dtMin = 0.0, dtMax = (double)xCenters.size() - 1.000001;

    std::vector<std::pair<double, double>> seeds;
    for (size_t i=0;i<newPE.size();++i) seeds.emplace_back(newDT[i], newPE[i]);
    std::sort(seeds.begin(), seeds.end());

    ModelSelectOptions msopt;
    msopt.closeBins = 10;
    msopt.useBIC = false;
    
    FitResult r;
    bool tried_select = (seeds.size() >= 2) && (std::fabs(seeds[1].first - seeds[0].first) <= msopt.closeBins);
    if (tried_select) r = model_select_1_vs_2(prob, seeds, peMin, peMax, dtMin, dtMax, msopt);
    if (!r.ok) r = fit_n_pulses_bobyqa(prob, (int)newPE.size(), newPE, newDT, peMin, peMax, dtMin, dtMax, 200);
    if (!r.ok) return false;

    auto prune = [&](const std::vector<double>& P, const std::vector<double>& D){
        std::vector<double> p2, d2;
        for (size_t i=0;i<P.size();++i) if (P[i] >= 5.0) { p2.push_back(P[i]); d2.push_back(D[i]); }
        return std::make_pair(std::move(p2), std::move(d2));
    };

    {
        auto pr = prune(r.PEs, r.DTs);
        if (pr.first.empty()) return false;
        r.PEs.swap(pr.first); r.DTs.swap(pr.second);
    }

    for (int it=0; it<3; ++it) {
        auto r2 = fit_n_pulses_bobyqa(prob, (int)r.PEs.size(), r.PEs, r.DTs, peMin, peMax, dtMin, dtMax, 200);
        if (!r2.ok) break;
        auto pr2 = prune(r2.PEs, r2.DTs);
        if ((int)pr2.first.size() == (int)r.PEs.size()) { r = r2; break; } // 收斂
        if (pr2.first.empty()) return false;
        r.PEs.swap(pr2.first); r.DTs.swap(pr2.second);
    }

    fittedPEs = std::move(r.PEs);
    fittedDTs = std::move(r.DTs);

    vector<double> expected(hist.size(), 0.0);
    if (fixedExpected) {
        for (size_t j = 0; j < expected.size(); ++j) expected[j] = std::max(0.0, (*fixedExpected)[j]);
    }

    for (size_t m = 0; m < fittedPEs.size(); ++m) {
        int dt_idx = (int)std::llround(fittedDTs[m]);
        dt_idx = std::clamp(dt_idx, 0, (int)pdfLookup.size()-1);
        const auto& row = pdfLookup[dt_idx];
        for (size_t j = 0; j < expected.size(); ++j)
            expected[j] += fittedPEs[m] * row[j];
    }

    int nObs = 0; for (int c : hist) nObs += c;
    double nExp = 0.0; for (double v : expected) nExp += v;

    auto poissonLogL = [&](const std::vector<int>& y, const std::vector<double>& mu){
        double s = 0.0;
        for (size_t j=0;j<y.size();++j){
            const int    k = y[j];
            const double lam = std::max(1e-12, mu[j]);
            s += k * std::log(lam) - lam - std::lgamma(k + 1.0);
        }
        return s;
    };
    double logL = poissonLogL(hist, expected);

    WindowStat ws;
    ws.windowIndex   = curWindowIndex_;
    ws.startTimeUs   = curWindowStartUs_;
    ws.binWidthUs    = curWindowBinWidth_;
    ws.nPulsesChosen = (int)fittedPEs.size();
    ws.logL          = logL;
    ws.nObserved     = nObs;
    ws.nExpected     = nExp;

    windowStats_.push_back(ws);

    return !fittedPEs.empty();
}