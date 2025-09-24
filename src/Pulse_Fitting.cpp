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

Pulse_Fitting::Pulse_Fitting(const EventList& events) : events_(events) { extractTimes(events); } // copy event.realtime to peTimes_us_ (us)

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
    cout << "Event size: " << peTimes_us_.size() << endl; // total PE hits loaded

    vector<double> signalTimes = applyTimeWindow(peTimes_us_, startAfterUs_, stopAfterUs_);
    vector<double> backgroundTimes;

    if (backgroundAfterUs_ > 0)
        backgroundTimes = applyTimeWindow(peTimes_us_, backgroundAfterUs_, backgroundAfterUs_ + 60e6);

    cout << "SignalTime PE Event size: " << signalTimes.size() << "  |  ";
    cout << "Background PE Event size: " << backgroundTimes.size() << endl;
    
    const double bg_duration_sec = 60.0;
    peBackgroundRate_ = std::max(0.0, backgroundTimes.size() / bg_duration_sec);
    fitting_bg_ = false;

    const int MAX_ITERS = 10;
    const double tolRel = 1e-2;

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

void Pulse_Fitting::extractTimes(const EventList& events) {
    // linearize event.realtime (s) -> vector of times (us)
    peTimes_us_.clear();
    peChans_.clear();
    peTimes_us_.reserve(events.size());
    peChans_.reserve(events.size());
    
    for (const auto& e : events) {
        peTimes_us_.push_back(e.realtime * 1e6);
        peChans_.push_back(e.channel);
    }
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
    // grow a window starting at 'startIdx' until an inter-hit gap > minGap_us_
    int N = static_cast<int>(times.size());
    if (startIdx >= N - 1) {
        return {0.0, startIdx + 1, 0.0, 0.0};
    }

    double start = times[startIdx];
    int j = startIdx + 1;

    if (use_coinc_) {
        bool armed = false;
        while (j < N) {
            double dt_seed = times[j] - times[startIdx];
            if (dt_seed > coinc_win_us_) break;
            if (events_[j].channel != events_[startIdx].channel) {
                armed = true;
                ++j;
                break;
            }
            ++j;
        }
        if (!armed) {
            const double t = times[startIdx];
            return {0.0, startIdx + 1, t, t};
        }
        int last = j -  1;
        while (j < N) {
            double dt_step = times[j] - times[last];
            if (dt_step >= minGap_us_) break;
            last = j;
            ++j;
        }
    } else {
        while (j < N && (times[j] - times[j - 1]) <= minGap_us_) {
            ++j;
        }
    }

    // previously: double end = times[j - 1]; // August 6th 2025
    // pad by minGap_us_ so the histogram includes explicit zeros after the last hit // August 10th 2025
    const double t_last = times[j - 1];
    const double nextHit = (j < N ? times[j] : numeric_limits<double>::infinity());
    const double end_pad = t_last + minGap_us_;
    const double cap_by_next = std::isfinite(nextHit) ? nextafter(nextHit, -numeric_limits<double>::infinity()) : end_pad;
    const double end = std::min(end_pad, cap_by_next);

    const double windowWidth = end - start;
    return make_tuple(windowWidth, j, start, end);
}

bool Pulse_Fitting::makeHistogram(const vector<double>& times, int i, double binWidth_us,
                                  double& windowWidth, int& j, double& startTime, double& endTime,
                                  vector<int>& hist, vector<double>& xCenters) 
{
    // compute [startTime, endTime] window and bin hits into 'hist' with given binWidth_us
    tie(windowWidth, j, startTime, endTime) = movingWindow(times, i);
    if (windowWidth < binWidth_us) return false;

    int nBins = static_cast<int>(ceil(windowWidth / binWidth_us));
    if (nBins < 1) return false;

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
            i = j;
            continue;
        }

        if (xCenters.size() < 2) {
            binWidth_us = fineBinWidth_us_;
            if (!makeHistogram(data_us, i, binWidth_us, windowWidth, j, startTime, endTime, hist, xCenters)) {
                i = j;
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

vector<int> Pulse_Fitting::findCoincidenceSeeds(double startUs, double endUs, double binWidthUs, int pre_bins) const
{
    std::vector<int> seeds;
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
            const double dt = tj - t0;
            if (dt > coinc_win_us_) break;
            if (!armed && peChans_[j] != ch0) armed = true;
            ++total;
        }

        const int seed_idx = pre_bins + (int)std::floor((t0 - startUs) / binWidthUs);
        const bool accepted = (armed && total >= coinc_seed_pe_min_);
        if (accepted) {
            seeds.push_back(seed_idx);
            prev_time_us = t0 + seeding_window_us_; // us
        }
    }

    std::sort(seeds.begin(), seeds.end());
    seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

    return seeds;
}

bool Pulse_Fitting::fitPulses(const vector<int>& hist, const vector<double>& xCenters,
                              const vector<vector<double>>& pdfLookup,
                              vector<double>& fittedPEs, vector<double>& fittedDTs,
                              const double binWidth_us,
                              const vector<double>* fixedExpected) 
{
    const int pre_bins = static_cast<int>(std::llround(shiftUs_ / binWidth_us)); // Since PDF is shifted by +5 us
    std::vector<int> peaks;

    const bool dbg = curWindowIndex_ == debug_window_index_ && segmentId_ == debug_segment_id_ && debug_;

    vector<double> peGuess;
    vector<double> dtGuess;

    const int W = std::max(1, (int)std::ceil(seeding_window_us_ / (binWidth_us + 1e-12)));

    if (use_coinc_) {
        const double startUs = curWindowStartUs_;
        const double endUs = curWindowStartUs_ + std::max(0, (int)xCenters.size() - pre_bins) * binWidth_us;
        peaks = findCoincidenceSeeds(startUs, endUs, binWidth_us, pre_bins);
        if (dbg) {
            std::vector<int> dbg_grad_peaks = findGradientPeaks(hist, gradThr_);
            std::cerr << "[COINC] window#" << curWindowIndex_
                    << " start=" << startUs << "us"
                    << " end=" << endUs << "us"
                    << " duration=" << (endUs - startUs) << "us"
                    << " binWidth=" << binWidth_us
                    << " pre_bins=" << pre_bins
                    << " nbins=" << xCenters.size()
                    << "\n";

            std::cerr << "[";
            for (auto b : xCenters) std::cerr << b << ' ';
            std::cerr << "]\n";

            std::cerr << "[";
            for (auto b : hist) std::cerr << b << ' ';
            std::cerr << "]\n";

            std::cerr << "  [GRAD] seed bins: ";
            for (int b : dbg_grad_peaks) std::cerr << b << ' ';
            std::cerr << "\n";
            std::cerr << "  seed bins: ";
            for (int b : peaks) std::cerr << b << ' ';
            std::cerr << "\n";
        }
    } else {
        peaks = findGradientPeaks(hist, gradThr_);
        peGuess.push_back(seedPE_);
        dtGuess.push_back(0.0);
    }

    for (int p : peaks) {
        if (p < pre_bins + guardBin_) {
            continue;
        };
        int idx = p - pre_bins;

        if (dbg) {
            std::cerr << "  [PEAK] raw peak at bin " << p << " -> dt guess at bin " << idx << "\n";
        }

        if (idx < 0 || idx >= (int)hist.size()) continue;

        int end = min(p + W + 1, static_cast<int>(hist.size())); // p here not idx because there are pre-bins, +1 because accumulate is exclusive at the end
        int sum = accumulate(hist.begin() + p, hist.begin() + end, 0); // p here not idx because there are pre-bins

        if (dbg) {
            std::cerr << "    [PEAK] sum from bin " << p << " to " << end - 1 << " = " << sum << "\n";
        }

        // if (sum >= static_cast<int>(max(1.0, std::ceil(pe_min_thresh_ / 2.0)))) {
        if (sum >= (int)pe_min_thresh_) {
            peGuess.push_back(sum);
            dtGuess.push_back(idx);
        }
    }

    if (peGuess.empty()) return false;
    
    if (dbg) {
        std::cerr << "  [INIT] raw pe/dt_bin/dt_us guesses (size=" << peGuess.size() << "): ";
        for (size_t i=0;i<peGuess.size();++i) std::cerr << "(" << peGuess[i] << "," << dtGuess[i] << "," << dtGuess[i] * binWidth_us << ") ";
        std::cerr << "\n";
    }

    FitProblem prob;
    prob.observed = &hist;
    prob.pdfLookup = &pdfLookup;
    prob.fixedExpected = fixedExpected;
    prob.nTime = (int)hist.size();

    // --- NEW: background rate configuration ---
    prob.bg_rate_hz = peBackgroundRate_;
    prob.bin_width_sec = binWidth_us * 1e-6;
    prob.fit_bg = fitting_bg_;
    // --- end NEW (August 31, 2025) ---

    double peMin = 1.0, peMax = 300.0;
    double dtMin = 0.0, dtMax = (double)xCenters.size() - 1.000001;

    std::vector<std::pair<double, double>> seeds;
    seeds.reserve(peGuess.size());
    for (size_t i=0;i<peGuess.size();++i) seeds.emplace_back(dtGuess[i], peGuess[i]);
    std::sort(seeds.begin(), seeds.end(), [](auto a, auto b){ return a.first < b.first; });

    auto clusters = cluster_seeds(seeds, binWidth_us, cluster_close_us_);

    kSelectOptions kopt;
    kopt.criterion = kSelectOptions::Criterion::AICc;
    kopt.use_local_T = true;
    kopt.lrt_min_delta = 6.0; // set 0.0 to disable
    kopt.maxEval = 200;

    size_t ci = 0;
    std::vector<double> chosenPEs, chosenDTs;
    const int pad_bins = std::max(1, (int)std::llround(2.0 / (binWidth_us + 1e-12))); // pad by 2 us on each side

    for (const auto& cl : clusters) {
        if (dbg) {
            std::cerr << "\n[CLUSTER] #" << ci << "  size=" << cl.size() << "\n";
            std::cerr << "  initial seeds (dt, pe): ";
            for (const auto& s : cl) std::cerr << "(" << s.first << ", " << s.second << ") ";
            std::cerr << "\n";
        }

        const double cmin = std::max(dtMin, cl.front().first - pad_bins);
        const double cmax = std::min(dtMax, cl.back().first + pad_bins);

        if (cl.size() == 1) {
            const double pe_final = std::max(1.0, cl[0].second);
            const double dt_final = std::clamp(cl[0].first, cmin, cmax);

            chosenPEs.push_back(pe_final);
            chosenDTs.push_back(dt_final);
            if (dbg) {
                std::cerr << "  [FINAL] singleton -> (dt, pe) = ("
                        << dt_final << ", " << pe_final << ")\n";
            }
            ++ci;
            continue;
        }

        // const int t0 = (int)std::ceil(cmin);
        // const int t1 = (int)std::floor(cmax);
        // FitProblem subprob = make_subproblem(prob, t0, t1);
        
        // std::vector<std::pair<double,double>> cl_local; cl_local.reserve(cl.size());
        // for (auto& s : cl) cl_local.push_back({ s.first - t0, s.second });  

        // const double dtMin_loc = 0.0; const double dtMax_loc = (double)(t1 - t0);

        // auto [best, kstart] = select_k_for_cluster(
        //     subprob, cl_local, peMin, peMax, dtMin_loc, dtMax_loc, kopt
        // );

        auto [best, kstart] = select_k_for_cluster(
            prob, cl, peMin, peMax, dtMin, dtMax, kopt
        );
        
        if (best.ok) {
            // for (auto& d : best.DTs) d += t0;
            chosenPEs.insert(chosenPEs.end(), best.PEs.begin(), best.PEs.end());
            chosenDTs.insert(chosenDTs.end(), best.DTs.begin(), best.DTs.end());
            if (dbg) {
                std::cerr << "  [FINAL] fit OK   k=" << kstart
                        << "  nll=" << best.nll
                        << "  -> pulses:\n";
                for (size_t i = 0; i < best.PEs.size(); ++i) {
                    std::cerr << "    (" << best.DTs[i] << ", " << best.PEs[i] << ")\n";
                }
            }
        } else {
            double pe_sum=0, dtw=0;
            for (auto& p:cl){ pe_sum += p.second; dtw += p.first * p.second; }
            const double pe_final = std::max(1.0, pe_sum);
            const double dt_final = std::clamp(cl.front().first, cmin, cmax);

            chosenPEs.push_back(pe_final);
            chosenDTs.push_back(dt_final);

            if (dbg) {
                const double dt_centroid = (pe_sum > 0.0) ? (dtw / pe_sum) : cl.front().first;
                std::cerr << "  [FINAL] fallback (fit failed)\n"
                        << "    pe_sum=" << pe_sum
                        << "    dt_first=" << cl.front().first
                        << "    dt_centroid=" << dt_centroid << " (not used)\n"
                        << "    -> (dt, pe) = (" << dt_final << ", " << pe_final << ")\n";
            }
        }
        ++ci;
    }
    
    if (dbg) {
        std::cerr << "  [SEL] chosen for global fit (size=" << chosenPEs.size() << "): ";
        for (size_t i=0;i<chosenPEs.size();++i) std::cerr << "(" << chosenPEs[i] << "," << chosenDTs[i] << ") ";
        std::cerr << "\n";
    }

    auto r = fit_global_from_selections(prob, chosenPEs, chosenDTs, peMin, peMax, dtMin, dtMax, 200);
    if (!r.ok || r.PEs.empty()) return false;

    if (dbg) {
        if (r.ok) {
            std::cerr << "  [GLOBAL] after fit_global_from_selections: ";
            for (size_t i=0;i<r.PEs.size();++i) std::cerr << "(" << r.PEs[i] << "," << r.DTs[i] << ") ";
            std::cerr << "\n";
        } else {
            std::cerr << "  [GLOBAL] fit_global_from_selections failed\n";
        }
    }

    auto prune = [&](const std::vector<double>& P, const std::vector<double>& D){
        std::vector<double> p2, d2;
        for (size_t i=0;i<P.size();++i) if (P[i] >= pe_min_thresh_) { p2.push_back(P[i]); d2.push_back(D[i]); }
        return std::make_pair(std::move(p2), std::move(d2));
    };

    {
        auto pr = prune(r.PEs, r.DTs);
        if (pr.first.empty()) return false;
        r.PEs.swap(pr.first); r.DTs.swap(pr.second);
    }

    for (int it=0; it<3; ++it) {
        auto r2 = fit_n_pulses_bobyqa(prob, (int)r.PEs.size(), r.PEs, r.DTs, peMin, peMax, dtMin, dtMax, 200);

        if (dbg) {
            if (r2.ok) {
                std::cerr << "  [FINAL LOOPING] after fit_n_pulses_bobyqa iteration " << it << ": ";
                for (size_t i=0;i<r2.PEs.size();++i) std::cerr << "(" << r2.PEs[i] << "," << r2.DTs[i] << ") ";
                std::cerr << "\n";
            } else {
                std::cerr << "  [FINAL LOOPING] fit_global_from_selections failed\n";
            }
        }

        if (!r2.ok) break;
        auto pr2 = prune(r2.PEs, r2.DTs);
        if ((int)pr2.first.size() == (int)r.PEs.size()) { r = r2; break; }
        if (pr2.first.empty()) return false;
        r.PEs.swap(pr2.first); r.DTs.swap(pr2.second);
    }

    fittedPEs = std::move(r.PEs);
    fittedDTs = std::move(r.DTs);

    if (dbg) {
        std::cerr << "  [FINAL] fitted pe/dt: ";
        for (size_t i=0;i<fittedPEs.size();++i) std::cerr << "(" << fittedPEs[i] << "," << fittedDTs[i] << ") ";
        std::cerr << "\n\n";
    }

    vector<double> expected(hist.size(), 0.0);
    if (fixedExpected) {
        const auto& fix = *fixedExpected;
        for (size_t j = 0; j < expected.size(); ++j){
            double v = (j < (int)fix.size() ? fix[j] : 0.0);
            expected[j] = std::isfinite(v) && v >= 0 ? v : 0.0;
            if (prob.bg_rate_hz > 0.0 && prob.bin_width_sec > 0.0 && prob.fit_bg) {
                const double b = prob.bg_rate_hz * prob.bin_width_sec;
                if (std::isfinite(b) && b > 0.0) {
                    expected[j] += b;
                }
            }
        }
    }

    for (size_t m = 0; m < fittedPEs.size(); ++m) {
        int dt_idx = (int)std::llround(fittedDTs[m]);
        dt_idx = std::clamp(dt_idx, 0, (int)pdfLookup.size()-1);
        const auto& row = pdfLookup[dt_idx];
        for (size_t j = 0; j < expected.size(); ++j)
            expected[j] += fittedPEs[m] * row[j]; // density -> mass distribution (Sept 1, 2025)
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
    ws.binWidthUs    = curWindowBinWidth_us_;
    ws.nPulsesChosen = (int)fittedPEs.size();
    ws.logL          = logL;
    ws.nObserved     = nObs;
    ws.nExpected     = nExp;

    windowStats_.push_back(ws);

    // == NEW: capture large deviations (lower-right outliers) ==
    if (captureOutliers_) {
        const double ratio = (nObs > 0 ? nExp / std::max(1.0, (double)nObs) : 0.0);
        if (nObs >= outlierMinObs_ && ratio <= outlierRatioLow_) {
            OutlierRecord rec;
            rec.windowIndex = curWindowIndex_;
            rec.startTimeUs = curWindowStartUs_;
            rec.binWidthUs = curWindowBinWidth_us_;
            rec.nObserved = nObs;
            rec.nExpected = nExp;
            rec.ratioExpOverObs = ratio;
            rec.hist = hist;
            rec.totalExpected = expected;
            outliers_.push_back(std::move(rec));
        }
    }
    // == end NEW ==

    return !fittedPEs.empty();
}