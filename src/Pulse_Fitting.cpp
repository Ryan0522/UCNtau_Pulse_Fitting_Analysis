#include "Pulse_Fitting.h"
#include "PDF_Global.h"
#include <numeric>
#include <algorithm>
#include <nlopt.hpp>
#include <iostream>
#include <limits>

using namespace std;

std::vector<double> makeLogFactorialTable(int max_k) {
    // precompute log(k!) for small k (speed Poisson logL at low counts)
    std::vector<double> table(max_k);
    for (int k = 0; k < max_k; ++k) {
        table[k] = std::lgamma(k + 1.0);
    }
    return table;
}

PDFParams pdfParams_ = {
    1.09453333e+03,
    5.32077446e+03,
    9.93074362e+03,
    3.65357381e-01,
    2.77520732e+00,
    2.30165740e+01,
    -4.84253484e-03
};

const int MAX_K = 1000;
std::vector<double> log_fact_table = makeLogFactorialTable(MAX_K);

Pulse_Fitting::Pulse_Fitting(const EventList& events, double binWidth, double minGap)
    : binWidth_(binWidth), minGap_(minGap), fineBinWidth_(0.25),
      startAfterUs_(0), stopAfterUs_(1e12), backgroundAfterUs_(-1),
      peBackgroundRate_(0), eventBackgroundRate_(0) {extractTimes(events);} // copy event.realtime to peTimes_ (us)

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

    return true;
}

void Pulse_Fitting::fitRegion(const vector<double>& data_us,
                              vector<tuple<double, double, int, double, bool>>& output) 
{
    // slide over data, window by window, fit pulses per window
    int i = 0;
    int N = static_cast<int>(data_us.size());
    int windowCount = 0;
    
    while (i < N) {
        vector<int> hist;
        vector<double> xCenters;
        double windowWidth, startTime, endTime;
        int j;

        if (!makeHistogram(data_us, i, binWidth_, windowWidth, j, startTime, endTime, hist, xCenters)) {
            i = j;
            continue;
        }

        if (xCenters.size() < 2) {
            if (!makeHistogram(data_us, i, fineBinWidth_, windowWidth, j, startTime, endTime, hist, xCenters)) {
                i = j;
                continue;
            }
        }

        vector<vector<double>> pdfLookup = generatePDFLookup(xCenters); // shifted PDFs cache

        vector<double> fittedPEs, fittedDTs;

        bool success = fitPulses(hist, xCenters, pdfLookup, fittedPEs, fittedDTs);
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

            cout << "PE: " << endl;
            for (const auto& val : fittedPEs) cout << val << ", ";
            cout << "\n" << endl;
            cout << "DT: " << endl;
            for (const auto& val : fittedDTs) cout << val << ", ";
            cout << "\n" << endl;

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

        for (size_t k = 0; k < fittedPEs.size(); ++k) {
            double pulse_time_us = startTime + fittedDTs[k] * (xCenters[1] - xCenters[0]);
            output.emplace_back(pulse_time_us, fittedPEs[k], windowCount, windowWidth, fittedPEs.size() > 1); // store result
            // cout << (double)j/(double)N << ", " << pulse_time_us / 1e6 << ", " << fittedPEs[k] << ", " << endl;
        }

        windowCount++;
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

double Pulse_Fitting::poissonLogLikelihood(const vector<int>& observed, const vector<double>& expected) {
    // logL = sum_k [ k*log(lam) - lam - log(k!) ]; Stirling for large k
    double logL = 0.0;
    for (size_t i = 0; i < observed.size(); ++i) {
        double lam = expected[i] + 1e-10;
        int k = observed[i];
        if (k < 20) {
            logL += k * log(lam) - lam - log_fact_table[k];
        } else {
            logL += k * log(lam) - lam - (k * log(k) - k + 0.5 * log(2 * M_PI * k));
        }
    }
    return logL;
}

double Pulse_Fitting::negLogLikelihood(const vector<double>& params, const vector<int>& observed,
                                        const vector<vector<double>>& pdfLookup, int nPulses) 
{
    // params = [PE_0..PE_{n-1}, dt_0..dt_{n-1}] ; expected = sum_i PE_i * shiftedPDF(dt_i)
    vector<double> expected(observed.size(), 0.0);
    for (int i = 0; i < nPulses; ++i) {
        double PE = params[i];
        int dt = static_cast<int>(params[nPulses + i]);
        for (size_t j = 0; j < observed.size(); ++j) {
            expected[j] += PE * pdfLookup[dt][j];
        }
    }
    return -poissonLogLikelihood(observed, expected);
}

vector<int> Pulse_Fitting::findGradientPeaks(const vector<int>& hist, double thresholdFactor, int ignoreIdx) {
    // simple gradient-based seed find; thresholdFactor in units of grad "std"
    if ((int)hist.size() <= ignoreIdx + 2) {
        return {};
    }

    vector<double> grad(hist.size() - ignoreIdx);
    for (size_t i = ignoreIdx; i + 1 < hist.size(); ++i) {
        grad[i - ignoreIdx] = static_cast<double>(hist[i + 1] - hist[i - 1]) / 2.0;
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
                              vector<double>& fittedPEs, vector<double>& fittedDTs) 
{
    // seed candidates from gradient; then NLOpt (bounded) to fit PE, dt
    const int minPE = 5;
    const int window = 5;
    const int ignoreIdx = 3;

    // Use separate gradient peak detection function
    vector<int> peaks = findGradientPeaks(hist, 2.0, ignoreIdx);

    vector<double> peGuess = {20.0};
    vector<double> dtGuess = {0.0};
    for (int p : peaks) {
        int idx = p + ignoreIdx;
        int start = idx;
        int end = min(idx + window, static_cast<int>(hist.size()));
        int sum = accumulate(hist.begin() + start, hist.begin() + end, 0);

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

    int nPulses = static_cast<int>(newPE.size());
    vector<double> params;
    params.insert(params.end(), newPE.begin(), newPE.end());
    params.insert(params.end(), newDT.begin(), newDT.end());

    vector<double> lb, ub;
    for (int i = 0; i < nPulses; ++i) {
        lb.push_back(1.0);
        ub.push_back(300.0);
    }
    for (int i = 0; i < nPulses; ++i) {
        lb.push_back(0.0);
        ub.push_back(static_cast<double>(xCenters.size() - 1));
    }

    nlopt::opt opt(nlopt::LN_BOBYQA, params.size()); // derivative-free local
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    auto objective = [&](const vector<double> &x, vector<double> &grad) {
        return negLogLikelihood(x, hist, pdfLookup, nPulses);
    };

    opt.set_min_objective([](const vector<double> &x, vector<double> &grad, void *data) -> double {
        return (*static_cast<decltype(objective)*>(data))(x, grad);
    }, &objective);

    opt.set_xtol_rel(1e-4); // relative tolerance
    opt.set_maxeval(200); // iteration cap

    double minf;
    try {
        nlopt::result result = opt.optimize(params, minf);
    } catch (exception& e) {
        cerr << "NLopt failed: " << e.what() << endl;
        return false;
    }

    fittedPEs.assign(params.begin(), params.begin() + nPulses);
    fittedDTs.assign(params.begin() + nPulses, params.end());

    vector<double> finalPEs, finalDTs;
    for (size_t i = 0; i < fittedPEs.size(); ++i) {
        if (fittedPEs[i] >= 5) {
            finalPEs.push_back(fittedPEs[i]);
            finalDTs.push_back(fittedDTs[i]);
        }
    }

    if (finalPEs.empty()) {
        return false;
    }

    if (finalPEs.size() < fittedPEs.size()) {
        // drop sub-threshold pulses and re-fit the reduced model
        int refinedN = static_cast<int>(finalPEs.size());
        vector<double> refinedParams = finalPEs;
        refinedParams.insert(refinedParams.end(), finalDTs.begin(), finalDTs.end());

        lb.clear(); ub.clear();
        for (int i = 0; i < refinedN; ++i) {
            lb.push_back(1.0); ub.push_back(300.0);
        }
        for (int i = 0; i < refinedN; ++i) {
            lb.push_back(0.0); ub.push_back(static_cast<double>(xCenters.size() - 1));
        }

        nlopt::opt opt2(nlopt::LN_BOBYQA, refinedParams.size());
        opt2.set_lower_bounds(lb);
        opt2.set_upper_bounds(ub);
        auto refinedObj = [&](const vector<double> &x, vector<double> &grad) {
            return negLogLikelihood(x, hist, pdfLookup, refinedN);
        };

        opt2.set_min_objective([](const vector<double> &x, vector<double> &grad, void *data) -> double {
            return (*static_cast<decltype(refinedObj)*>(data))(x, grad);
        }, &refinedObj);

        opt2.set_xtol_rel(1e-4);
        opt2.set_maxeval(200);

        try {
            double refinedMinf;
            opt2.optimize(refinedParams, refinedMinf);
            fittedPEs.assign(refinedParams.begin(), refinedParams.begin() + refinedN);
            fittedDTs.assign(refinedParams.begin() + refinedN, refinedParams.end());
        } catch (exception& e) {
            cerr << "Refined NLopt failed: " << e.what() << endl;
            return false;
        }
    } else {
        fittedPEs = finalPEs;
        fittedDTs = finalDTs;
    }

    return true;
}