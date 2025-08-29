#include "FitCore.h"
#include <nlopt.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>

std::vector<double> g_log_fact = [](){
    const int MAX_K = 1000;
    std::vector<double> table(MAX_K);
    for (int k = 0; k < MAX_K; ++k) {
        table[k] = std::lgamma(k + 1.0);
    }
    return table;
}();

std::vector<double> makeLogFactorialTable(int max_k) {
    std::vector<double> t(max_k);
    for (int k = 0; k < max_k; ++k) t[k] = std::lgamma(k + 1.0);
    return t;
}

static inline double poissonLogL(const std::vector<int>& k, const std::vector<double>& lam) {
    // logL = sum_k [ k*log(lam) - lam - log(k!) ]; Stirling for large k
    double logL = 0.0;
    for (size_t i = 0; i < k.size(); ++i) {
        double L = lam[i] + 1e-12;
        int    c = k[i];
        if (c < (int)g_log_fact.size()) {
            logL += c * std::log(L) - L - g_log_fact[c];
        } else {
            // Stirling
            logL += c * std::log(L) - L - (c*std::log(c) - c + 0.5*std::log(2*M_PI*c));
        }
    }
    return logL;
}

static inline void dump_params_bounds(const std::vector<double>& params,
                                const std::vector<double>& lb,
                                const std::vector<double>& ub,
                                int nPulses,
                            const char* tag) 
{
    std::cerr << "\n=== NLopt debug dump: " << tag << " ===\n";
    std::cerr << "nPulses=" << nPulses
              << "  nParams=" << params.size() << " (PE[" << nPulses
              << "], DT[" << nPulses << "])\n";
    for (size_t i = 0; i < params.size(); ++i) {
        std::string kind = (i < (size_t)nPulses) ? "PE" : "DT";
        std::cerr << "  [" << std::setw(2) << i << "] " << kind
                  << " : " << params[i]
                  << "  bounds: [" << lb[i] << ", " << ub[i] << "]\n";
    }
    std::cerr << "=== end dump ===\n\n";
}

static double negLogL_core(const std::vector<double>& params,
                           const FitProblem& prob,
                           int nPulses)
{
    const auto& y  = *prob.observed;
    const auto& LU = *prob.pdfLookup;
    const int T = prob.nTime;
    const int R = (int)LU.size();
    if (T <= 0 || R < 2) return 1e300;

    std::vector<double> expv(T, 0.0);

    if (prob.fixedExpected) {
        const auto& fix = *prob.fixedExpected;
        for (int i = 0; i < T; ++i) {
            double v = (i < (int)fix.size() ? fix[i] : 0.0);
            expv[i] = std::isfinite(v) && v >= 0 ? v : 0.0;
        }
    }

    // params = [PE_0..PE_{n-1}, dt_0..dt_{n-1}]
    for (int i = 0; i < nPulses; ++i) {
        double PE = params[i];
        double dt = params[nPulses + i];
        if (!(PE >= 0.0) || !std::isfinite(PE) || !std::isfinite(dt)) return 1e300;

        if (dt < 0.0) dt = 0.0;
        if (dt > (double)(R - 1) - 1e-9) dt = (double)(R - 1) - 1e-9;

        int i0 = (int)std::floor(dt);
        double a = dt - i0; // [0,1)
        const auto& r0 = LU[i0];
        const auto& r1 = LU[i0 + 1];

        for (int t = 0; t < T; ++t) {
            double row = (1.0 - a) * r0[t] + a * r1[t];
            double add = PE * row;
            if (!std::isfinite(add) || add < 0) return 1e300;
            expv[t] += add;
        }
    }
    for (double& v : expv) if (v < 0.0) v = 1e-12;

    return -poissonLogL(y, expv);
}

FitResult fit_n_pulses_bobyqa(
    const FitProblem& prob,
    int nPulses,
    const std::vector<double>& initPE,
    const std::vector<double>& initDT,
    double peMin, double peMax,
    double dtMin, double dtMax,
    int maxEval)
{
    FitResult out;
    if (nPulses <= 0) return out;

    std::vector<double> x; x.reserve(2*nPulses);
    x.insert(x.end(), initPE.begin(), initPE.begin()+nPulses);
    x.insert(x.end(), initDT.begin(), initDT.begin()+nPulses);

    std::vector<double> lb, ub; lb.reserve(2*nPulses); ub.reserve(2*nPulses);
    for (int i=0;i<nPulses;++i){ lb.push_back(peMin); ub.push_back(peMax); }
    for (int i=0;i<nPulses;++i){ lb.push_back(dtMin); ub.push_back(dtMax); }

    nlopt::opt opt(nlopt::LN_BOBYQA, (int)x.size());
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    auto obj = [&](const std::vector<double>& xx, std::vector<double>& g){
        (void)g;
        return negLogL_core(xx, prob, nPulses);
    };
    
    opt.set_min_objective([](const std::vector<double> &xx, std::vector<double> &gg, void *data) -> double {
        return (*static_cast<decltype(obj)*>(data))(xx, gg);
    }, &obj);

    opt.set_xtol_rel(1e-4);
    opt.set_maxeval(maxEval);

    double fmin=0.0;
    try {
        nlopt::result result = opt.optimize(x, fmin);
    } catch (const std::exception& e) {
        dump_params_bounds(x, lb, ub, nPulses, "on-error");
        std::cerr << "[NLopt] optimize() failed: " << e.what() << "\n";
        throw std::runtime_error("Fatal NLopt failure, exiting.");
        return out;
    }

    out.ok = true;
    out.nll = fmin;
    out.PEs.assign(x.begin(), x.begin()+nPulses);
    out.DTs.assign(x.begin()+nPulses, x.end());
    return out;
}

static inline double score(double nll, int k, int T, bool useBIC) {
    if (!useBIC) return nll;
    return 2.0*nll + k*std::log(std::max(1, T));
}

FitResult model_select_1_vs_2(
    const FitProblem& prob,
    const std::vector<std::pair<double,double>>& seeds, // sorted
    double peMin, double peMax, double dtMin, double dtMax,
    const ModelSelectOptions& opt)
{
    FitResult best;

    if (seeds.size() < 2) return best;
    
    double ddt = std::fabs(seeds[1].first - seeds[0].first);
    if (ddt > opt.closeBins) return best; // too far

    std::vector<double> pe2 = {seeds[0].second, seeds[1].second};
    std::vector<double> dt2 = {seeds[0].first,  seeds[1].first};

    double pe1 = pe2[0] + pe2[1];
    double dt1 = seeds[0].first;

    auto r1 = fit_n_pulses_bobyqa(prob, 1, {pe1}, {dt1}, peMin, peMax, dtMin, dtMax, 200);
    auto r2 = fit_n_pulses_bobyqa(prob, 2, pe2, dt2, peMin, peMax, dtMin, dtMax, 200);

    int T = prob.nTime;
    if (r1.ok) {
        double s = score(r1.nll, 2, T, opt.useBIC);
        best = r1; best.nll = s;
    }
    if (r2.ok) {
        double s = score(r2.nll, 4, T, opt.useBIC);
        if (!best.ok || s < best.nll) { best = r2; best.nll = s; }
    }

    best.ok = (r1.ok || r2.ok);
    return best;
}