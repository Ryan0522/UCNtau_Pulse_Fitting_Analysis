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
            if (prob.bg_rate_hz > 0.0 && prob.bin_width_sec > 0.0 && prob.fit_bg) {
                const double b = prob.bg_rate_hz * prob.bin_width_sec;
                if (std::isfinite(b) && b > 0.0) {
                    expv[i] += b;
                }
            }
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

std::vector<std::vector<std::pair<double, double>>> 
cluster_seeds(const std::vector<std::pair<double,double>>& s, int closeBins) {
    std::vector<std::vector<std::pair<double,double>>> clusters;
    if (s.empty()) return clusters;
    std::vector<std::pair<double,double>> cur; cur.push_back(s[0]);
    for (size_t i=1;i<s.size();++i){
        if (std::fabs(s[i].first - s[i-1].first) <= closeBins) cur.push_back(s[i]);
        else { clusters.push_back(cur); cur.clear(); cur.push_back(s[i]); }
    }
    clusters.push_back(cur);
    return clusters;
}

static void init_k_weighted(
    const std::vector<std::pair<double, double>>& cl, 
    int K,
    std::vector<double>& pe0,
    std::vector<double>& dt0)
{
    K = std::max(1, K);
    pe0.clear(); dt0.clear();
    if (cl.empty()) { pe0.assign(K, 1.0); dt0.assign(K, 0.0); return; }

    double wsum = 0.0;
    for (const auto& s : cl) wsum += std::max(1.0, s.second);

    pe0.assign(K, std::max(1.0, wsum / K));
    dt0.resize(K);

    const double step = wsum / K;

    // two-pointer cumulative scan (using acc. weight to split dts)
    // following cl, find each index for (i+0.5)*step 
    double acc = 0.0; // accumulated weight
    size_t j = 0;// cl index
    for (int i=0; i<K; ++i) {
        const double target = (i + 0.5) * step;
        while (j + 1 < cl.size() && acc + std::max(1.0, cl[j].second) < target) {
            acc += std::max(1.0, cl[j].second);
            ++j;
        }
        dt0[i] = cl[j].first; // use dt0 near where the weight cut it
    }
}

std::pair<FitResult,int> select_k_for_cluster(
    const FitProblem& prob, const std::vector<std::pair<double,double>>& cl,
    double peMin, double peMax, double dtMin, double dtMax,
    const KSelectOptions& opt)
{
    FitResult best; int best_k = 0;
    const int m = (int)cl.size();
    for (int k=1;k<=m;++k){
        std::vector<double> pe0, dt0; init_k_weighted(cl, k, pe0, dt0);
        FitResult r = fit_n_pulses_bobyqa(prob, k, pe0, dt0, peMin, peMax, dtMin, dtMax, opt.maxEval);
        if (!r.ok) continue;
        double s = opt.useBIC ? (2.0*r.nll + (2*k)*std::log(std::max(1, prob.nTime))) : r.nll;
        double s_best = best.ok ? (opt.useBIC ? (2.0*best.nll + (2*best_k)*std::log(std::max(1, prob.nTime))) : best.nll) : std::numeric_limits<double>::infinity();
        if (!best.ok || s < s_best) { best=r; best_k=k; }
    }
    return {best, best_k};
}

FitResult fit_global_from_selections(
    const FitProblem& prob,
    const std::vector<double>& initPE,
    const std::vector<double>& initDT,
    double peMin, double peMax, double dtMin, double dtMax,
    int maxEval)
{
    const int n = (int)initPE.size();
    if (n==0) return {};
    return fit_n_pulses_bobyqa(prob, n, initPE, initDT, peMin, peMax, dtMin, dtMax, maxEval);
}