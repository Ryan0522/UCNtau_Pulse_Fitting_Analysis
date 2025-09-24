#include "FitCore.h"
#include <nlopt.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>

std::vector<double> g_log_fact = [](){
    const int MAX_K = 5000;
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
                                const FitProblem& prob,
                                const std::vector<std::vector<double>>& LU, // prob.pdfLookup deref
                                const char* tag) 
{
    std::cerr << "\n=== NLopt debug dump: " << tag << " ===\n";
    std::cerr << "nPulses=" << nPulses
              << "  nParams=" << params.size()
              << " (PE[" << nPulses << "], DT[" << nPulses << "])\n";

    for (size_t i = 0; i < params.size(); ++i) {
        const bool atL = std::abs(params[i] - lb[i]) <= 1e-12 * std::max(1.0, std::abs(lb[i]));
        const bool atU = std::abs(params[i] - ub[i]) <= 1e-12 * std::max(1.0, std::abs(ub[i]));
        const char* kind = (i < (size_t)nPulses) ? "PE" : "DT";
        std::cerr << "  [" << std::setw(2) << i << "] " << kind
                  << " : " << params[i]
                  << "  bounds: [" << lb[i] << ", " << ub[i] << "]"
                  << (atL ? "  (AT-LOWER)" : "")
                  << (atU ? "  (AT-UPPER)" : "")
                  << "\n";
    }

    const auto& y = *prob.observed;
    const int T   = prob.nTime;
    std::vector<double> expv(T, 0.0);

    if (prob.fixedExpected) {
        const auto& fix = *prob.fixedExpected;
        for (int t = 0; t < T; ++t) {
            double v = (t < (int)fix.size() ? fix[t] : 0.0);
            expv[t] = (std::isfinite(v) && v >= 0) ? v : 0.0;
            if (prob.fit_bg && prob.bg_rate_hz > 0.0 && prob.bin_width_sec > 0.0) {
                expv[t] += prob.bg_rate_hz * prob.bin_width_sec;
            }
        }
    }

    const int R = (int)LU.size();
    for (int k = 0; k < nPulses; ++k) {
        double PE = params[k];
        double dt = params[nPulses + k];
        if (!(PE >= 0.0) || !std::isfinite(PE) || !std::isfinite(dt)) continue;

        if (dt < 0.0) dt = 0.0;
        if (dt > (double)(R - 1) - 1e-9) dt = (double)(R - 1) - 1e-9;
        int i0 = (int)std::floor(dt);
        double a = dt - i0;
        const auto& r0 = LU[i0];
        const auto& r1 = LU[i0 + 1];

        for (int t = 0; t < T; ++t) {
            double row = (1.0 - a) * r0[t] + a * r1[t];
            expv[t] += PE * row;
        }
    }

    auto sum_vec = [](const auto& v){ double s=0; for (auto x: v) s+=x; return s; };
    auto minmax  = [](const auto& v){ double mn=1e300,mx=-1e300; for(auto x:v){mn=std::min(mn,x); mx=std::max(mx,x);} return std::pair{mn,mx}; };

    const double ysum   = sum_vec(y);
    const double msum   = sum_vec(expv);
    const auto [emin, emax] = minmax(expv);

    std::cerr << "bin_width_sec=" << prob.bin_width_sec
              << "  bg_rate_hz=" << prob.bg_rate_hz
              << "  fit_bg=" << (prob.fit_bg ? "true" : "false") << "\n";
    std::cerr << "observed: sum=" << ysum << "  nTime=" << T << "\n";
    std::cerr << "model   : sum=" << msum << "  min/max=(" << emin << ", " << emax << ")\n";

    std::cerr << "=== end dump ===\n\n";
}

FitProblem make_subproblem(const FitProblem& prob, int t0, int t1) {
    const auto& y = *prob.observed;
    const auto& LU = *prob.pdfLookup;
    const int T = prob.nTime;

    t0 = std::max(0, t0);
    t1 = std::min(T-1, t1);
    const int L = (t1 - t0 + 1);
    if (L <= 0) return {}; // invalid

    auto* y_slice = new std::vector<int>(y.begin() + t0, y.begin() + t1 + 1);
    auto* LU_slice = new std::vector<std::vector<double>>(LU.size());
    for (size_t r = 0; r < LU.size(); ++r) {
        const auto& row = LU[r];
        (*LU_slice)[r].assign(row.begin() + t0, row.begin() + t1 + 1);
    }

    const std::vector<double>* fix_ptr = nullptr;
    auto* fix = new std::vector<double>();
    const auto& F = *prob.fixedExpected;
    const int Fn = (int)F.size();
    const int a = std::min(std::max(0, t0), Fn);
    const int b = std::min(std::max(0, t1+1), Fn);
    fix->assign(F.begin() + a, F.begin() + b);
    fix_ptr = fix;

    FitProblem sub;
    sub.observed      = y_slice;
    sub.pdfLookup     = LU_slice;
    sub.fixedExpected = fix_ptr;
    sub.nTime         = L;
    sub.bg_rate_hz    = prob.bg_rate_hz;
    sub.bin_width_sec = prob.bin_width_sec;
    sub.fit_bg        = prob.fit_bg;
    return sub;
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
            double add = PE * row; // density -> mass distribution (Sept 1st, 2025)
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
        dump_params_bounds(x, lb, ub, nPulses, prob, *prob.pdfLookup, "on-error");
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
cluster_seeds(const std::vector<std::pair<double,double>>& s, double binWidth_us, double cluster_close_us) {
    std::vector<std::vector<std::pair<double,double>>> clusters;
    if (s.empty()) return clusters;
    std::vector<std::pair<double,double>> cur; cur.push_back(s[0]);
    for (size_t i=1;i<s.size();++i){
        if (std::fabs((s[i].first - s[i-1].first) * binWidth_us) <= cluster_close_us) cur.push_back(s[i]);
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
    if (cl.empty()) { pe0.assign(K,1.0); dt0.assign(K,0.0); return; }

    struct S { double dt, w; };
    std::vector<S> v; v.reserve(cl.size());
    for (auto& p : cl) v.push_back({ p.first, std::max(1.0, p.second) });

    if (v.size() > (size_t)K) {
        std::nth_element(v.begin(), v.begin()+K, v.end(),
                         [](const S& a, const S& b){ return a.w > b.w; });
        v.resize(K);
    }

    std::sort(v.begin(), v.end(), [](const S& a, const S& b){ return a.dt < b.dt; });
    const double eps = 0.25;
    for (int i = 1; i < (int)v.size(); ++i)
        if (v[i].dt <= v[i-1].dt) v[i].dt = v[i-1].dt + eps;

    while ((int)v.size() < K) v.push_back({ v.empty()?0.0 : v[0].dt + eps*v.size(),
                                            v.empty()?1.0 : v[0].w });

    pe0.resize(K); dt0.resize(K);
    for (int i = 0; i < K; ++i) { dt0[i] = v[i].dt; pe0[i] = std::max(1.0, v[i].w); }
}

std::pair<FitResult,int> select_k_for_cluster(
    const FitProblem& prob, const std::vector<std::pair<double,double>>& cl,
    double peMin, double peMax, double dtMin, double dtMax,
    const kSelectOptions& opt)
{
    FitResult best; int best_k = 0;
    const int m = (int)cl.size();

    auto score = [&](double nll, int k, int T_eff) {
        const int p = 2*k; // # parameters, PE and DT per pulse
        switch (opt.criterion) {
            case kSelectOptions::Criterion::NLL:
                return nll;
            case kSelectOptions::Criterion::AIC:
                return 2.0*nll + 2.0*p;
            case kSelectOptions::Criterion::AICc: {
                // fallback to AIC when T_eff too small
                if (T_eff <= p + 1) return 2.0*nll + 2.0*p;
                return 2.0*nll + 2.0*p + (2.0*p*(p+1)) / (T_eff - p - 1);
            }
            case kSelectOptions::Criterion::BIC:
                return 2.0*nll + p*std::log(std::max(1, T_eff));
            case kSelectOptions::Criterion::SoftBIC:
                return 2.0*nll + opt.use_bic_lambda * p*std::log(std::max(1, T_eff));
        }
        return std::numeric_limits<double>::infinity();
    };

    double prev_nll = std::numeric_limits<double>::infinity();
    bool have_prev = false;

    for (int k=1;k<=m;++k){
        std::vector<double> pe0, dt0; 
        init_k_weighted(cl, k, pe0, dt0);
        const int T_eff = opt.use_local_T 
            ? std::max(1, (int)std::floor(dtMax) - (int)std::ceil(dtMin) + 1)
            : std::max(1, prob.nTime);

        FitResult r = fit_n_pulses_bobyqa(prob, k, pe0, dt0, peMin, peMax, dtMin, dtMax, opt.maxEval);
        if (!r.ok) continue;
        
        bool pass_lrt = true;
        if (opt.lrt_min_delta > 0.0 && have_prev) {
            const double delta = 2.0*(prev_nll - r.nll); // improvement
            pass_lrt = (delta >= opt.lrt_min_delta);
        }
        
        prev_nll = r.nll;
        have_prev = true;

        if (!pass_lrt) continue;

        const double s = score(r.nll, k, T_eff);
        const double s_best = best.ok 
            ? score(best.nll, best_k, T_eff) 
            : std::numeric_limits<double>::infinity();

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