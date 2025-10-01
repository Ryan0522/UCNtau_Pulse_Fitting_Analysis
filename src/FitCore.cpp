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

static inline void dbg_print_vec(const char* tag, const std::vector<double>& v, int maxn = 16) {
    std::cerr << "    " << tag << " [n=" << v.size() << "]: ";
    int n = (int)v.size();
    for (int i = 0; i < std::min(n, maxn); ++i) std::cerr << v[i] << (i+1<std::min(n,maxn) ? ' ' : '\0');
    if (n > maxn) std::cerr << " ...";
    std::cerr << "\n";
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
    const int T = prob.nTime;

    t0 = std::max(0, t0);
    t1 = std::min(T-1, t1);
    if (t1 < t0) return {}; // invalid
    const int L = (t1 - t0 + 1);

    auto* y_slice = new std::vector<int>(y.begin() + t0, y.begin() + t1 + 1);
    const auto* LU_global = prob.pdfLookup;

    const std::vector<double>* F_slice_ptr = nullptr;
    if (prob.fixedExpected) {
        auto* F_slice = new std::vector<double>(
            prob.fixedExpected->begin() + t0, prob.fixedExpected->begin() + t1 + 1);
        F_slice_ptr = F_slice;
    }

    FitProblem sub;
    sub.observed      = y_slice;
    sub.pdfLookup     = LU_global;
    sub.fixedExpected = F_slice_ptr;
    sub.nTime         = L;
    sub.bg_rate_hz    = prob.bg_rate_hz;
    sub.bin_width_sec = prob.bin_width_sec;
    sub.fit_bg        = prob.fit_bg;
    sub.windowIndex   = prob.windowIndex;
    sub.segmentId     = prob.segmentId;
    sub.t0_offset = t0;
    return sub;
}

double frac_interp(const std::vector<double>& frac, double dt_local)
{
    if (frac.empty()) return 1.0;
    const int n = (int)frac.size();
    if (n == 1) return frac[0];

    double x = dt_local;
    if (x <= 0.0) return frac[0];
    if (x >= n-1) return frac[n-1];

    int i = (int)std::floor(x);
    int j = i + 1;
    double w = x - i;
    return (1.0 - w) * frac[i] + w * frac[j];
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

    if (prob.bg_rate_hz > 0.0 && prob.bin_width_sec > 0.0 && prob.fit_bg) {
        const double b = prob.bg_rate_hz * prob.bin_width_sec;
        if (std::isfinite(b) && b > 0.0) {
            for (int i = 0; i < T; ++i) expv[i] += b;
        }
    }

    // params = [PE_0..PE_{n-1}, dt_0..dt_{n-1}]
    for (int i = 0; i < nPulses; ++i) {
        double PE = params[i];
        double dt_loc = params[nPulses + i];
        if (!(PE >= 0.0) || !std::isfinite(PE) || !std::isfinite(dt_loc)) return 1e300;

        double dt = dt_loc + (double)prob.t0_offset;

        if (dt < 0.0) dt = 0.0;
        if (dt > (double)(R - 1) - 1e-9) dt = (double)(R - 1) - 1e-9;

        int i0 = (int)std::floor(dt);
        double a = dt - i0; // [0,1)
        const auto& r0 = LU[i0];
        const auto& r1 = LU[i0 + 1];
        const int C0 = (int)r0.size();
        const int C1 = (int)r1.size();

        for (int b = 0; b < T; ++b) {
            const int g = prob.t0_offset + b;
            if (g < 0 || g >= C0 || g >= C1) continue; // outside support

            double row = (1.0 - a) * r0[g] + a * r1[g];
            double add = PE * row; // density -> mass distribution (Sept 1st, 2025)
            if (!std::isfinite(add) || add < 0) return 1e300;
            expv[b] += add;
        }
    }
    for (double& v : expv) if (v < 0.0) v = 1e-12;

    return -poissonLogL(y, expv);
}

static inline double clamp_inside(double v, double lo, double hi) {
    if (!(lo < hi)) return v;
    const double eps = 0.5; // Magic number
    if (v <= lo) return std::nextafter(lo + eps, hi);
    if (v >= hi) return std::nextafter(hi - eps, lo);
    return v;
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
    for (int i = 0; i < nPulses; ++i) {
        x.push_back( clamp_inside(initPE[i], peMin, peMax) );
    }
    for (int i = 0; i < nPulses; ++i) {
        x.push_back( clamp_inside(initDT[i], dtMin, dtMax) );
    }

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
        std::cerr << "[NLopt] optimize() failed in window "
              << prob.windowIndex << " (segment " << prob.segmentId << "): "
              << e.what() << "\n";
        // throw std::runtime_error("Fatal NLopt failure, exiting.");
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

    if (opt.debug) {
        std::cerr << "[kSel] win#" << opt.dbg_window_index
                  << " cluster#" << opt.dbg_cluster_index
                  << "  size=" << m
                  << "  dt in [" << dtMin << "," << dtMax << "]"
                  << "  pe in [" << peMin << "," << peMax << "]"
                  << "  T_eff=" << (opt.use_local_T ? (int)std::floor(dtMax) - (int)std::ceil(dtMin) + 1
                                                   : prob.nTime)
                  << "\n";
        std::cerr << "       seeds(dt,pe): ";
        for (auto& s : cl) std::cerr << "(" << s.first << "," << s.second << ") ";
        std::cerr << "\n";
    }

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

        if (opt.debug) {
            std::cerr << "  [kSel] try k=" << k << "  (T_eff=" << T_eff << ")\n";
            dbg_print_vec("dt0", dt0);
            dbg_print_vec("pe0", pe0);
        }

        FitResult r = fit_n_pulses_bobyqa(prob, k, pe0, dt0, peMin, peMax, dtMin, dtMax, opt.maxEval);
        
        if (!r.ok) {
            if (opt.debug) std::cerr << "    [kSel] fit failed @k=" << k << "\n";
            continue;
        }

        if (opt.debug) {
            // at-bounds flags
            int pe_atU=0, dt_atL=0, dt_atU=0;
            for (double v : r.PEs) if (std::abs(v - peMax) <= 1e-9*std::max(1.0,peMax)) ++pe_atU;
            for (double d : r.DTs) {
                if (std::abs(d - dtMin) <= 1e-9*std::max(1.0,std::abs(dtMin))) ++dt_atL;
                if (std::abs(d - dtMax) <= 1e-9*std::max(1.0,std::abs(dtMax))) ++dt_atU;
            }
            std::cerr << "    [kSel] nll=" << r.nll
                      << "  PEs_atUpper=" << pe_atU
                      << "  DTs_atLower=" << dt_atL
                      << "  DTs_atUpper=" << dt_atU << "\n";
            dbg_print_vec("DT", r.DTs);
            dbg_print_vec("PE", r.PEs);
        }
        
        bool pass_lrt = true;
        if (opt.lrt_min_delta > 0.0 && have_prev) {
            const double delta = 2.0*(prev_nll - r.nll); // improvement
            if (opt.debug) std::cerr << "    [kSel] LRT delta=" << delta
                                     << "  min=" << opt.lrt_min_delta
                                     << (delta >= opt.lrt_min_delta ? "  (keep)\n" : "  (reject)\n");
            pass_lrt = (delta >= opt.lrt_min_delta);
        }
        
        prev_nll = r.nll;
        have_prev = true;

        if (opt.debug) {
            const int p = 2 * T_eff;
            const double T = (opt.use_local_T ? std::max(1, prob.nTime) : std::max(1, prob.nTime));
            const double AIC  = r.nll + p;
            const double AICc = (T > p + 1) ? AIC + (2.0*p*(p+1)) / (T - p - 1) : std::numeric_limits<double>::infinity();
            const double BIC  = r.nll + 0.5 * p * std::log(T);
            const double SBIC = r.nll + 0.5 * opt.use_bic_lambda * p * std::log(T);
            std::cerr << "    [kSel] k_eff=" << T_eff
                    << "  NLL=" << r.nll
                    << "  AIC=" << AIC
                    << "  AICc=" << AICc
                    << "  BIC=" << BIC
                    << "  SoftBIC(λ=" << opt.use_bic_lambda << ")=" << SBIC << "\n";
        }  

        if (!pass_lrt) continue;

        const double s = score(r.nll, k, T_eff);
        const double s_best = best.ok 
            ? score(best.nll, best_k, T_eff) 
            : std::numeric_limits<double>::infinity();

        if (!best.ok || s < s_best) { best=r; best_k=k; }
    }

    if (opt.debug) {
        if (best.ok) {
            std::cerr << "  [kSel] BEST k=" << best_k << "  nll=" << best.nll << "\n";
            dbg_print_vec("DT*", best.DTs);
            dbg_print_vec("PE*", best.PEs);
        } else {
            std::cerr << "  [kSel] no valid k found\n";
        }
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