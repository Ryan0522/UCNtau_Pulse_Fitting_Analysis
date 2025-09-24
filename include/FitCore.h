#pragma once
#include <vector>
#include <utility>
#include <memory>

struct kSelectOptions {
    enum class Criterion { NLL, AIC, AICc, BIC, SoftBIC };
    Criterion criterion = Criterion::AICc;
    double use_bic_lambda = 0.7;
    bool use_local_T = true;
    double lrt_min_delta = 0.0;
    int maxEval = 200;

    bool debug = false;
    int dbg_window_index = -1;
    int dbg_cluster_index = -1;
};

struct FitProblem {
    const std::vector<int>* observed;
    const std::vector<std::vector<double>>* pdfLookup;
    const std::vector<double>* fixedExpected;
    int nTime;
    double bg_rate_hz = 0.0;
    double bin_width_sec = 0.0;
    bool fit_bg = false;

    int windowIndex = -1;
    int segmentId = -1;
};

struct SubproblemWithFrac {
    FitProblem prob;
    std::vector<double> frac_in_slice;
    std::shared_ptr<std::vector<std::vector<double>>> owned_pdf;
};

struct FitResult {
    bool ok = false;
    std::vector<double> PEs;
    std::vector<double> DTs;
    double nll = 1e300;
};

std::vector<double> makeLogFactorialTable(int max_k);
extern std::vector<double> g_log_fact;

FitResult fit_n_pulses_bobyqa(
    const FitProblem& prob,
    int nPulses,
    const std::vector<double>& initPE,
    const std::vector<double>& initDT,
    double peMin, double peMax,
    double dtMin, double dtMax,
    int maxEval = 200
);

FitProblem make_subproblem(const FitProblem& prob, int t0, int t1);

SubproblemWithFrac make_subproblem_masscomp(const FitProblem& full, int t0, int t1);
// Linear interpolation of f for fractional dt (clamped to [0, n-1])
double frac_interp(const std::vector<double>& frac, double dt_local);

std::vector<std::vector<std::pair<double, double>>> cluster_seeds(const std::vector<std::pair<double, double>>& seeds_sorted_by_dt, double binWidth_us, double cluster_close_us);

std::pair<FitResult,int> select_k_for_cluster(
    const FitProblem& prob,
    const std::vector<std::pair<double,double>>& cluster,
    double peMin, double peMax, double dtMin, double dtMax,
    const kSelectOptions& opt
);

FitResult fit_global_from_selections(
    const FitProblem& prob,
    const std::vector<double>& initPE,
    const std::vector<double>& initDT,
    double peMin, double peMax, double dtMin, double dtMax,
    int maxEval = 300
);