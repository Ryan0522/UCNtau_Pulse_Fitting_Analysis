#pragma once
#include <vector>
#include <utility>

struct FitProblem {
    const std::vector<int>* observed;
    const std::vector<std::vector<double>>* pdfLookup;
    const std::vector<double>* fixedExpected;
    int nTime;
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
    int maxEval = 300
);

struct ModelSelectOptions {
    int closeBins = 10;
    bool useBIC = false;
};

FitResult model_select_1_vs_2(
    const FitProblem& prob,
    const std::vector<std::pair<double, double>>& seeds_sorted_by_dt,
    double peMin, double peMax, double dtMin, double dtMax,
    const ModelSelectOptions& opt
); 