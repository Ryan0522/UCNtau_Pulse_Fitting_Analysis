#pragma once
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <utility>   // std::pair
#include <limits>
#include <cmath>
#include <random>
#include <iostream>
#include <json.hpp>

using json = nlohmann::json;

namespace lifefit{

enum class DaggerMode { Segment, Overall };

struct TauFit {
    double tau;
    double dtau;
    int n_runs;
};

// Fetch hold_time from params for a run (int or string key).
bool get_hold_time(const json& params, const int run, double& ht_out);

// STL-friendly objective: sum of squares of y*exp(ht/tau).
double ssq_scaled(const std::vector<std::pair<double,double>>& ht_y, double tau);

// 1-D optimizer for tau (golden-section).
double minimize_tau_brent(const std::vector<std::pair<double,double>>& ht_y,
                           double tau_lo = 0.1, double tau_hi = 1e5,
                           int max_iter = 200, double rel_tol = 1e-6);

// Build (ht, y_norm) for a segment and fit tau.
TauFit fit_tau(const std::map<int, std::map<std::string,double>>& counts_by_run_seg,
               const json& params,
               DaggerMode mode,
               const std::string& seg = "",
               const std::vector<std::string>& seg_list_for_overall = {"12","34","56","78"});

} // namespace lifefit