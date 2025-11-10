#pragma once
#include <vector>
#include <cstddef>
#include <limits>
#include <cmath>
#include <algorithm>

struct PulseCandidate {
    double t_us;
    double amp_pe;
};

struct FitResult {
    std::vector<PulseCandidate> pulses;
    std::vector<double> mu_total;
    double final_nll;
};

FitResult runGreedyLRT(
    const std::vector<double>& k_hist,            // observed binned counts (same length as bin_edges-1)
    const std::vector<double>& bin_edges_us,      // bin edges [µs], length = nbins+1
    const std::vector<double>& pmf_unit,          // unit pulse PMF shape (area-normalized template at 1 PE total)
    const std::vector<double>& seed_times_us,     // candidate pulse seed times [µs since window start]
    double bin_width_us         = 0.5,            // grid step when scanning around each cluster (=bin_width_us)
    double max_offset_us        = 10.0,           // ± search half-width cap for each cluster
    double delta_null_cut       = 10,             // LRT threshold 
    double dmin_us              = 2.0,            // minimum time spacing between pulses
    double amp_min_cut_pe       = 5.0,            // reject pulses with fitted amp below this
    const std::vector<double>* fixedExpected = nullptr,
    double background_per_bin = 0.0,
    bool debug = false
);

// ================== helper functions ==================
double nll_poisson_binned(
    const std::vector<double>& k,
    const std::vector<double>& mu,
    const std::vector<double>* fixedExpected = nullptr,
    double background_per_bin = 0.0,
    double eps = 1e-12
);

std::vector<double> expected_hist(
    const std::vector<double>& pmf_unit,
    const std::vector<double>& dts_us,
    const std::vector<double>& pes,
    const std::vector<double>& bin_edges_us
);

std::vector<double> unit_shape_for_t(
    double t_us,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& bin_edges_us,
    double bin_width_us
);

double init_amp_for_new_pulse(
    const std::vector<double>& k,
    const std::vector<double>& mu_current,
    const std::vector<double>& s_new,
    double amp_lo,
    double amp_hi,
    double eps = 1e-12
);

std::vector<std::vector<double>> build_S_from_times(
    const std::vector<double>& times_us,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& bin_edges_us,
    double bin_width_us
);

std::vector<double> S_times_a(
    const std::vector<std::vector<double>>& S,
    const std::vector<double>& amps
);

std::vector<double> refit_all_amplitudes_lbfgsb(
    const std::vector<double>& k,
    const std::vector<std::vector<double>>& S,
    const std::vector<double>& a_init,
    double amp_lo,
    double amp_hi,
    int max_iter = 20,
    double ftol = 1e-10,
    double eps = 1e-12
);

std::vector<double> rebalance_close_pairs(
    const std::vector<double>& times_us,
    const std::vector<double>& amps_in,
    double d_merge_us,
    double amp_lo,
    double amp_hi
);

struct AddProposal {
    double t_us;
    double delta_nll;
    std::vector<double> times_out;
    std::vector<double> amps_out;
    std::vector<double> mu_out;
    std::vector<std::vector<double>> S_out;
    bool valid;
};

AddProposal best_ta_for_cluster(
    const std::vector<double>& k,
    const std::vector<std::vector<double>>& S_curr,
    const std::vector<double>& times_curr,
    const std::vector<double>& amps_curr,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& bin_edges_us,
    double bin_width_us,
    double t_left,
    double t_right,
    double t_step_us,
    double amp_lo,
    double amp_hi,
    double dmin_us,
    double eps = 1e-12
);

void backward_prune(
    const std::vector<double>& k,
    std::vector<std::vector<double>>& S,
    std::vector<double>& times_us,
    std::vector<double>& amps,
    double delta_null_cut,
    double amp_lo,
    double amp_hi,
    double eps = 1e-12
);

std::vector<std::vector<double>> clusterSeedsDBSCAN1D(
    const std::vector<double>& seed_times_us,
    double eps_us,
    int min_samples
);

struct ClusterBound {
    double t_center_us;
    double tL_us;
    double tR_us;
};

std::vector<ClusterBound> buildClusterBounds(
    const std::vector<std::vector<double>>& clusters,
    double bin_left_us,
    double bin_right_us,
    double pulse_shape_shift_us,
    double max_offset_us
);