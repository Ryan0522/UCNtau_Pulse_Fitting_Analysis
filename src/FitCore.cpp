#include "FitCore.h"

#include <cmath>
#include <cfloat>
#include <numeric>
#include <limits>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ================== Global Variables for BG and Pileup ==================
static const std::vector<double>* fixedExpected_ = nullptr;
static double background_per_bin_ = 0.0;
static bool debug_ = false;

// ================== Poisson NLL ==================
double nll_poisson_binned(
    const std::vector<double>& k,
    const std::vector<double>& mu,
    const std::vector<double>* fixedExpected,
    double background_per_bin,
    double eps
) {
    assert(k.size() == mu.size());
    double nll = 0.0;
    for (size_t b = 0; b < k.size(); ++b) {
        double mub = mu[b];
        if (fixedExpected && b < fixedExpected->size()) 
            mub += std::max(0.0, (*fixedExpected)[b]);
        mub += std::max(0.0, background_per_bin);

        if (mub < eps) mub = eps;

        nll += mub - k[b]*std::log(mub);
    }
    return nll;
}

// ================== Expected Histogram (with fractional bins) ==================
static inline void _accumulate_template(
    std::vector<double>& mu,
    const std::vector<double>& tpl,
    double weight,
    double start_frac,
    double bin_width_us
){
    // start_frac = (t0 - window_t0)/bin_width
    const int i0 = (int)std::floor(start_frac);
    const double frac = start_frac - i0;

    const int nBins = (int)mu.size();
    const int L = (int)tpl.size();

    for (int j = 0; j < L; ++j) {
        int left_idx  = i0 + j;
        int right_idx = left_idx + 1;

        // contribution to left_idx
        if (left_idx >= 0 && left_idx < nBins) {
            mu[left_idx] += weight * (1.0 - frac) * tpl[j] * bin_width_us;
        }
        // contribution to right_idx
        if (right_idx >= 0 && right_idx < nBins) {
            mu[right_idx] += weight * frac * tpl[j] * bin_width_us;
        }
    }
}

std::vector<double> expected_hist(
    const std::vector<double>& pmf_unit,
    const std::vector<double>& dts_us,
    const std::vector<double>& pes,
    const std::vector<double>& bin_edges_us
){
    const size_t nbins = bin_edges_us.size() - 1;
    std::vector<double> mu(nbins, 0.0);

    // enforce pmf_unit normalized
    double norm = 0.0;
    for (double v : pmf_unit) norm += v;
    if (norm <= 0.0) return mu;

    const double bin_width_us = bin_edges_us[1] - bin_edges_us[0];
    const double window_t0 = bin_edges_us.front();

    for (size_t i = 0; i < dts_us.size(); ++i) {
        double t0 = dts_us[i];
        double pe = pes[i];
        if (pe <= 0.0) continue;

        double start_frac = (t0 - window_t0) / bin_width_us;
        _accumulate_template(mu, pmf_unit, pe / norm, start_frac, bin_width_us);
    }
    return mu;
}

// ================== Unit Shape for t ==================
std::vector<double> unit_shape_for_t(
    double t_us,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& bin_edges_us,
    double bin_width_us
) {
    const size_t nbins = bin_edges_us.size() - 1;
    std::vector<double> s(nbins, 0.0);

    double norm = 0.0;
    for (double v : pmf_unit) norm += v;
    if (norm <= 0.0) return s;

    double start_frac = (t_us - bin_edges_us.front()) / bin_width_us;
    _accumulate_template(s, pmf_unit, 1.0 / norm, start_frac, bin_width_us);
    return s;
}

// ================== Initialize Amplitude for New Pulse ==================
double init_amp_for_new_pulse(
    const std::vector<double>& k,
    const std::vector<double>& mu_current,
    const std::vector<double>& s_new,
    double amp_lo,
    double amp_hi,
    double eps
) {
    double num = 0.0;
    double den = eps;
    for (size_t b = 0; b < k.size(); ++b) {
        double resid = k[b] - mu_current[b];
        if (resid < 0.0) resid = 0.0;
        num += s_new[b] * resid;
        den += s_new[b] * s_new[b];
    }
    double a = (den > 0.0 ? num/den : 0.0);
    if (a < amp_lo) a = amp_lo;
    if (a > amp_hi) a = amp_hi;
    return a;
}

// ================== build S from times (cache) ==================
std::vector<std::vector<double>> build_S_from_times(
    const std::vector<double>& times_us,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& bin_edges_us,
    double bin_width_us
){
    const size_t nbins = bin_edges_us.size() - 1;
    const size_t J = times_us.size();

    std::vector<std::vector<double>> S;
    S.reserve(J);

    for (size_t j = 0; j < J; ++j) {
        S.push_back( unit_shape_for_t(times_us[j], pmf_unit, bin_edges_us, bin_width_us) );
    }

    std::vector<std::vector<double>> S_out;
    S_out.assign(nbins, std::vector<double>(J, 0.0));
    for (size_t j = 0; j < J; ++j) {
        for (size_t b = 0; b < nbins; ++b) {
            S_out[b][j] = S[j][b];
        }
    }
    return S_out;
}

// ================== S @ a ==================
std::vector<double> S_times_a(
    const std::vector<std::vector<double>>& S,
    const std::vector<double>& amps
){
    const size_t nbins = S.size();
    const size_t J = amps.size();
    std::vector<double> mu(nbins, 0.0);
    for (size_t b = 0; b < nbins; ++b) {
        double sum_b = 0.0;
        for (size_t j = 0; j < J; ++j) {
            sum_b += S[b][j] * amps[j];
        }
        mu[b] = sum_b;
    }
    return mu;
}

// ================== refit all amplitudes lbfgsb ==================
std::vector<double> refit_all_amplitudes_lbfgsb(
    const std::vector<double>& k,
    const std::vector<std::vector<double>>& S,
    const std::vector<double>& a_init,
    double amp_lo,
    double amp_hi,
    int max_iter,
    double ftol,
    double eps
) {
    const size_t J = a_init.size();
    const size_t nbins = k.size();
    std::vector<double> a = a_init;

    auto clamp_amps = [&](std::vector<double>& arr) {
        for (double& v : arr) {
            if (v < amp_lo) v = amp_lo;
            if (v > amp_hi) v = amp_hi;
        }
    };
    clamp_amps(a);

    double prev_f = std::numeric_limits<double>::infinity();

    for (int it = 0; it < max_iter; ++it) {
        // mu = S @ a
        std::vector<double> mu(nbins, 0.0);
        for (size_t b = 0; b < nbins; ++b) {
            double sum_b = 0.0;
            for (size_t j = 0; j < J; ++j) sum_b += S[b][j]*a[j];
            if (sum_b < eps) sum_b = eps;
            mu[b] = sum_b;
        }

        // f = sum ( mu - k log mu )
        double f = 0.0;
        for (size_t b = 0; b < nbins; ++b) {
            f += mu[b] - k[b]*std::log(mu[b]);
        }

        // check improvement
        if (std::fabs(prev_f - f) < ftol) break;
        prev_f = f;

        // grad g_j = d/d a_j f = sum_b S[b][j]*(1 - k[b]/mu[b])
        std::vector<double> g(J, 0.0);
        for (size_t j = 0; j < J; ++j) {
            double gj = 0.0;
            for (size_t b = 0; b < nbins; ++b) {
                double r = 1.0 - k[b]/mu[b];
                gj += S[b][j] * r;
            }
            g[j] = gj;
        }

        // simple step: a_new = a - alpha*g, with crude line search
        double alpha = 1e-2;
        bool accepted = false;
        for (int ls = 0; ls < 10; ++ls) {
            std::vector<double> a_try(J);
            for (size_t j = 0; j < J; ++j) {
                a_try[j] = a[j] - alpha*g[j];
            }
            clamp_amps(a_try);

            // evaluate
            std::vector<double> mu_try(nbins, 0.0);
            for (size_t b = 0; b < nbins; ++b) {
                double sum_b = 0.0;
                for (size_t j = 0; j < J; ++j) sum_b += S[b][j]*a_try[j];
                if (sum_b < eps) sum_b = eps;
                mu_try[b] = sum_b;
            }

            double f_try = 0.0;
            for (size_t b = 0; b < nbins; ++b) {
                f_try += mu_try[b] - k[b]*std::log(mu_try[b]);
            }

            if (f_try <= f) {
                a = a_try;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            // couldn't improve
            break;
        }
    }

    clamp_amps(a);
    return a;
}

// ================== Re-balance Amplitudes for Close Pairs ==================
std::vector<double> rebalance_close_pairs(
    const std::vector<double>& times_us,
    const std::vector<double>& amps_in,
    double d_merge_us,
    double amp_lo,
    double amp_hi
) {
    if (times_us.empty()) return amps_in;

    std::vector<int> order(times_us.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(),
        [&](int a, int b){return times_us[a] < times_us[b];});
    
    std::vector<std::vector<int>> groups;
    std::vector<int> cur;
    cur.push_back(order[0]);
    for (size_t ii = 1; ii < order.size(); ++ii) {
        int prev = order[ii - 1];
        int now = order[ii];
        if (std::fabs(times_us[now] - times_us[prev]) < d_merge_us) {
            cur.push_back(now);
        } else {
            groups.push_back(cur);
            cur.clear();
            cur.push_back(now);
        }
    }
    groups.push_back(cur);

    std::vector<double> amps_out = amps_in;
    for (auto& g : groups) {
        if (g.size() >= 2) {
            double mean_a = 0.0;
            for (int idx : g) mean_a += amps_out[idx];
            mean_a /= g.size();
            if (mean_a < amp_lo) mean_a = amp_lo;
            if (mean_a > amp_hi) mean_a = amp_hi;
            for (int idx : g) amps_out[idx] = mean_a;
        }
    }
    return amps_out;
}

// utility
static bool valid_time_candidate(
    double t,
    const std::vector<double>& accepted_times,
    double dmin_us
) {
    for (double tj : accepted_times) {
        if (std::fabs(t - tj) < dmin_us) return false;
    }
    return true;
}

// ================== Best T_a for cluter ==================
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
    double eps
) {
    AddProposal best;
    best.t_us = 0.0;
    best.delta_nll = 0.0;
    best.valid = false;

    // baseline model
    std::vector<double> amps0 = amps_curr;
    std::vector<std::vector<double>> S0 = S_curr;
    std::vector<double> mu0;
    if (!S0.empty() && !amps0.empty()) {
        mu0 = S_times_a(S0, amps0);
    } else {
        mu0.assign(k.size(), 0.0);
    }
    double nll0 = nll_poisson_binned(k, mu0, fixedExpected_, background_per_bin_, eps);

    if (debug_) {
        std::cerr << "[TA]  scanning t ∈ [" << t_left << "," << t_right
                  << "] step=" << t_step_us
                  << " with " << times_curr.size()
                  << " existing pulses\n";
    }

    for (double t = t_left; t <= t_right + 0.5 * t_step_us; t += t_step_us) {
        if (!valid_time_candidate(t, times_curr, dmin_us)) continue;
        std::vector<double> s_new = unit_shape_for_t(t, pmf_unit, bin_edges_us, bin_width_us);
        double a_init_new = init_amp_for_new_pulse(k, mu0, s_new, amp_lo, amp_hi, eps);

        std::vector<double> times1 = times_curr;
        std::vector<double> amps1 = amps_curr;
        times1.push_back(t);
        amps1.push_back(a_init_new);

        amps1 = rebalance_close_pairs(times1, amps1, 5.0, amp_lo, amp_hi); // 5.0 us to rebalance
        std::vector<std::vector<double>> S1 = S0;
        {
            const size_t nbins = s_new.size();
            const size_t Jold = S0.empty() ? 0 : S0[0].size();
            const size_t Jnew = Jold + 1;
            if (S1.empty()) {
                S1.assign(nbins, std::vector<double>(Jnew, 0.0));
                for (size_t b = 0; b < nbins; ++b) {
                    S1[b][Jold] = s_new[b];
                }
            } else {
                for (size_t b = 0; b < nbins; ++b) {
                    S1[b].push_back(s_new[b]);
                }
            }
        }
        std::vector<double> a1 = refit_all_amplitudes_lbfgsb(
            k, S1, amps1, amp_lo, amp_hi
        );

        std::vector<double> mu1 = S_times_a(S1, a1);
        double nll1 = nll_poisson_binned(k, mu1, fixedExpected_, background_per_bin_, eps);
        double d = nll0 - nll1;

        if (debug_) {
            std::cerr << "    [TA] t=" << std::fixed << std::setprecision(2) << t
                      << " µs  a_init=" << a_init_new
                      << "  nll0=" << nll0
                      << "  nll1=" << nll1
                      << "  ΔNLL=" << d << "\n";
            // Print the full candidate times and amplitudes after refit
            std::cerr << "         times1=[";
            for (size_t i = 0; i < times1.size(); ++i) {
                if (i) std::cerr << ", ";
                std::cerr << std::fixed << std::setprecision(3) << times1[i];
                if (i + 1 == times1.size()) std::cerr << "*";  // mark the newly added one
            }
            std::cerr << "]  a1=[";
            for (size_t i = 0; i < a1.size(); ++i) {
                if (i) std::cerr << ", ";
                std::cerr << std::fixed << std::setprecision(3) << a1[i];
                if (i + 1 == a1.size()) std::cerr << "*";  // mark the newly added one
            }
            std::cerr << "]\n";
        }

        if (!best.valid || d > best.delta_nll) {
            best.valid = true;
            best.t_us = t;
            best.delta_nll = d;
            best.times_out = times1;
            best.amps_out = a1;
            best.mu_out = mu1;
            best.S_out = S1;
        }
    }

    if (debug_) {
        if (best.valid)
            std::cerr << "  [TA]  best_t=" << best.t_us
                      << "  ΔNLL_best=" << best.delta_nll
                      << "  amp_best=" << (best.amps_out.empty() ? 0.0 : best.amps_out.back()) << "\n";
        else
            std::cerr << "  [TA]  no valid pulse found in range\n";
    }
    
    return best;
}

// ================== Backward Check for Greedy ==================
void backward_prune(
    const std::vector<double>& k,
    std::vector<std::vector<double>>& S,
    std::vector<double>& times_us,
    std::vector<double>& amps,
    double delta_null_cut,
    double amp_lo,
    double amp_hi,
    double eps
){
    bool changed = true;
    while (changed && (int)times_us.size() > 1) {
        changed = false;

        // baseline
        std::vector<double> mu_base = S_times_a(S, amps);
        double nll_base = nll_poisson_binned(k, mu_base, fixedExpected_, background_per_bin_, eps);

        for (int j = (int)times_us.size() - 1; j >= 0; --j) {
            // try removing pulse j
            std::vector<double> times_try(times_us);
            times_try.erase(times_try.begin() + j);

            // drop column j from S
            std::vector<std::vector<double>> S_try;
            if (!times_try.empty()) {
                const size_t nbins = S.size();
                const size_t Jold  = S[0].size();
                const size_t Jnew  = Jold - 1;
                S_try.assign(nbins, std::vector<double>(Jnew, 0.0));
                for (size_t b = 0; b < nbins; ++b) {
                    size_t col2 = 0;
                    for (size_t col = 0; col < Jold; ++col) {
                        if ((int)col == j) continue;
                        S_try[b][col2++] = S[b][col];
                    }
                }
            }

            std::vector<double> amps_try;
            if (!times_try.empty()) {
                amps_try.reserve(times_try.size());
                for (size_t jj = 0; jj < amps.size(); ++jj) {
                    if ((int)jj == j) continue;
                    amps_try.push_back(amps[jj]);
                }
                amps_try = refit_all_amplitudes_lbfgsb(
                    k, S_try, amps_try, amp_lo, amp_hi
                );
            }

            std::vector<double> mu_try;
            if (!times_try.empty()) {
                mu_try = S_times_a(S_try, amps_try);
            } else {
                mu_try.assign(k.size(), 0.0);
            }

            double nll_try = nll_poisson_binned(k, mu_try, fixedExpected_, background_per_bin_, eps);
            double loss = nll_try - nll_base;
            if (debug_) {
                std::cerr << "[PRUNE] pulse#" << j
                    << " loss=" << loss
                    << " threshold=" << delta_null_cut
                    << (loss < delta_null_cut ? " → remove\n" : " → keep\n");
            }

            if (loss < delta_null_cut) {
                // remove it permanently
                times_us = times_try;
                amps     = amps_try;
                S        = S_try;
                changed  = true;
                break;
            }
        }
    }
}

// ================== Cluster Seed DB-scan ==================
std::vector<std::vector<double>> clusterSeedsDBSCAN1D(
    const std::vector<double>& seed_times_us,
    double eps_us,
    int min_samples
) {
    std::vector<double> sorted = seed_times_us;
    std::sort(sorted.begin(), sorted.end());

    std::vector<std::vector<double>> clusters;
    std::vector<double> cur;
    for (size_t i = 0; i < sorted.size(); ++i) {
        if (cur.empty()) {
            cur.push_back(sorted[i]);
        } else {
            if ( std::fabs(sorted[i] - cur.back()) <= eps_us ) {
                cur.push_back(sorted[i]);
            } else {
                if ((int)cur.size() >= min_samples) {
                    clusters.push_back(cur);
                }
                cur.clear();
                cur.push_back(sorted[i]);
            }
        }
    }
    if (!cur.empty() && (int)cur.size() >= min_samples) {
        clusters.push_back(cur);
    }
    return clusters;
}

// ================== Building CLuster Bounds ==================
std::vector<ClusterBound> buildClusterBounds(
    const std::vector<std::vector<double>>& clusters,
    double bin_left_us,
    double bin_right_us,
    double pulse_shape_shift_us,
    double max_offset_us
) {
    std::vector<double> centers;
    centers.reserve(clusters.size());
    for (auto const& c : clusters) {
        double tmin = c.empty() ? 0.0 : c.front();
        centers.push_back(tmin);
    }

    std::vector<int> order(centers.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(),
        [&](int a, int b){ return centers[a] < centers[b]; });

    std::vector<ClusterBound> out;
    out.resize(order.size());

    if (debug_) {
        std::cerr << "\n[CB] Building cluster bounds"
                  << " | n_clusters=" << clusters.size()
                  << " | bin_left=" << bin_left_us
                  << " | bin_right=" << bin_right_us
                  << " | shift=" << pulse_shape_shift_us
                  << " | max_offset=" << max_offset_us << "\n";
    }

    for (size_t ii = 0; ii < order.size(); ++ii) {
        int idx = order[ii];

        double t_center = centers[idx];

        double left_neighbor  = (ii > 0)               ? centers[ order[ii-1] ] : -std::numeric_limits<double>::infinity();
        double right_neighbor = (ii+1 < order.size())  ? centers[ order[ii+1] ] :  std::numeric_limits<double>::infinity();

        double half_gap_left  = std::isfinite(left_neighbor)
                              ? 0.5 * (t_center - left_neighbor)
                              : max_offset_us;
        double half_gap_right = std::isfinite(right_neighbor)
                              ? 0.5 * (right_neighbor - t_center)
                              : max_offset_us;

        double left_span  = std::min(max_offset_us,  half_gap_left);
        double right_span = std::min(max_offset_us, half_gap_right);

        double center_for_search = t_center - pulse_shape_shift_us;

        double tL = std::max(bin_left_us,  center_for_search - left_span);
        double tR = std::min(bin_right_us, center_for_search + right_span);

        ClusterBound cb;
        cb.t_center_us = t_center;
        cb.tL_us       = tL;
        cb.tR_us       = tR;
        out[ii] = cb;

        if (debug_) {
            std::cerr << std::fixed << std::setprecision(3)
                      << "  [CB] cluster#" << ii
                      << " center=" << t_center
                      << " | neighbors=("
                      << (std::isfinite(left_neighbor)  ? std::to_string(left_neighbor)  : "-inf") << ","
                      << (std::isfinite(right_neighbor) ? std::to_string(right_neighbor) : "inf") << ")"
                      << " | half_gaps=(" << half_gap_left << "," << half_gap_right << ")"
                      << " | span=(" << left_span << "," << right_span << ")"
                      << " | shifted_center=" << center_for_search
                      << " | bound=[" << tL << "," << tR << "]\n";
        }
    }
    return out;
}

// ================== Run Greedy LRT ==================
FitResult runGreedyLRT(
    const std::vector<double>& k_hist,
    const std::vector<double>& bin_edges_us,
    const std::vector<double>& pmf_unit,
    const std::vector<double>& seed_times_us,
    double bin_width_us,
    double max_offset_us,
    double delta_null_cut,
    double dmin_us,
    double amp_min_cut_pe,
    const std::vector<double>* fixedExpected,
    double background_per_bin,
    bool debug
) {
    fixedExpected_ = fixedExpected;
    background_per_bin_ = background_per_bin;
    debug_ = debug;

    FitResult result;
    result.final_nll = 0.0;

    const double bin_left_us = bin_edges_us.front();
    const double bin_right_us = bin_edges_us.back();

    std::vector<std::vector<double>> clusters = clusterSeedsDBSCAN1D(seed_times_us, /*eps_us=*/1.0, /*min_samples=*/1);
    std::vector<ClusterBound> clusterBounds = buildClusterBounds(clusters, bin_left_us, bin_right_us, /*pulse_shape_shift_us=*/5.0, max_offset_us);

    if (debug_) {
        std::cerr << "\n[LRT] Starting Greedy Fit"
                  << " | nbins=" << k_hist.size()
                  << " | seeds=" << seed_times_us.size()
                  << " | clusters=" << clusters.size() << "\n";

        for (size_t i = 0; i < clusters.size(); ++i) {
            const auto& c = clusters[i];
            std::cerr << "  cluster#" << i << " size=" << c.size() << " times(us): ";
            for (double t : c) std::cerr << std::fixed << std::setprecision(2) << t << " ";
            std::cerr << "\n";
            std::cerr << "    -> bound [" << clusterBounds[i].tL_us
                      << "," << clusterBounds[i].tR_us << "]\n";
        }
    }

    std::vector<double> curr_times_us;
    std::vector<double> curr_amps;
    std::vector<std::vector<double>> S_curr;

    std::vector<int> remaining(clusterBounds.size());
    for (size_t i = 0; i < remaining.size(); ++i) remaining[i] = (int)i;

    int iter = 0;
    while (true) {
        ++iter;
        double best_delta = 0.0;
        int best_idx = -1;
        AddProposal best_prop;

        for (int ci : remaining) {
            const ClusterBound& B = clusterBounds[ci];
            AddProposal prop = best_ta_for_cluster(
                k_hist, S_curr, curr_times_us, curr_amps,
                pmf_unit, bin_edges_us, bin_width_us,
                B.tL_us, B.tR_us, bin_width_us,
                /*amp_lo=*/0.0, /*amp_hi=*/300.0, dmin_us
            );
            if (prop.valid && prop.delta_nll > best_delta) {
                best_delta = prop.delta_nll;
                best_idx = ci;
                best_prop = prop;
            }
        }

        if (best_idx < 0 || best_delta < delta_null_cut) {
            if (debug_)
                std::cerr << "[LRT]  iteration#" << iter
                          << "  no ΔNLL above threshold (" << delta_null_cut << ") → stop\n";
            break;
        }

        if (!best_prop.amps_out.empty()) {
            double newest_amp = best_prop.amps_out.back();
            if (newest_amp < amp_min_cut_pe) {
                
                if (debug_) {
                    std::cerr << "[LRT]  iteration#" << iter
                            << "  rejected pulse amp=" << newest_amp
                            << " < " << amp_min_cut_pe << " PE\n";
                }
                
                remaining.erase(std::remove(remaining.begin(), remaining.end(), best_idx), remaining.end());
                
                if (debug_) {
                    std::cerr << "[LRT]  remaining clusters after removing #" << best_idx << ": ";
                    if (remaining.empty()) {
                        std::cerr << "(none)";
                    } else {
                        for (int ri : remaining) std::cerr << ri << " ";
                    }
                    std::cerr << "\n";
                }

                continue;
            }
        }

        curr_times_us = best_prop.times_out;
        curr_amps = best_prop.amps_out;
        S_curr = best_prop.S_out;

        if (debug_) {
            std::cerr << "[LRT]  iteration#" << iter
                      << "  accepted cluster#" << best_idx
                      << "  ΔNLL=" << std::fixed << std::setprecision(3) << best_delta
                      << "  pulses_now=" << curr_times_us.size() << "\n";
            for (size_t i = 0; i < curr_times_us.size(); ++i)
                std::cerr << "      (" << curr_times_us[i]
                          << " µs, " << curr_amps[i] << " PE)\n";
        }

        remaining.erase(std::remove(remaining.begin(), remaining.end(), best_idx), remaining.end());

        if (debug_) {
            std::cerr << "[LRT]  remaining clusters after removing #" << best_idx << ": ";
            if (remaining.empty()) {
                std::cerr << "(none)";
            } else {
                for (int ri : remaining) std::cerr << ri << " ";
            }
            std::cerr << "\n";
        }

        if (!curr_times_us.empty()) {
            backward_prune(
                k_hist, S_curr, curr_times_us, curr_amps,
                delta_null_cut, /*amp_lo=*/0.0, /*amp_hi=*/300.0
            );
            if (debug_) {
                std::cerr << "[LRT]  after backward_prune:"
                          << " pulses_now=" << curr_times_us.size();
                if (!curr_times_us.empty()) {
                    std::cerr << " -> ";
                    for (size_t i = 0; i < curr_times_us.size(); ++i)
                        std::cerr << "(" << curr_times_us[i]
                                  << " µs," << curr_amps[i] << " PE) ";
                }
                std::cerr << "\n";
            }
        }
    }

    std::vector<double> mu_final;
    if (!curr_times_us.empty()) {
        mu_final = S_times_a(S_curr, curr_amps);
    } else {
        mu_final.assign(k_hist.size(), 0.0);
    }
    double nll_final = nll_poisson_binned(k_hist, mu_final, fixedExpected_, background_per_bin_);

    result.mu_total = std::move(mu_final);
    result.final_nll = nll_final;
    
    result.pulses.reserve(curr_times_us.size());
    for (size_t i = 0; i < curr_times_us.size(); ++i) {
        PulseCandidate pc;
        pc.t_us   = curr_times_us[i];
        pc.amp_pe = curr_amps[i];
        result.pulses.push_back(pc);
    }

    if (debug_) {
        std::cerr << "[LRT] Done."
                  << "  total_pulses=" << result.pulses.size()
                  << "  final_NLL=" << std::fixed << std::setprecision(4)
                  << result.final_nll << "\n";
    }

    if (debug_) {
        // ---- Build expected including fixedExpected + background ----
        std::vector<double> expected = result.mu_total;
        if (fixedExpected_) {
            for (size_t i = 0; i < expected.size() && i < fixedExpected_->size(); ++i)
                expected[i] += std::max(0.0, (*fixedExpected_)[i]);
        }
        if (background_per_bin_ > 0.0) {
            for (double& e : expected) e += background_per_bin_;
        }

        // ---- Compute simple stats ----
        double sum_obs = std::accumulate(k_hist.begin(), k_hist.end(), 0.0);
        double sum_exp = std::accumulate(expected.begin(), expected.end(), 0.0);
        double max_obs = *std::max_element(k_hist.begin(), k_hist.end());
        double max_exp = *std::max_element(expected.begin(), expected.end());

        std::cerr << "\n[DBG] Histogram sanity check:"
                << " nbins=" << k_hist.size()
                << " sum_obs=" << sum_obs
                << " sum_exp=" << sum_exp
                << " max_obs=" << max_obs
                << " max_exp=" << max_exp
                << " background_per_bin=" << background_per_bin_
                << "\n";

        // ---- Print compressed view (not every bin) ----
        const size_t step = std::max<size_t>(1, k_hist.size() / 20);
        std::cerr << "  bin   obs   exp\n";
        for (size_t i = 0; i < k_hist.size(); i += step) {
            std::cerr << "  "
                    << std::setw(3) << i << "   "
                    << std::setw(6) << k_hist[i] << "   "
                    << std::setw(8) << std::fixed << std::setprecision(2)
                    << expected[i] << "\n";
        }

        // ---- Optional: print a small ASCII overlay ----
        std::cerr << "\n  [ASCII overlay: '*'=obs, '+'=expected]\n";
        const int scale_height = 20;
        double scale = max_obs > 0 ? (scale_height / max_obs) : 1.0;
        for (size_t i = 0; i < k_hist.size(); ++i) {
            int h_obs = (int)std::round(k_hist[i] * scale);
            int h_exp = (int)std::round(expected[i] * scale);
            std::cerr << std::setw(3) << i << " | ";
            for (int h = 0; h < scale_height; ++h) {
                if (h == h_obs && h == h_exp)
                    std::cerr << '#'; // both match
                else if (h == h_obs)
                    std::cerr << '*';
                else if (h == h_exp)
                    std::cerr << '+';
                else
                    std::cerr << ' ';
            }
            std::cerr << "\n";
        }
    }

    return result;
}