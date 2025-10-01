#include "Lifetime_Fit.h"

namespace lifefit {

bool get_hold_time(const json& params, const int run, double& ht_out) {
    std::string z = std::to_string(run);
    if (params.contains(z) && params[z].contains("hold_time")) {
        ht_out = params[z]["hold_time"].get<double>();
        return true;
    }
    return false;
}

double ssq_scaled(const std::vector<std::pair<double, double>>& ht_y, double tau) {
    if (tau <= 0.0 || ht_y.size() < 2) return std::numeric_limits<double>::infinity();

    const size_t n = ht_y.size();
    double sum = 0.0, sum2 = 0.0;
    for (const auto& [ht, y] : ht_y) {
        const double s = y * std::exp(ht / tau);
        sum += s;
        sum2 += s * s;
    }
    const double mean = sum / static_cast<double>(n);
    const double var = sum2 - static_cast<double>(n) * mean * mean;
    return (var > 0.0) ? std::sqrt(var) : 0.0;
}

// Bounded Brent minimization for a scalar function (no derivatives)
template <typename F>
static double brent_minimize_bounded(F f, double ax, double bx,
                                     int max_iter = 100, double tol = 1e-6)
{
    // Adapted from classic Brent (as in scipy.optimize._minimize_scalar)
    const double CGOLD = 0.3819660112501051; // 1 - 1/phi
    const double ZEPS  = std::numeric_limits<double>::epsilon() * 1e-3;

    double a = std::min(ax, bx);
    double b = std::max(ax, bx);
    double x = a + CGOLD * (b - a);
    double w = x, v = x;
    double fx = f(x), fw = fx, fv = fx;

    double d = 0.0;        // movement on step before last
    double e = 0.0;        // movement on last step

    for (int iter = 0; iter < max_iter; ++iter) {
        const double xm = 0.5 * (a + b);
        const double tol1 = tol * std::fabs(x) + ZEPS;
        const double tol2 = 2.0 * tol1;

        // Convergence check
        if (std::fabs(x - xm) <= (tol2 - 0.5 * (b - a))) break;

        double p = 0.0, q = 0.0, r = 0.0;
        if (std::fabs(e) > tol1) {
            // parabolic fit
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);
            const double etemp = e;
            e = d;

            // accept parabolic step if in bounds and not too small
            if (std::fabs(p) < std::fabs(0.5 * q * etemp) &&
                p > q * (a - x) && p < q * (b - x))
            {
                d = p / q;
                double u = x + d;
                if (u - a < tol2 || b - u < tol2) d = (x < xm ? tol1 : -tol1);
            } else {
                // golden
                e = (x < xm) ? (b - x) : (a - x);
                d = CGOLD * e;
            }
        } else {
            // golden
            e = (x < xm) ? (b - x) : (a - x);
            d = CGOLD * e;
        }

        const double u = (std::fabs(d) >= tol1) ? (x + d) : (x + (d > 0 ? tol1 : -tol1));
        const double fu = f(u);

        // Update brackets
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return x;
}

double minimize_tau_brent(const std::vector<std::pair<double,double>>& ht_y,
                           double tau_lo, double tau_hi,
                           int max_iter, double rel_tol)
{
    auto f = [&](double tau){ return ssq_scaled(ht_y, tau); };
    return brent_minimize_bounded(f, tau_lo, tau_hi, max_iter, rel_tol);
}

static double jackknife_tau_se(const std::vector<std::pair<double,double>>& ht_y)
{
    const size_t n = ht_y.size();
    if (n < 3) return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> tau_loos; tau_loos.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        std::vector<std::pair<double,double>> sample; sample.reserve(n-1);
        sample.insert(sample.end(), ht_y.begin(), ht_y.begin() + i);
        sample.insert(sample.end(), ht_y.begin() + i + 1, ht_y.end());
        double t = minimize_tau_brent(sample);
        if (std::isfinite(t)) tau_loos.push_back(t);
    }
    if (tau_loos.size() < 2) return std::numeric_limits<double>::quiet_NaN();

    double mean = 0.0; for (double t : tau_loos) mean += t; mean /= tau_loos.size();
    double accum = 0.0; for (double t : tau_loos) { double d = t - mean; accum += d*d; }
    // classic jackknife SE
    return std::sqrt( (tau_loos.size() - 1.0) / tau_loos.size() * accum );
}

static bool pick_y_for_run(const std::map<std::string,double>& segmap,
                           DaggerMode mode,
                           const std::string& seg,
                           const std::vector<std::string>& seg_list_for_overall,
                           double& y_out)
{
    if (mode == DaggerMode::Segment) {
        auto it = segmap.find(seg);
        if (it == segmap.end()) return false;
        y_out = it->second;
        return std::isfinite(y_out) && y_out > 0.0;
    } else { // Overall
        double sum = 0.0;
        bool any = false;
        for (const auto& s : seg_list_for_overall) {
            auto it = segmap.find(s);
            if (it == segmap.end()) continue;
            if (std::isfinite(it->second) && it->second > 0.0) {
                sum += it->second;
                any = true;
            }
        }
        if (!any) return false;
        y_out = sum;
        return true;
    }
}

TauFit fit_tau(const std::map<int, std::map<std::string,double>>& counts_by_run_seg,
               const json& params,
               DaggerMode mode,
               const std::string& seg,
               const std::vector<std::string>& seg_list_for_overall)
{
    std::vector<std::pair<double,double>> ht_y; // (hold_time, y)
    ht_y.reserve(counts_by_run_seg.size());

    for (const auto& [run, segmap] : counts_by_run_seg) {
        double y = 0.0;
        if (!pick_y_for_run(segmap, mode, seg, seg_list_for_overall, y)) continue;

        double ht = 0.0;
        if (!get_hold_time(params, run, ht)) {
            std::cerr << "[fit_tau] Missing hold_time for run " << run << "\n";
            continue;
        }
        ht_y.emplace_back(ht, y);
    }

    TauFit out{std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN(),
               static_cast<int>(ht_y.size())};

    if (ht_y.size() < 3) {
        std::cerr << "[fit_tau] Not enough runs for fit (n=" << ht_y.size() << ")\n";
        return out;
    }

    const double tau_hat = minimize_tau_brent(ht_y);
    const double dtau    = jackknife_tau_se(ht_y);
    out.tau  = tau_hat;
    out.dtau = dtau;
    return out;
}

}