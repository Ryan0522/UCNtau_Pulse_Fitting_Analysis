import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

RUNINFO_CSV = 'runinfo_2022_all.csv'
ANALYSIS_DIR = './output/results/'
GRAPH_DIR = './output/graphs/'
GOOD_RUNS_TXT = './config/2022runlist.txt'

SEGMENTS = ['12', '34', '56', '78']

SEG_TO_FILL = {
    '12': 'fillUCN12',
    '34': 'fillUCN34',
    '56': 'fillUCN1112',
    '78': 'fillUCN1314'
}

def load_good_runs(path):
    good = set()
    if not os.path.exists(path):
        print(f"Warning: Good runs file {path} does not exist.")
        return good
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                good.add(int(s))
            except ValueError:
                pass
    return good

GOOD_RUNS = load_good_runs(GOOD_RUNS_TXT)

# -----------------------------
# Lifetime fitting (profile A)
# -----------------------------
def fit_tau_profiled(counts_raw, ts):
    """
    Fit tau for y_i vs t_i with model A * exp(-t_i/tau),
    profiling A analytically for each tau.

    counts_raw: array-like of background-corrected, normalized counts (>=0 preferred)
    ts:         array-like of hold times (seconds)
    Returns (tau, dtau)
    """
    counts = np.asarray(counts_raw, dtype=float)
    ts = np.asarray(ts, dtype=float)
    errors = np.sqrt(np.maximum(counts, 1e-9))

    w = 1.0 / (errors**2)

    def chi2_for_tau(tau):
        if tau <= 0:
            return np.inf, np.nan
        m = np.exp(-ts / tau)
        num = np.sum(w * counts * m)
        den = np.sum(w * m * m)
        if den <= 0:
            return np.inf, np.nan
        Ahat = num / den
        chi2 = np.sum(w * (counts - Ahat*m)**2)
        return chi2, Ahat

    res = opt.minimize(lambda x: chi2_for_tau(x[0])[0],
                       x0=np.array([800.0]), method='Powell',
                       bounds = [(1.0, 1e6)])
    tau_best = float(res.x[0])
    chi2_min, _ = chi2_for_tau(tau_best)
    
    # Scan around best to get Δχ²=1 band
    span = max(50.0, 0.15 * tau_best)
    taus = np.linspace(max(1.0, tau_best - span), tau_best + span, 2001)
    chi2_vals = np.array([chi2_for_tau(tau)[0] for tau in taus])
    ok = chi2_vals < (np.nanmin(chi2_vals) + 1.0)
    if not np.any(ok):
        return tau_best, np.nan

    left_idx  = np.argmax(ok)
    right_idx = len(ok) - np.argmax(ok[::-1]) - 1
    tau_left, tau_right = taus[left_idx], taus[right_idx]
    dtau = 0.5 * (tau_right - tau_left)

    return tau_best, dtau

runinfo = pd.read_csv(RUNINFO_CSV)
hold_times = np.sort(runinfo['Holding Time'].unique())

results = {}
fillucn_sum = {}

for _, r in runinfo.iterrows():
    run = int(r['Run Number'])
    if GOOD_RUNS and run not in GOOD_RUNS:
        continue

    hold_t = int(r['Holding Time'])
    per_seg_fill = {
        seg: float(r.get(SEG_TO_FILL[seg], 0.0))
        for seg in SEGMENTS
    }

    f = os.path.join(ANALYSIS_DIR, f'PulseAnalysis_{run}.csv')
    if not os.path.exists(f):
        print(f"Missing {f}, skipping run {run}.")
        continue
    
    df = pd.read_csv(f)
    df.columns = df.columns.str.strip()

    for seg in df['Segment'].unique():
        seg = str(seg).strip()
        if seg not in SEGMENTS:
            print(f"Unknown segment {seg} in run {run}, skipping.")
            continue
        
        key = (hold_t, seg)
        if key not in results:
            results[key] = {'times': [], 'PE': [], 'bg_flag': []}
            fillucn_sum[key] = 0.0
        
        mask = df['Segment'] == seg
        results[key]['times'].extend(df.loc[mask, 'Time (us)'].values.tolist())
        results[key]['PE'].extend(df.loc[mask, 'PE'].values.tolist())
        results[key]['bg_flag'].extend(df.loc[mask, 'Event'].values.tolist())
        
        fillucn_sum[key] += per_seg_fill[seg]

plt.figure(figsize=(12, 5), dpi=160)
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

for hold_t in hold_times:
    all_pe = []
    for seg in SEGMENTS:
        key = (hold_t, seg)
        if key not in results:
            print(f"Missing results for hold time {hold_t}s, segment {seg}, skipping.")
            continue
        all_pe.extend(results[key]['PE'])
    if len(all_pe) == 0:
        continue
    ax1.hist(all_pe, bins=np.arange(0, 200, 1), histtype='step', label=f'Hold {hold_t}s')
    ax2.hist(all_pe, bins=np.arange(0, 200, 1), histtype='step', label=f'Hold {hold_t}s')
ax1.set_title('PE distribution (linear)')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax2.set_title('PE distribution (log)')
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3)
ax2.legend()
plt.tight_layout()
plt.savefig(os.path.join(GRAPH_DIR, 'pe_dist.png'))
plt.close()

thresholds = np.arange(5, 20)
lifetime_by_seg = {}

plt.figure(figsize=(7, 5), dpi=160)
for i, seg in enumerate(SEGMENTS):
    lifetimes = []
    dlifetimes = []
    for thresh in thresholds:
        counts_per_hold = []

        for hold_t in hold_times:
            key = (hold_t, seg)
            if key not in results:
                print(f"Missing results for hold time {hold_t}s, segment {seg}, skipping.")
                continue
            pe = np.array(results[key]['PE'])
            bg_flag = np.array(results[key]['bg_flag'])

            sig_n = np.sum((pe > thresh) & (bg_flag == 1))
            bg_n = np.sum((pe > thresh) & (bg_flag == 0))

            corrected = sig_n - bg_n
            
            denom = max(fillucn_sum.get(key, 0.0), 1.0)
            norm_count = corrected / denom

            counts_per_hold.append((hold_t, norm_count))

        if len(counts_per_hold) < 2:
            lifetimes.append(np.nan)
            dlifetimes.append(np.nan)
            continue

        counts_per_hold.sort(key=lambda x: x[0])
        ts = np.array([x[0] for x in counts_per_hold], dtype=float)
        ys = np.array([x[1] for x in counts_per_hold], dtype=float)

        tau, dtau = fit_tau_profiled(ys, ts)
        lifetimes.append(tau)
        dlifetimes.append(dtau)

    lifetimes_by_seg[seg] = (lifetimes, dlifetimes)
    plt.errorbar(thresholds, lifetimes, yerr=dlifetimes, fmt='o-', label=f'Segment {seg}')

plt.xlabel('PE Threshold')
plt.ylabel('Lifetime τ (s)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(GRAPH_DIR, 'lifetime_vs_threshold.png'))
plt.close()