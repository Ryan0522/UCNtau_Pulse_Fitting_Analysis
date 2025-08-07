import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def calculate_lifetime_2021(counts, ts):
    errors = np.sqrt(counts)

    def chisqfunc(params):
        model = params[0] * np.exp(-(ts) / params[1])
        chisq = np.sum(((counts - model) / errors)**2 / (6)) # # data point - # dof + 1
        return chisq

    x0 = np.array([1000000, 800])
    result = opt.minimize(chisqfunc, x0, method='Powell')
    taus = np.arange(result.x[1] - 20, result.x[1] + 20, 0.1)
    x0s = [np.array([result.x[0], tau]) for tau in taus]
    points = [chisqfunc(x0) for x0 in x0s]
    min_point = np.min(points)
    tau_left = taus[np.where(points < (min_point + 1))[0][0]]
    tau_right = taus[np.where(points < (min_point + 1))[0][-1]]
    error = (tau_right - tau_left) / (2 * np.sqrt(len(counts)))  
    return result.x[1], error

runinfo_file = 'runinfo_2022_all.csv'
analysis_dir = './output/results/'

runinfo = pd.read_csv(runinfo_file)
run_numbers = runinfo['Run Number'].values
hold_times = runinfo['Holding Time'].values

unique, ht_counts = np.unique(hold_times, return_counts=True)

fillucn_cols = [col for col in runinfo.columns if col.startswith('fillUCN')]

segment_labels = ['12', '34', '56', '78']
holdtime_labels = [20, 50, 100, 200, 1550]
results = {}

for run, hold_t in zip(run_numbers, hold_times):
    analysis_file = os.path.join(analysis_dir, f'PulseAnalysis_{run}.csv')
    if not os.path.exists(analysis_file):
        # print(f"Missing {analysis_file}, skipping.")
        continue

    df = pd.read_csv(analysis_file)
    df.columns = df.columns.str.strip()

    for seg in df['Segment'].unique():
        key = (hold_t, str(seg))
        seg_mask = df['Segment'] == seg
        times = df.loc[seg_mask, 'Time (us)'].values
        PE = df.loc[seg_mask, 'PE'].values
        bg = df.loc[seg_mask, 'Event'].values
        if key not in results:
            results[key] = {'times': [], 'PE': [], 'bg_flag': []}
        results[key]['times'].extend(times)
        results[key]['PE'].extend(PE)
        results[key]['bg_flag'].extend(bg)

thresholds = np.arange(5, 20)
all_lifetimes = {}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=200)
for hold_t in holdtime_labels:
    PE_values = []
    for seg in segment_labels:
        key = (hold_t, seg)
        PE_values.extend(results[key]['PE'])
    ax1.hist(PE_values, bins=np.arange(200), histtype='step', label=f'Hold Time {hold_t}s')
    ax2.hist(PE_values, bins=np.arange(200), histtype='step', label=f'Hold Time {hold_t}s')
ax1.grid()
ax1.legend(loc='best')
ax2.set_yscale('log')
ax2.grid()
ax2.legend(loc='best')
plt.savefig(f'./output/graphs/pe_dist.png')
plt.close()

PMT_Channel = ['12', '34', '1112', '1314']
for i, seg in enumerate(segment_labels):
    lifetimes = []
    dlifetimes = []
    for thresh in thresholds:
        counts = []
        for (hold_t, c) in zip(holdtime_labels, ht_counts):
            key = (hold_t, seg)
            all_PE = np.array(results[key]['PE'])
            all_bg_flag = np.array(results[key]['bg_flag'])
            signal_events = all_PE[(all_PE > thresh) & (all_bg_flag == 1)]
            bg_events = all_PE[(all_PE > thresh) & (all_bg_flag == 0)]
            bg_rate = len(bg_events) / 60.0
            corrected_count = len(signal_events) - bg_rate * 60
            
            fillucn_col = f'fillUCN{PMT_Channel[i]}'
            fill_ucn_value = runinfo[(runinfo['Run Number'] == run)][fillucn_col].values
            normalized_count = corrected_count / fill_ucn_value[0]
            counts.append(normalized_count / c)

        tau, dtau = calculate_lifetime_2021(counts, np.array(holdtime_labels))
        lifetimes.append(tau)
        dlifetimes.append(dtau)

    all_lifetimes[seg] = (lifetimes, dlifetimes)
    plt.errorbar(thresholds, lifetimes, yerr=dlifetimes, fmt='o-', color=f'C{i}', label=f"Segment {seg}")

plt.ylabel("Lifetime (s)")
plt.xlabel("PE Threshold")
plt.grid()
plt.legend(loc='best')
plt.savefig(f'./output/graphs/lifetime_vs_threshold.png')
plt.close()