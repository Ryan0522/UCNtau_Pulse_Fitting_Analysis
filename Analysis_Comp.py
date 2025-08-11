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
    # load whitelist of run numbers (one per line); empty set => use all runs
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

runinfo = pd.read_csv(RUNINFO_CSV)
hold_times = np.sort(runinfo['Holding Time'].unique()) # all unique hold times in dataset

results = {}
fillucn_sum = {}

for _, r in runinfo.iterrows(): # aggregate all runs passing filter
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
        if str(seg) not in SEGMENTS:
            print(f"Unknown segment {seg} in run {run}, skipping.")
            continue
        
        key = (hold_t, seg)
        if key not in results:
            results[key] = {'times': [], 'PE': [], 'bg_flag': []}
            fillucn_sum[key] = 0.0
        
        mask = (df['Segment'].tolist() == seg)
        # if run == '26594' or run == 26594:
        #     print(df.head)
        #     print(seg)
        #     print(df.loc[mask])
        results[key]['times'].extend(df.loc[mask, 'Time (us)'].values)
        results[key]['PE'].extend(df.loc[mask, 'PE'].values)
        results[key]['bg_flag'].extend(df.loc[mask, 'Event'].values)
        
        fillucn_sum[key] += per_seg_fill[str(seg)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=160) # PE spectra by hold time
for hold_t in hold_times:
    all_pe = []
    for seg in SEGMENTS:
        key = (hold_t, int(seg))
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
