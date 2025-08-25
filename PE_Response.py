import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

ANALYSIS_DIR = './output/tail/'
GRAPH_DIR = './output/graphs/'

df_list = []

for filepath in sorted(glob.glob(os.path.join(ANALYSIS_DIR, "*.csv"))):
    run_id = os.path.splitext(os.path.basename(filepath))[0]
    temp_df = pd.read_csv(filepath)
    temp_df["Run"] = run_id
    df_list.append(temp_df)

all_runs_df = pd.concat(df_list, ignore_index=True)

print(all_runs_df.head())
print("Total rows:", len(all_runs_df))

time_col = "Time(us)"
seg_cols = [c for c in all_runs_df.columns if c.startswith("Segment_")]
g = all_runs_df.groupby(time_col, as_index=False)[seg_cols].sum()
g = g.sort_values(time_col)
x = g[time_col].values
bin_w = 0.1

csv_outpath = os.path.join(ANALYSIS_DIR, "all_tail_response.csv")
g.to_csv(csv_outpath, index=False)
print("Saved CSV: ", csv_outpath)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), dpi=160, sharex=True)
colors = ["red", "blue", "green", "magenta"]

def draw_panel(ax, logy=False):
    for seg, color in zip(seg_cols, colors):
        y = g[seg].values
        ax.plot(x, y, color=color, linewidth=1.2, label=seg)
    if logy:
        ax.set_yscale("log")
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel("Time after pulse (µs)")
    ax.set_ylabel("Counts")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(title="Segment", loc="best", fontsize=12)

draw_panel(ax1, logy=False)
ax1.set_title("Summed Tail Response (linear)")

draw_panel(ax2, logy=True)
ax2.set_title("Summed Tail Response (log y)")

plt.tight_layout()
outpath = os.path.join(GRAPH_DIR, "PE_Response/summed_tail_response_line.png")
plt.savefig(outpath)
plt.close()
print("Saved:", outpath)