import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ANALYSIS_DIR = './output/results/'
OUT_CSV = os.path.join(ANALYSIS_DIR, "pe_distribution_by_segment.csv")

csv_files = [
    fp for fp in glob.glob(os.path.join(ANALYSIS_DIR, "*.csv"))
    if "PulseWindowStats" not in os.path.basename(fp)
]

usecols = ["Segment", " PE"]
df_list = []
for fp in sorted(csv_files):
    try:
        tmp = pd.read_csv(fp, usecols=usecols)
        df_list.append(tmp)
    except Exception as e:
        print(f"Skip {fp}: {e}")

df = pd.concat(df_list, ignore_index=True)

df["PE"] = pd.to_numeric(df[" PE"], errors="coerce").fillna(0).astype(int)

counts = df.groupby(["PE", "Segment"]).size().reset_index(name="count")

wide = counts.pivot(index="PE", columns="Segment", values="count").fillna(0).astype(int)

rename_map = {seg: f"Segment_{seg}" for seg in wide.columns}
wide = wide.rename(columns=rename_map)

wide.to_csv(OUT_CSV)
print(f"Saved PE distribution -> {OUT_CSV}")