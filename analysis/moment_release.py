#!/usr/bin/env python3
import json
import os
import sys

import numpy as np
import pandas as pd

cat = pd.read_csv(sys.argv[1])
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
cat["time"] = pd.to_datetime(cat.time, utc=True, format="mixed")
cat["m0"] = 10 ** (1.5 * cat.mw + 9.1)
cat["yr"] = cat.time.dt.year
D = 0.5
cat["cx"] = np.floor(cat.longitude / D) * D + D / 2
cat["cy"] = np.floor(cat.latitude / D) * D + D / 2
cells = cat.groupby(["cx", "cy"]).m0.sum().reset_index()
cells["lg"] = np.log10(cells.m0)
cells[["cx", "cy", "lg"]].to_csv(os.path.join(out, "moment_cells.xyz"), sep=" ", header=False, index=False,
                                 float_format="%.3f")
doms = ["HKP", "NEC", "PBZ", "WEST"]
years = np.arange(1900, 2026)
ann = pd.DataFrame({"year": years})
for d in doms:
    s = cat[cat.domain == d].groupby("yr").m0.sum().reindex(years, fill_value=0).cumsum()
    ann[d] = s.values
ann.to_csv(os.path.join(out, "moment_annual.txt"), sep=" ", header=False, index=False, float_format="%.4e")
tot = cat.m0.sum()
summ = {"total_Nm": float(tot), "mw_equiv_total": float((np.log10(tot) - 9.1) / 1.5), "domains": {}}
for d in doms:
    s = cat[cat.domain == d]
    top = s.nlargest(3, "m0")
    summ["domains"][d] = {"share": round(float(s.m0.sum() / tot), 3), "n": int(len(s)),
                          "top3_share": round(float(top.m0.sum() / s.m0.sum()), 3),
                          "top": [(str(t.time.date()), round(float(t.mw), 2)) for t in top.itertuples()],
                          "post2000_share": round(float(s[s.yr >= 2000].m0.sum() / s.m0.sum()), 3)}
cells_sorted = cells.sort_values("m0", ascending=False)
summ["n_cells"] = int(len(cells))
summ["top5_cells"] = cells_sorted.head(5)[["cx", "cy", "lg"]].values.tolist()
summ["top10pct_cells_share"] = round(float(cells_sorted.m0.head(max(1, len(cells) // 10)).sum() / tot), 3)
for name, date in (("E1", "2022-06-21"), ("E2", "2023-10-07"), ("E3", "2025-08-31"), ("E4", "2025-11-0")):
    e = cat[cat.time.astype(str).str.startswith(date)].nlargest(1, "mw")
    if len(e):
        summ[name] = {"mw": float(e.mw.iloc[0]), "domain": e.domain.iloc[0], "lon": float(e.longitude.iloc[0]),
                      "lat": float(e.latitude.iloc[0])}
json.dump(summ, open(os.path.join(out, "moment_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
