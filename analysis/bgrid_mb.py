#!/usr/bin/env python3
import json
import os
import sys

import numpy as np
import pandas as pd

N, RMAX, NMIN, DM = 150, 150.0, 50, 0.1
cat = pd.read_csv(sys.argv[1])
cat["time"] = pd.to_datetime(cat.time, utc=True, format="mixed")
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
d = cat[cat.nnd_bg & (cat.time.dt.year >= 1973) & (cat.magType == "mb")].copy()
d["m"] = np.round(d.mag.values, 1)

def dist_km(la0, lo0, la, lo):
    p0, p = np.radians(la0), np.radians(la)
    dl = np.radians(lo - lo0)
    a = np.sin((p - p0) / 2) ** 2 + np.cos(p0) * np.cos(p) * np.sin(dl / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))

def mc_maxc(m):
    bins = np.round(np.arange(m.min(), m.max() + DM, DM), 1)
    h = np.array([(m == b).sum() for b in bins])
    return float(bins[np.argmax(h)]) + 0.2

def bfit(m, mc):
    x = m[m >= mc - 1e-9]
    n = len(x)
    if n < NMIN:
        return None
    mean = x.mean()
    b = np.log10(np.e) / (mean - (mc - DM / 2))
    sb = 2.3 * b ** 2 * np.sqrt(((x - mean) ** 2).sum() / (n * (n - 1)))
    return b, sb, n

rows = []
for layer, sub in (("crust", d[d.depth < 70]), ("deep", d[d.depth >= 70])):
    la, lo, m = sub.latitude.values, sub.longitude.values, sub.m.values
    for x in np.arange(58.25, 76.0, 0.25):
        for y in np.arange(27.25, 41.0, 0.25):
            r = dist_km(y, x, la, lo)
            idx = np.argsort(r)[:N]
            idx = idx[r[idx] <= RMAX]
            if len(idx) < NMIN:
                continue
            mm = m[idx]
            mc = mc_maxc(mm)
            f = bfit(mm, mc)
            if f:
                rows.append((round(x, 2), round(y, 2), layer, round(f[0], 3), round(f[1], 3), mc, f[2],
                             round(float(r[idx].max()), 1)))
g = pd.DataFrame(rows, columns=["lon", "lat", "layer", "b", "sb", "mc", "n", "rmax_km"])
g.to_csv(os.path.join(out, "bgrid_mb.csv"), index=False)
summ = {}
for L in ("crust", "deep"):
    s = g[g.layer == L]
    summ[L] = {"nodes": int(len(s)), "events": int((d.depth < 70).sum() if L == "crust" else (d.depth >= 70).sum()),
               "b_min": float(s.b.min()), "b_med": float(s.b.median()), "b_max": float(s.b.max()),
               "sb_med": float(s.sb.median()), "mc_med": float(s.mc.median())}
    for name, q in (("low", s.nsmallest(5, "b")), ("high", s.nlargest(5, "b"))):
        summ[L][name] = q[["lon", "lat", "b", "sb", "n"]].values.tolist()
json.dump(summ, open(os.path.join(out, "bgrid_mb_summary.json"), "w"), indent=1)
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("low", "high")} for k, v in summ.items()}, indent=1))
for L in ("crust", "deep"):
    print(L, "low", summ[L]["low"][:3], "high", summ[L]["high"][:3])
