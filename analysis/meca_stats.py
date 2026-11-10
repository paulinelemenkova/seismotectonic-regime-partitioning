#!/usr/bin/env python3
import json
import sys

import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

def feats(lon, lat, dep, w=3.0):
    return np.c_[(lon - 66.0) * 111.32 * np.cos(np.radians(33.0)), (lat - 33.0) * 110.57, dep * w]

cat = pd.read_csv(sys.argv[1])
cols = "lon lat depth mrr mtt mff mrt mrf mtf exp date mw cls strike dip rake m0".split()
g = pd.read_csv(sys.argv[2], sep=r"\s+", comment="#", names=cols)

knn = KNeighborsClassifier(n_neighbors=15).fit(
    feats(cat.longitude, cat.latitude, cat.depth), cat.domain)
g["domain"] = knn.predict(feats(g.lon, g.lat, g.depth))

def p_axis(r):
    m = np.array([[r.mrr, r.mrt, r.mrf], [r.mrt, r.mtt, r.mtf], [r.mrf, r.mtf, r.mff]])
    w, v = np.linalg.eigh(m)
    up, south, east = v[:, 0]
    vec = np.array([-south, east, -up])
    if vec[2] < 0:
        vec = -vec
    return vec / np.linalg.norm(vec)

def az_pl(vec):
    n, e, d = vec
    if d < 0:
        n, e, d = -n, -e, -d
    return (np.degrees(np.arctan2(e, n)) % 360, np.degrees(np.arcsin(min(d, 1.0))))

out = {}
for D in ["HKP", "NEC", "PBZ", "WEST", "ALL"]:
    s = g if D == "ALL" else g[g.domain == D]
    if D != "HKP" and D != "ALL":
        s = s[s.depth <= 70]
    vecs = np.array([p_axis(r) for r in s.itertuples()])
    T = vecs.T @ vecs
    ew, ev = np.linalg.eigh(T)
    az, pl = az_pl(ev[:, -1])
    cnt = s.cls.value_counts()
    n = len(s)
    out[D] = {"n": int(n),
              **{k: [int(cnt.get(k, 0)), round(100 * cnt.get(k, 0) / max(n, 1), 1)]
                 for k in ("TF", "SS", "NF", "OB")},
              "p_az": round(float(az), 0), "p_pl": round(float(pl), 0),
              "p_conc": round(float(ew[-1] / ew.sum()), 2),
              "m0_sum_Nm": float(s.m0.sum()), "mw_max": float(s.mw.max())}
    print(D, out[D])
json.dump(out, open(sys.argv[3], "w"), indent=1)
