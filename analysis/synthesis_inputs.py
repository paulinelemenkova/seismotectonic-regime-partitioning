#!/usr/bin/env python3
import os
import sys

import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

cat_f, cmt_f, kos_f, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
c = pd.read_csv(cat_f)
d = 0.5
c["cx"] = np.floor(c.longitude / d) * d + d / 2
c["cy"] = np.floor(c.latitude / d) * d + d / 2
code = {"NEC": 1, "PBZ": 2, "WEST": 3}
sh = c[c.domain.isin(list(code))]
g = sh.groupby(["cx", "cy"]).domain.agg(lambda x: x.value_counts().idxmax()).reset_index()
g = g.merge(sh.groupby(["cx", "cy"]).size().reset_index(name="n"))
g = g[g.n >= 2]
with open(os.path.join(out, "dom.xyz"), "w") as f:
    for r in g.itertuples():
        f.write(f"{r.cx} {r.cy} {code[r.domain]}\n")
h = c[c.domain == "HKP"].groupby(["cx", "cy"]).size().reset_index(name="n")
with open(os.path.join(out, "hkp.xyz"), "w") as f:
    for r in h[h.n >= 5].itertuples():
        f.write(f"{r.cx} {r.cy} 1\n")
cols = "lon lat depth mrr mtt mff mrt mrf mtf exp date mw cls strike dip rake m0".split()
m = pd.read_csv(cmt_f, sep=r"\s+", comment="#", names=cols)

def feats(lon, lat, dep, w=3.0):
    return np.c_[(lon - 66.0) * 111.32 * np.cos(np.radians(33.0)), (lat - 33.0) * 110.57, dep * w]

m["dom"] = KNeighborsClassifier(n_neighbors=15).fit(
    feats(c.longitude, c.latitude, c.depth), c.domain).predict(feats(m.lon, m.lat, m.depth))
pos = {"HKP": (73.3, 35.5), "NEC": (73.7, 39.9), "PBZ": (66.9, 29.2), "WEST": (59.6, 36.8)}
with open(os.path.join(out, "meca.txt"), "w") as f:
    for D, (x, y) in pos.items():
        s = m[(m.dom == D) & ((m.depth <= 70) | (D == "HKP"))]
        T = np.array([(s[k] * 10.0 ** s.exp).sum() for k in ["mrr", "mtt", "mff", "mrt", "mrf", "mtf"]])
        T = T / np.abs(T).max()
        f.write(f"{x} {y} 10 " + " ".join(f"{v:.4f}" for v in T) + " 26\n")
k = np.loadtxt(kos_f)
with open(os.path.join(out, "cross.txt"), "w") as f:
    for r in k:
        mx = max(abs(r[6]), abs(r[7]))
        f.write(f"{r[0]} {r[1]} {r[6] / mx:.3f} {r[7] / mx:.3f} {r[8]:.1f}\n")
