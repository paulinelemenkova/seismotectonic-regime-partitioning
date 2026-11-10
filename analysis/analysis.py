#!/usr/bin/env python3
import json, os, sys
import numpy as np, pandas as pd
from sklearn.cluster import HDBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KNeighborsClassifier

CC, GEM, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(42)
R = {}

def dist_km(a, b, c, d):
    p = np.radians
    x = np.sin(p(c - a) / 2) ** 2 + np.cos(p(a)) * np.cos(p(c)) * np.sin(p(d - b) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(x))

cc, gem = pd.read_csv(CC), pd.read_csv(GEM)
for d in (cc, gem):
    d["time"] = pd.to_datetime(d["time"], format="ISO8601", utc=True)
    d["magType"] = d["magType"].fillna("unk").astype(str).str.lower()
tcc = cc["time"].values.astype("datetime64[ms]").astype(np.int64) / 1000.0
pairs, used = [], set()
for _, g in gem.iterrows():
    j = np.where(np.abs(tcc - g["time"].value / 1e9) <= 16)[0]
    if len(j):
        dd = dist_km(g.latitude, g.longitude, cc.latitude.values[j], cc.longitude.values[j])
        if dd.min() <= 50:
            k = j[np.argmin(dd)]
            used.add(k)
            pairs.append((cc.loc[k, "mag"], cc.loc[k, "magType"], g["mag"]))
pr = pd.DataFrame(pairs, columns=["m", "t", "mw"]).dropna()
gem["source"], rest = "ISC-GEM", cc.drop(index=list(used)).copy()
rest["source"] = "ComCat"
cat = pd.concat([gem, rest], ignore_index=True).sort_values("time").reset_index(drop=True)
ts = cat["time"].values.astype("datetime64[ms]").astype(np.int64) / 1000.0
drop = set()
for i in range(len(cat) - 1):
    j = i + 1
    while j < len(cat) and ts[j] - ts[i] <= 16:
        if dist_km(cat.latitude[i], cat.longitude[i], cat.latitude[j], cat.longitude[j]) <= 50:
            keep_i = (cat.source[i] == "ISC-GEM") or (cat.source[j] != "ISC-GEM" and cat.mag[i] >= cat.mag[j])
            drop.add(j if keep_i else i)
        j += 1
R["n_residual_duplicates"] = len(drop)
cat = cat.drop(index=list(drop)).reset_index(drop=True)
R["n_comcat"], R["n_iscgem"], R["n_matched"] = len(cc), len(gem), len(pr)

mb = pr[pr.t == "mb"]
x, y = mb.m.values, mb.mw.values
eta = 0.2
sxx, syy, sxy = np.var(x), np.var(y), np.cov(x, y, bias=True)[0, 1]
beta = (syy - eta * sxx + np.sqrt((syy - eta * sxx) ** 2 + 4 * eta * sxy ** 2)) / (2 * sxy)
alpha = y.mean() - beta * x.mean()
res = y - (0.85 * x + 1.03)
R["mb_pairs"] = {"n": int(len(mb)), "mb_min": float(x.min()), "mb_max": float(x.max()),
                 "gor_slope": round(float(beta), 3), "gor_intercept": round(float(alpha), 3),
                 "scordilis_resid_mean": round(float(res.mean()), 3),
                 "scordilis_resid_sd": round(float(res.std()), 3)}

def to_mw(m, t, rel="scordilis"):
    if t.startswith("mw"):
        return m
    if t in ("mb", "mblg", "mb_lg"):
        return 0.85 * m + 1.03 if rel == "scordilis" else beta * m + alpha
    if t.startswith("ms"):
        return 0.67 * m + 2.07 if m <= 6.1 else 0.99 * m + 0.08
    return m

cat["mw"] = np.round([to_mw(m, t) for m, t in zip(cat.mag, cat.magType)], 2)
cat["mw_gor"] = np.round([to_mw(m, t, "gor") for m, t in zip(cat.mag, cat.magType)], 2)
cat["mw_conv"] = ~cat.magType.str.startswith("mw")
cat = cat[cat.mw >= 4.0].reset_index(drop=True)
cat["depth"] = cat["depth"].fillna(33.0)
cat["fixed_depth"] = cat.depth.round(1).isin([10.0, 33.0, 35.0])
R["n_final"] = len(cat)
R["source_counts"] = {k: int(v) for k, v in cat.source.value_counts().items()}
R["n_converted"] = int(cat.mw_conv.sum())
R["n_fixed_depth"] = int(cat.fixed_depth.sum())
R["n_deep70"] = int((cat.depth >= 70).sum())

t = cat.time.values.astype("datetime64[s]").astype(np.int64) / 86400.0 / 365.25 + 1970.0
lat, lon, dep, mw = (cat[c].values for c in ("latitude", "longitude", "depth", "mw"))
dep_ = dep.copy()

gk = np.zeros(len(cat), bool)
for i in np.argsort(-mw):
    if gk[i]:
        continue
    L = 10 ** (0.1238 * mw[i] + 0.983)
    T = (10 ** (0.032 * mw[i] + 2.7389) if mw[i] >= 6.5 else 10 ** (0.5409 * mw[i] - 0.547)) / 365.25
    c = np.where((np.abs(t - t[i]) <= T) & (mw <= mw[i]))[0]
    c = c[c != i]
    if len(c):
        r = np.sqrt(dist_km(lat[i], lon[i], lat[c], lon[c]) ** 2 + (dep_[i] - dep_[c]) ** 2)
        gk[c[r <= L]] = True
cat["gk_main"] = ~gk

DF, BB = 1.6, 1.0
le = np.full(len(cat), np.nan)
for j in range(1, len(cat)):
    i = np.arange(j)
    dt = t[j] - t[i]
    ok = dt > 0
    if ok.any():
        i, dt = i[ok], dt[ok]
        r = np.maximum(np.sqrt(dist_km(lat[j], lon[j], lat[i], lon[i]) ** 2 + (dep_[j] - dep_[i]) ** 2), 1.0)
        le[j] = np.log10(np.min(dt * r ** DF * 10 ** (-BB * mw[i])))
v = le[np.isfinite(le)].reshape(-1, 1)
gm = GaussianMixture(2, random_state=0).fit(v)
mu = np.sort(gm.means_.ravel())
grid = np.linspace(mu[0], mu[1], 4000)
lo = int(np.argmin(gm.means_.ravel()))
thr = float(grid[np.argmin(np.abs(gm.predict_proba(grid.reshape(-1, 1))[:, lo] - 0.5))])
cat["log_eta"] = le
cat["nnd_bg"] = ~(cat.log_eta < thr)
R["nnd"] = {"threshold_log_eta": round(thr, 3), "modes": [round(float(m), 3) for m in mu],
            "n_background": int(cat.nnd_bg.sum())}
R["gk"] = {"n_mainshocks": int(cat.gk_main.sum())}

def feats(d, w=3.0):
    return np.c_[(d.longitude - 66.0) * 111.32 * np.cos(np.radians(33.0)),
                 (d.latitude - 33.0) * 110.57, d.depth * w]

bg = cat[cat.nnd_bg]
lab = HDBSCAN(min_cluster_size=80, min_samples=20, copy=True).fit_predict(feats(bg))
cl = []
for L in sorted(set(lab) - {-1}):
    s = bg[lab == L]
    cl.append({"label": int(L), "n": int(len(s)),
               "lon": round(float(s.longitude.median()), 2), "lat": round(float(s.latitude.median()), 2),
               "depth_med": round(float(s.depth.median()), 1),
               "fixed_share": round(float(s.fixed_depth.mean()), 2)})
R["hdbscan"] = {"n_clusters": len(cl), "n_noise": int((lab == -1).sum()), "clusters": cl}

def group(c):
    if c["depth_med"] >= 70:
        return "HKP"
    if c["lon"] < 64.5:
        return "WEST"
    if c["lat"] < 35.2 or c["lon"] > 72.5 and c["lat"] < 36.0:
        return "PBZ"
    return "NEC"

gmap = {c["label"]: group(c) for c in cl}
for c in cl:
    c["domain"] = gmap[c["label"]]
lab_dom = np.array([gmap.get(L, "NOISE") for L in lab])
knn = KNeighborsClassifier(n_neighbors=15).fit(feats(bg)[lab_dom != "NOISE"], lab_dom[lab_dom != "NOISE"])
cat["domain"] = knn.predict(feats(cat))
cat["hdb_label"] = -2
cat.loc[bg.index, "hdb_label"] = lab
DOMS = ["HKP", "NEC", "PBZ", "WEST"]
NAMES = {"HKP": "Intermediate-depth domain",
         "NEC": "North-eastern crustal domain", "PBZ": "Plate-boundary domain",
         "WEST": "Western interior domain"}

def mc_maxc(m, dm=0.1):
    m = np.round(m, 1)
    b = np.arange(m.min() - dm / 2, m.max() + dm, dm)
    h, e = np.histogram(m, b)
    return round(float(e[np.argmax(h)] + dm / 2 + 0.2), 2)

def bfit(m, mc, dm=0.1, nboot=1000):
    m = np.round(m, 1)
    x = m[m >= mc - 1e-9]
    n = len(x)
    if n < 30:
        return None
    b = np.log10(np.e) / (x.mean() - (mc - dm / 2))
    sb = 2.3 * b ** 2 * np.sqrt(((x - x.mean()) ** 2).sum() / (n * (n - 1)))
    bs = np.array([np.log10(np.e) / (RNG.choice(x, n).mean() - (mc - dm / 2)) for _ in range(nboot)])
    if nboot == 0:
        bs = np.array([b, b])
    return {"n": int(n), "mc": mc, "b": round(float(b), 3), "sb": round(float(sb), 3),
            "lo": round(float(np.percentile(bs, 2.5)), 3), "hi": round(float(np.percentile(bs, 97.5)), 3),
            "a": round(float(np.log10(n) + b * mc), 3), "boot": bs}

T0, T1, MC = 1964, 2026.0, 5.5
TC = T1 - T0
st = {}
for D in DOMS + ["ALL"]:
    d = cat if D == "ALL" else cat[cat.domain == D]
    dn = d[d.nnd_bg & (d.time.dt.year >= T0)]
    sub = dn[np.round(dn.mw, 1) >= MC]
    f = bfit(dn.mw.values, MC)
    out = {k: v for k, v in f.items() if k != "boot"}
    out["native_share"] = round(float((~sub.mw_conv).mean()), 3)
    out["mc_maxc_native"] = mc_maxc(dn[~dn.mw_conv].mw.values)
    rate = f["n"] / TC
    out["rate_mc"] = round(rate, 3)
    bl = f["boot"]
    for M in (6.0, 7.0):
        lam = rate * 10 ** (-f["b"] * (M - MC))
        lam_b = rate * 10 ** (-bl * (M - MC))
        k = str(int(M))
        out[f"ret_m{k}_yr"] = round(float(1 / lam), 1)
        out[f"ret_m{k}_ci"] = [round(float(np.percentile(1 / lam_b, q)), 1) for q in (2.5, 97.5)]
        out[f"p_m{k}_50yr"] = round(float(1 - np.exp(-lam * 50)), 3)
        out[f"p_m{k}_50yr_ci"] = [round(float(np.percentile(1 - np.exp(-lam_b * 50), q)), 3) for q in (2.5, 97.5)]
    out["n_obs_m6"] = int((d.mw >= 6).sum())
    out["n_obs_m7"] = int((d.mw >= 7).sum())
    out["n_obs_m6_1964"] = int(((d.mw >= 6) & (d.time.dt.year >= 1964)).sum())
    out["n_all"] = int(len(d))
    out["n_shallow"] = int((d.depth < 70).sum())
    dd = d[~d.fixed_depth].depth
    out["depth_q"] = [round(float(dd.quantile(q)), 1) for q in (0.1, 0.5, 0.9)]
    out["m0_sum"] = float((10 ** (1.5 * d.mw + 9.1)).sum())
    out["lon_med"], out["lat_med"] = round(float(d.longitude.median()), 2), round(float(d.latitude.median()), 2)
    sens = {}
    for lab_, sb_ in [("raw", d[d.time.dt.year >= T0]), ("gk", d[d.gk_main & (d.time.dt.year >= T0)]),
                      ("nnd_nofixed", dn[~dn.fixed_depth])]:
        g = bfit(sb_.mw.values, MC, nboot=300)
        sens[lab_] = None if g is None else [g["mc"], g["b"], g["sb"], g["n"]]
    for mc_ in (5.3, 5.8):
        g = bfit(dn.mw.values, mc_, nboot=300)
        sens[f"mc{mc_}"] = None if g is None else [g["mc"], g["b"], g["sb"], g["n"]]
    d73 = d[d.nnd_bg & (d.time.dt.year >= 1973)]
    for lab_, col in (("conv_scordilis", "mw"), ("conv_gor", "mw_gor")):
        g = bfit(d73[col].values, mc_maxc(d73[col].values), nboot=300)
        sens[lab_] = None if g is None else [g["mc"], g["b"], g["sb"], g["n"]]
    out["sensitivity"] = sens
    st[D] = out
    pd.DataFrame({"m": np.round(dn.mw.values, 1)}).to_csv(os.path.join(OUT, f"fmd_{D}.csv"), index=False)
tot = sum(st[D]["m0_sum"] for D in DOMS)
for D in DOMS:
    st[D]["m0_share"] = round(st[D]["m0_sum"] / tot, 3)
R["domains"] = st
R["names"] = NAMES

boots = {}
for D in DOMS:
    d = cat[(cat.domain == D) & cat.nnd_bg & (cat.time.dt.year >= T0)]
    boots[D] = bfit(d.mw.values, MC)["boot"]
R["b_pairwise_p"] = {f"{a}>{b}": round(float(np.mean(boots[a] > boots[b])), 3)
                     for i, a in enumerate(DOMS) for b in DOMS[i + 1:]}

h = cat[(cat.domain == "HKP") & ~cat.fixed_depth & cat.nnd_bg].depth.values.reshape(-1, 1)
g2 = GaussianMixture(2, random_state=0).fit(h)
R["hkp_gmm"] = sorted([[round(float(m), 1), round(float(np.sqrt(v)), 1), round(float(w), 2)]
                       for m, v, w in zip(g2.means_.ravel(), g2.covariances_.ravel(), g2.weights_)])

mt = cat[cat.time.dt.year >= 1964].sort_values("time")
rows = []
for s0 in range(0, len(mt) - 400, 100):
    w = mt.iloc[s0:s0 + 400]
    bsmc = [mc_maxc(RNG.choice(w.mw.values, len(w))) for _ in range(100)]
    rows.append((w.time.iloc[200].year + w.time.iloc[200].dayofyear / 365.25, mc_maxc(w.mw.values),
                 np.percentile(bsmc, 16), np.percentile(bsmc, 84)))
pd.DataFrame(rows, columns=["year", "mc", "lo", "hi"]).to_csv(os.path.join(OUT, "mc_time.csv"), index=False)

gb = []
dn = cat[cat.nnd_bg & (cat.time.dt.year >= 1973)]
for lo_ in np.arange(58.25, 76.0, 0.25):
    for la_ in np.arange(27.25, 41.0, 0.25):
        for layer, sub in (("crust", dn[dn.depth < 70]), ("deep", dn[dn.depth >= 70])):
            r = dist_km(la_, lo_, sub.latitude.values, sub.longitude.values)
            idx = np.argsort(r)[:150]
            idx = idx[r[idx] <= 150]
            if len(idx) < 80:
                continue
            m = sub.mw.values[idx]
            mc = mc_maxc(m)
            f = bfit(m, mc, nboot=0) if (np.round(m, 1) >= mc).sum() >= 50 else None
            if f:
                gb.append((lo_, la_, layer, f["b"], f["sb"], mc, f["n"]))
pd.DataFrame(gb, columns=["lon", "lat", "layer", "b", "sb", "mc", "n"]).to_csv(os.path.join(OUT, "bgrid.csv"), index=False)

cat.to_csv(os.path.join(OUT, "catalogue_final.csv"), index=False)
json.dump(R, open(os.path.join(OUT, "stats.json"), "w"), indent=1, default=str)
print(json.dumps({k: R[k] for k in ("n_final", "source_counts", "n_converted", "mb_pairs", "nnd", "gk")}, default=str))
for D in DOMS + ["ALL"]:
    s = st[D]
    print(D, s["n_all"], "b", s["b"], "+-", s["sb"], [s["lo"], s["hi"]], "n", s["n"], "nat", s["native_share"],
          "ret6", s["ret_m6_yr"], s["ret_m6_ci"], "ret7", s["ret_m7_yr"], s["ret_m7_ci"], "P7_50", s["p_m7_50yr"], s["p_m7_50yr_ci"],
          "obs6", s["n_obs_m6"], s["n_obs_m6_1964"], "obs7", s["n_obs_m7"], "share", s.get("m0_share"), "dq", s["depth_q"], "nsh", s["n_shallow"])
    print("   sens", s["sensitivity"])
print(R["b_pairwise_p"], R["hkp_gmm"])
print([(c["label"], c["n"], c["domain"], c["lon"], c["lat"], c["depth_med"], c["fixed_share"]) for c in cl])
