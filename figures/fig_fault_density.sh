#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
GDRIVE_ID=${GDRIVE_ID:-1TLTNVIpkC4p7Mff4VkBo8mUVRfA3Yh6v}
W=60; E=75.5; S=29; N=38.8; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work

curl -s -f -m 600 -L -o afg.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o afg.zip -x '__MACOSX/*'
P=$(find . -name 'fltafg.shp' | head -1)
ogr2ogr -f OGR_GMT flt.gmt "$P" -t_srs EPSG:4326
curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
ogr2ogr -f OGR_GMT qf.gmt "$(find . -name 'faults.shp' | head -1)" -t_srs EPSG:4326
gmt coast -EAF -M -Di > af.txt
curl -s -m 300 "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=$S&maxlatitude=$N&minlongitude=$W&maxlongitude=$E&minmagnitude=4&maxdepth=40&starttime=1900-01-01&endtime=2025-12-31T23:59:59" \
  | awk -F, 'NR>1{print $3, $2, $4, $5}' > eq_all.txt
gmt select eq_all.txt -Faf.txt > eq_afg.txt
NEQ=$(wc -l < eq_afg.txt)
python3 - $NEQ <<'PY'
import sys, numpy as np
rng = np.random.default_rng(42)
n = 20 * int(sys.argv[1])
x = rng.uniform(60.4, 75.0, n); y = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(29.3)), np.sin(np.radians(38.5)), n)))
np.savetxt("rnd_box.txt", np.c_[x, y], fmt="%.4f")
PY
gmt select rnd_box.txt -Faf.txt | head -n $((NEQ * 5)) > rnd.txt
python3 - <<'PY'
import urllib.request
B = "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
files = ["jan76_dec20.ndk"] + [f"NEW_MONTHLY/{y}/{m}{y % 100:02d}.ndk" for y in range(2021, 2026)
                               for m in "jan feb mar apr may jun jul aug sep oct nov dec".split()]
out = open("cmt_sdr.txt", "w")
for f in files:
    try:
        t = urllib.request.urlopen(urllib.request.Request(B + f, headers={"User-Agent": "afg-paper/1.0"}),
                                   timeout=180).read().decode("latin-1").splitlines()
    except Exception:
        continue
    for i in range(0, len(t) - 4, 5):
        c = t[i + 2].split()
        if not t[i + 2].startswith("CENTROID") or int(t[i][5:9]) > 2025:
            continue
        lat, lon, dep = float(c[3]), float(c[5]), float(c[7])
        if dep > 40 or not (60 <= lon <= 75.5 and 29 <= lat <= 38.8):
            continue
        v = t[i + 4].split()
        out.write(f"{lon} {lat} {v[11]} {v[14]}\n")
out.close()
PY
gmt select cmt_sdr.txt -Faf.txt > cmt_afg.txt

python3 - <<'PY'
import json, numpy as np
from scipy import stats
from scipy.spatial import cKDTree
R = 6371.0
def read_lines(fn):
    segs, cur = [], []
    for ln in open(fn):
        if ln.startswith("#"):
            continue
        if ln.startswith(">"):
            if len(cur) > 1: segs.append(np.array(cur))
            cur = []; continue
        p = ln.split()
        if len(p) >= 2: cur.append((float(p[0]), float(p[1])))
    if len(cur) > 1: segs.append(np.array(cur))
    return segs
def pieces(segs):
    out = []
    for s in segs:
        lo1, la1, lo2, la2 = s[:-1, 0], s[:-1, 1], s[1:, 0], s[1:, 1]
        p1, p2 = np.radians(la1), np.radians(la2); dl = np.radians(lo2 - lo1)
        a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
        L = 2 * R * np.arcsin(np.sqrt(a))
        az = np.degrees(np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))) % 180
        out.append(np.c_[(lo1 + lo2) / 2, (la1 + la2) / 2, L, az])
    return np.vstack(out)
F = pieces(read_lines("flt.gmt")); Q = pieces(read_lines("qf.gmt"))
d = 0.25; xs = np.arange(60, 75.5, d); ys = np.arange(29, 38.75, d)
G = np.zeros((len(ys), len(xs)))
ix = np.floor((F[:, 0] - 60) / d).astype(int); iy = np.floor((F[:, 1] - 29) / d).astype(int)
ok = (ix >= 0) & (ix < len(xs)) & (iy >= 0) & (iy < len(ys))
np.add.at(G, (iy[ok], ix[ok]), F[ok, 2])
area = (R ** 2) * np.radians(d) * (np.sin(np.radians(ys + d)) - np.sin(np.radians(ys)))
G = G / area[:, None]
with open("dens.xyz", "w") as f:
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            if G[j, i] > 0: f.write(f"{x + d / 2:.3f} {y + d / 2:.3f} {G[j, i]:.5f}\n")
def rose(az, w, fn):
    h, _ = np.histogram(az % 180, bins=np.arange(0, 181, 10), weights=w)
    h = h / h.max()
    with open(fn, "w") as f:
        for k, v in enumerate(h): f.write(f"{v:.4f} {5 + 10 * k}\n")
    return int(np.argmax(h)) * 10 + 5
pk_f = rose(F[:, 3], F[:, 2], "rose_flt.txt"); pk_q = rose(Q[:, 3], Q[:, 2], "rose_qf.txt")
c = np.loadtxt("cmt_afg.txt", ndmin=2)
st = np.r_[c[:, 2], c[:, 3]] % 180
pk_c = rose(st, np.ones_like(st), "rose_cmt.txt")
def proj(lon, lat):
    return np.c_[np.radians(lon - 67.5) * R * np.cos(np.radians(lat)), np.radians(lat) * R]
vx = []
for sgm in read_lines("qf.gmt"):
    for (x1, y1), (x2, y2) in zip(sgm[:-1], sgm[1:]):
        L = np.hypot((x2 - x1) * 111.2 * np.cos(np.radians(y1)), (y2 - y1) * 111.2)
        k = max(2, int(np.ceil(L / 0.5)) + 1)
        t = np.linspace(0, 1, k); vx.append(np.c_[x1 + t * (x2 - x1), y1 + t * (y2 - y1)])
vx = np.vstack(vx); tree = cKDTree(proj(vx[:, 0], vx[:, 1]))
pe = np.loadtxt("eq_afg.txt", ndmin=2); pr = np.loadtxt("rnd.txt", ndmin=2)
de = tree.query(proj(pe[:, 0], pe[:, 1]))[0]; dr = tree.query(proj(pr[:, 0], pr[:, 1]))[0]
ks = stats.ks_2samp(de, dr)
for name, arr in (("eq", de), ("rnd", dr)):
    s = np.sort(arr); np.savetxt(f"cdf_{name}.txt", np.c_[s, np.arange(1, len(s) + 1) / len(s)], fmt="%.4f")
S = {"fault_length_km": float(F[:, 2].sum()), "qfault_length_km": float(Q[:, 2].sum()),
     "dens_max": float(G.max()), "dens_mean_nonzero": float(G[G > 0].mean()),
     "dens_p90": float(np.percentile(G[G > 0], 90)), "n_cells": int((G > 0).sum()),
     "peak_strike_flt": pk_f, "peak_strike_qf": pk_q, "peak_strike_cmt": pk_c, "n_cmt": int(len(c)),
     "n_eq": int(len(de)), "n_rnd": int(len(dr)),
     "med_eq": float(np.median(de)), "med_rnd": float(np.median(dr)),
     "frac_eq_10km": float((de <= 10).mean()), "frac_rnd_10km": float((dr <= 10).mean()),
     "frac_eq_25km": float((de <= 25).mean()), "frac_rnd_25km": float((dr <= 25).mean()),
     "ks_D": float(ks.statistic), "ks_p": float(ks.pvalue)}
def dens_at(p):
    i = np.clip(np.floor((p[:, 0] - 60) / d).astype(int), 0, len(xs) - 1); j = np.clip(np.floor((p[:, 1] - 29) / d).astype(int), 0, len(ys) - 1)
    return G[j, i]
S["dens_med_eq"] = float(np.median(dens_at(pe[:, :2])))
S["dens_med_rnd"] = float(np.median(dens_at(pr[:, :2])))
json.dump(S, open("fault_stats.json", "w"), indent=1); print(json.dumps(S, indent=1))
open("n.txt", "w").write(f"{S['n_eq']} {S['n_rnd']} {S['n_cmt']} {S['med_eq']:.1f} {S['med_rnd']:.1f} {S['ks_D']:.2f} {S['ks_p']:.1e}\n")
PY
read NEQ NRND NCMT MEQ MRND KSD KSP < n.txt
mkdir -p "$OUT"; cp fault_stats.json "$OUT/"

gmt xyz2grd dens.xyz -R60/75.5/29/38.75 -I0.25 -r -Gdens.nc
gmt makecpt -Cjet -T0/0.6/0.05 > dens.cpt
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt

gmt begin $OUT/fig_fault_density pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8.5p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray40 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y11.0c
gmt grdimage dens.nc -Cdens.cpt -Q -t25 -nn
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot qf.gmt -W0.5p,black
gmt plot af.txt -W1.5p,black
gmt plot eq_afg.txt -Sc0.09c -Gwhite -W0.35p,black
printf "60.15 38.65 (a) Fault density and shallow epicentres\n" | gmt text -F+f9p,Helvetica-Bold+jTL -Gwhite@15 -C0.06c/0.05c
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.1c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
gmt colorbar -DJBC+w10c/0.3c+o-3.0c/0.75c+h+ef -Cdens.cpt -Bxa0.1f0.05+l"Fault density (km km@+-2@+)"
cat > leg.txt <<EOF
N 1
S 0.25c c 0.12c white 0.35p,black 0.55c Epicentre, depth < 40 km ($NEQ)
S 0.25c - 0.45c - 0.6p,black 0.55c Quaternary fault (USGS)
EOF
gmt legend leg.txt -DJBR+jTR+w5.4c+o0/0.55c -F+gwhite+p0.3p,gray50 --FONT_ANNOT_PRIMARY=7.5p,Helvetica
gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 6.5p,Helvetica MAP_GRID_PEN_PRIMARY 0.25p,gray70
gmt rose rose_flt.txt -R0/1/0/360 -JX3.3c -X0.3c -Y-6.9c -A10 -T -L-,-,-,N -G0/70/255@20 -W0.4p,black \
    -Bxg0.25 -Byg30 -B+t"(b) Mapped faults"
gmt rose rose_qf.txt -R0/1/0/360 -JX3.3c -X3.9c -A10 -T -L-,-,-,N -G220/20/60@20 -W0.4p,black \
    -Bxg0.25 -Byg30 -B+t"(c) Quaternary faults"
gmt rose rose_cmt.txt -R0/1/0/360 -JX3.3c -X3.9c -A10 -T -L-,-,-,N -G240/160/0@15 -W0.4p,black \
    -Bxg0.25 -Byg30 -B+t"(d) CMT nodal planes ($NCMT)"
gmt set FONT_ANNOT_PRIMARY 7p,Helvetica MAP_GRID_PEN_PRIMARY 0.25p,gray80
gmt basemap -R0/100/0/1 -JX4.0c/3.6c -X4.9c -Y-0.1c -Bxa25f5g25+l"Distance to Quaternary fault (km)" \
    -Bya0.2f0.1g0.2+l"Cumulative fraction" -BWSen+t"(e) Epicentre-fault distance"
gmt plot cdf_rnd.txt -W1.4p,gray40,6_3:0p
gmt plot cdf_eq.txt -W1.6p,220/20/60
cat > leg2.txt <<EOF
S 0.2c - 0.4c - 1.6p,220/20/60 0.5c Epicentres, $MEQ km
S 0.2c - 0.4c - 1.4p,gray40,6_3:0p 0.5c Random, $MRND km
L 6.5p,Helvetica L (medians); KS D = $KSD
EOF
gmt legend leg2.txt -DjBR+w2.7c+o0.06c -F+gwhite@10+p0.3p,gray50 --FONT_ANNOT_PRIMARY=6.5p,Helvetica
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-13.0c -Y-2.9c <<EOF
0 1.55 Depicts: (a) length of USGS mapped faults per 0.25-degree cell with ComCat epicentres (M >= 4, depth < 40 km, 1900-2025) inside the study area;
0 1.15 (b-d) length-weighted strikes of mapped and Quaternary faults and Global CMT nodal-plane strikes (<= 40 km), 10-degree sectors, normalised;
0 0.75 (e) distance to the nearest Quaternary fault for $NEQ epicentres and $NRND uniformly random points. Data: USGS faults; USGS OFR 2007-1103; ComCat; Global CMT.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python (NumPy, SciPy). Source: authors.
EOF
gmt end
