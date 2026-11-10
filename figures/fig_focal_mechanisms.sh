#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
W=59; E=76; S=28; N=40; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work

python3 - "$W" "$E" "$S" "$N" <<'PY'
import sys, math, urllib.request
W, E, S, N = map(float, sys.argv[1:5])
B = "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
files = ["jan76_dec20.ndk"]
mon = "jan feb mar apr may jun jul aug sep oct nov dec".split()
for y in range(2021, 2026):
    files += [f"NEW_MONTHLY/{y}/{m}{y % 100:02d}.ndk" for m in mon]
txt = []
for f in files:
    try:
        r = urllib.request.urlopen(urllib.request.Request(B + f, headers={"User-Agent": "afg-paper/1.0"}),
                                   timeout=180).read().decode("latin-1")
        txt.append(r)
    except Exception as e:
        print("skip", f, e, file=sys.stderr)
lines = "".join(txt).splitlines()
rows = []
for i in range(0, len(lines) - 4, 5):
    l1, l2, l3, l4, l5 = lines[i:i + 5]
    if not l3.startswith("CENTROID"):
        continue
    year = int(l1[5:9])
    if year > 2025:
        continue
    c = l3.split(); lat, lon, dep = float(c[3]), float(c[5]), float(c[7])
    if not (W <= lon <= E and S <= lat <= N):
        continue
    m = l4.split(); ex = int(m[0]); mt = [float(m[k]) for k in (1, 3, 5, 7, 9, 11)]
    v = l5.split()
    tpl, taz = float(v[2]), float(v[3]); bpl = float(v[5]); ppl, paz = float(v[8]), float(v[9])
    m0 = float(v[10]) * 10 ** ex
    mw = 2.0 / 3.0 * (math.log10(m0) - 16.1)
    s2 = [math.sin(math.radians(x)) ** 2 for x in (tpl, ppl, bpl)]
    tot = sum(s2); s2 = [100 * x / tot for x in s2]
    if tpl >= 50: k = "TF"
    elif ppl >= 60: k = "NF"
    elif bpl >= 60: k = "SS"
    else: k = "OB"
    rows.append(dict(y=year, date=l1[5:15], lon=lon, lat=lat, dep=dep, mt=mt, ex=ex, mw=mw, k=k,
                     s2=s2, paz=paz % 180, sdr=" ".join(v[11:14]), m0=m0))
rows.sort(key=lambda r: -r["mw"])
with open("cmt_all.txt", "w") as f:
    f.write("# lon lat depth mrr mtt mff mrt mrf mtf exp date mw class strike dip rake m0_Nm\n")
    for r in rows:
        f.write(f"{r['lon']} {r['lat']} {r['dep']} " + " ".join(f"{x}" for x in r["mt"]) +
                f" {r['ex']} {r['date']} {r['mw']:.2f} {r['k']} {r['sdr']} {r['m0'] * 1e-7:.4e}\n")
for k in ("TF", "NF", "SS", "OB"):
    for d, sel in (("sh", lambda r: r["dep"] <= 70), ("in", lambda r: r["dep"] > 70)):
        sub = [r for r in rows if r["k"] == k and sel(r)]
        with open(f"meca_{d}_{k}.txt", "w") as f:
            for r in sub:
                f.write(f"{r['lon']} {r['lat']} {r['dep']} " + " ".join(f"{x}" for x in r["mt"]) + f" {r['ex']}\n")
        with open(f"tern_{d}_{k}.txt", "w") as f:
            for r in sub:
                f.write(f"{r['s2'][1]:.2f} {r['s2'][2]:.2f} {r['s2'][0]:.2f}\n")
with open("meca_in_all.txt", "w") as f:
    for r in rows:
        if r["dep"] > 70:
            f.write(f"{r['lon']} {r['lat']} {r['dep']} " + " ".join(f"{x}" for x in r["mt"]) + f" {r['ex']}\n")
sh = [r for r in rows if r["dep"] <= 70]
with open("paz_sh.txt", "w") as f:
    for r in sh:
        f.write(f"1 {r['paz']:.1f}\n")
cnt = [0] * 18
for r in sh:
    cnt[int(r["paz"] // 10) % 18] += 1
rmax = 5 * math.ceil(max(cnt) / 5 + 0.2)
nin = len(rows) - len(sh)
st = {k: (sum(1 for r in sh if r["k"] == k), sum(1 for r in rows if r["dep"] > 70 and r["k"] == k))
      for k in ("TF", "NF", "SS", "OB")}
print("total", len(rows), "shallow", len(sh), "intermediate", nin)
for k, v in st.items():
    print(k, "shallow", v[0], f"{100 * v[0] / max(len(sh), 1):.1f}%", "interm", v[1], f"{100 * v[1] / max(nin, 1):.1f}%")
print("years", min(r["y"] for r in rows), max(r["y"] for r in rows),
      "Mw", f"{min(r['mw'] for r in rows):.1f}-{max(r['mw'] for r in rows):.1f}")
print("largest shallow:", [(r["date"], round(r["mw"], 1), r["k"], round(r["dep"])) for r in sh[:8]])
open("n.txt", "w").write(f"{len(rows)} {len(sh)} {nin} {rmax} " +
                         " ".join(str(st[k][0]) for k in ("TF", "NF", "SS", "OB")) + "\n")
PY
read NALL NSH NIN RMAX NTF NNF NSS NOB < n.txt
cp cmt_all.txt "$OUT/gcmt_ndk_subset.txt"

curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
SHP=$(find . -name 'faults.shp' | head -1)
ogr2ogr -f OGR_GMT qf.gmt "$SHP" -t_srs EPSG:4326
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt coast -EAF -M -Di > af.txt
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt

C_TF=220/20/60; C_NF=0/70/255; C_SS=255/215/0; C_OB=0/225/235

gmt begin $OUT/fig_focal_mechanisms pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 9p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray55 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y12.5c
gmt coast -Di -W0.3p,gray30 -N1/0.6p,gray25,- -I1/0.35p,steelblue -C+l+gsteelblue -A50 -Bxa2f1g2 -Bya2f1g2 -BWSen
gmt plot qf.gmt -W0.45p,gray20
gmt plot af.txt -W1.4p,black
gmt meca meca_in_all.txt -Sm0.22c -Ggray30 -Ewhite -W0.15p,gray10
for k in TF NF SS OB; do
  eval col=\$C_$k
  gmt meca meca_sh_$k.txt -Sm0.38c -G$col -Ewhite -W0.25p,black
done
printf "59.15 39.85 (a)\n" | gmt text -F+f11p,Helvetica-Bold+jTL -Gwhite@20 -C0.06c/0.05c
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.2c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<EOF
H 9p,Helvetica-Bold Shallow (<= 70 km) mechanisms, class after Frohlich (1992)
N 4
S 0.25c c 0.3c $C_TF 0.3p,black 0.5c Thrust ($NTF)
S 0.25c c 0.3c $C_NF 0.3p,black 0.5c Normal ($NNF)
S 0.25c c 0.3c $C_SS 0.3p,black 0.5c Strike-slip ($NSS)
S 0.25c c 0.3c $C_OB 0.3p,black 0.5c Oblique ($NOB)
G 0.08c
N 2
S 0.25c c 0.22c gray30 0.2p,gray10 0.5c Intermediate depth, > 70 km ($NIN)
S 0.3c - 0.45c - 0.6p,gray20 0.65c Quaternary fault (USGS)
EOF
gmt legend leg.txt -DJBC+w17c+o0/0.8c -F+gwhite+p0.4p,gray40

gmt set FONT_ANNOT_PRIMARY 6.5p,Helvetica MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.25p,gray70
gmt ternary -R0/100/0/100/0/100 -JX4.2c -X1.3c -Y-8.3c -Bafg -L-/-/- \
    -B+t"(b) Shallow, $NSH events" tern_sh_TF.txt -Sc0.12c -G$C_TF -W0.15p,black
gmt ternary tern_sh_NF.txt -Sc0.13c -G$C_NF -W0.15p,black
gmt ternary tern_sh_SS.txt -Sc0.13c -G$C_SS -W0.15p,black
gmt ternary tern_sh_OB.txt -Sc0.13c -G$C_OB -W0.15p,black
gmt text -R0/4.2/0/3.637 -JX4.2c/3.637c -N -F+f8p,Helvetica-Bold+j <<'EOF'
0 -0.42 TC Thrust
4.2 -0.42 TC Normal
2.1 4.15 BC Strike-slip
EOF
gmt ternary -R0/100/0/100/0/100 -JX4.2c -X6.4c -Bafg -L-/-/- \
    -B+t"(c) Intermediate depth, $NIN events" tern_in_TF.txt -Sc0.10c -G$C_TF -W0.1p,black
gmt ternary tern_in_NF.txt -Sc0.10c -G$C_NF -W0.1p,black
gmt ternary tern_in_SS.txt -Sc0.10c -G$C_SS -W0.1p,black
gmt ternary tern_in_OB.txt -Sc0.10c -G$C_OB -W0.1p,black
gmt text -R0/4.2/0/3.637 -JX4.2c/3.637c -N -F+f8p,Helvetica-Bold+j <<'EOF'
0 -0.42 TC Thrust
4.2 -0.42 TC Normal
2.1 4.15 BC Strike-slip
EOF
gmt rose paz_sh.txt -R0/$RMAX/0/360 -JX3.8c -X5.9c -Y0.2c -A10 -T -Zu -L-,-,-,N -G$C_NF@15 -W0.4p,black \
    -Bxg5 -Byg30 -B+t"(d) P axes, shallow"
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-13.6c -Y-2.9c <<EOF
0 1.55 Depicts: (a) $NALL Global CMT focal mechanisms, 1976-2025 ($NSH at <= 70 km coloured by class; $NIN intermediate-depth in grey);
0 1.15 (b, c) sin@+2@+ of T-, P- and B-axis plunges; (d) P-axis azimuths in 10-degree sectors (axial).
0 0.75 Data: Global CMT Project (NDK); USGS OFR 2007-1103 faults; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
