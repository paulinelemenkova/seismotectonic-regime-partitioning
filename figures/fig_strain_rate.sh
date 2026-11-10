#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
W=59; E=76; S=28; N=40; R=$W/$E/$S/$N; J=M17c
RB=68.5/74/34.5/39; JB=M8.4c
H_S=30; H_I=230; MU=3.0e10; T=50; NMIN=3
mkdir -p work && cd work

python3 - $H_S $H_I $MU $T $NMIN <<'PY'
import sys, math, urllib.request
HS, HI, MU, T, NMIN = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
B = "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
files = ["jan76_dec20.ndk"]
mon = "jan feb mar apr may jun jul aug sep oct nov dec".split()
for y in range(2021, 2026):
    files += [f"NEW_MONTHLY/{y}/{m}{y % 100:02d}.ndk" for m in mon]
txt = []
for f in files:
    try:
        txt.append(urllib.request.urlopen(urllib.request.Request(B + f, headers={"User-Agent": "afg-paper/1.0"}),
                                          timeout=180).read().decode("latin-1"))
    except Exception as e:
        print("skip", f, e, file=sys.stderr)
L = "".join(txt).splitlines()
ev = []
for i in range(0, len(L) - 4, 5):
    l1, l3, l4 = L[i], L[i + 2], L[i + 3]
    if not l3.startswith("CENTROID") or int(l1[5:9]) > 2025:
        continue
    c = l3.split(); lat, lon, dep = float(c[3]), float(c[5]), float(c[7])
    if not (59 <= lon <= 76 and 28 <= lat <= 40):
        continue
    m = l4.split(); ex = int(m[0])
    mrr, mtt, mpp, mrt, mrp, mtp = (float(m[k]) * 10 ** ex * 1e-7 for k in (1, 3, 5, 7, 9, 11))
    ev.append((lon, lat, dep, mpp, mtt, -mtp))
RE = 6371.0e3

def kostrov(sel, x0, x1, y0, y1, d, H, tag):
    rows = []
    nx, ny = round((x1 - x0) / d), round((y1 - y0) / d)
    for j in range(ny):
        for i in range(nx):
            a, b = x0 + i * d, y0 + j * d
            s = [e for e in sel if a <= e[0] < a + d and b <= e[1] < b + d]
            if len(s) < NMIN:
                continue
            area = RE ** 2 * math.radians(d) * (math.sin(math.radians(b + d)) - math.sin(math.radians(b)))
            k = 1.0 / (2 * MU * area * H * 1e3 * T)
            ee, nn, en = (k * sum(e[q] for e in s) for q in (3, 4, 5))
            mean = 0.5 * (ee + nn); rad = math.sqrt((0.5 * (ee - nn)) ** 2 + en ** 2)
            e1, e2 = mean + rad, mean - rad
            th = 0.5 * math.degrees(math.atan2(2 * en, ee - nn))
            az2 = (90.0 - th + 90.0) % 180
            sc = math.sqrt(ee ** 2 + nn ** 2 + 2 * en ** 2)
            rows.append((a + d / 2, b + d / 2, len(s), ee, nn, en, e1, e2, az2, sc))
    with open(f"cells_{tag}.txt", "w") as f:
        f.write("# lon lat n e_ee e_nn e_en e1 e2 az_e2 scalar (1/yr)\n")
        for r in rows:
            f.write(" ".join(f"{v:.6g}" for v in r) + "\n")
    with open(f"grid_{tag}.xyz", "w") as f:
        for r in rows:
            f.write(f"{r[0]} {r[1]} {math.log10(r[9]):.4f}\n")
    with open(f"cross_{tag}.txt", "w") as f:
        for r in rows:
            mx = max(abs(r[6]), abs(r[7]))
            f.write(f"{r[0]} {r[1]} {r[6] / mx:.3f} {r[7] / mx:.3f} {r[8]:.1f}\n")
    with open(f"hist_{tag}.txt", "w") as f:
        for r in rows:
            f.write(f"{math.log10(r[9]):.4f}\n")
    lg = sorted(math.log10(r[9]) for r in rows)
    top = max(rows, key=lambda r: r[9])
    comp = sum(1 for r in rows if r[7] < 0 and abs(r[7]) >= abs(r[6]))
    print(tag, "cells", len(rows), "events", sum(r[2] for r in rows),
          "log10 min/med/max", f"{lg[0]:.2f} {lg[len(lg) // 2]:.2f} {lg[-1]:.2f}",
          "max cell", f"{top[0]},{top[1]} n={top[2]} rate={top[9]:.3e} e2az={top[8]:.0f}",
          "shortening-dominated", comp)
    return rows

sh = [e for e in ev if e[2] <= 70]
im = [e for e in ev if 70 < e[2] <= 300]
rs = kostrov(sh, 59, 76, 28, 40, 1.0, HS, "sh")
ri = kostrov(im, 68.5, 74, 34.5, 39, 0.5, HI, "in")
open("n.txt", "w").write(f"{len(sh)} {len(im)} {len(rs)} {len(ri)}\n")
PY
read NSH NIN CSH CIN < n.txt
mkdir -p "$OUT"; cp cells_sh.txt "$OUT/kostrov_cells_shallow.txt"; cp cells_in.txt "$OUT/kostrov_cells_intermediate.txt"

gmt xyz2grd grid_sh.xyz -R$R -I1 -r -Gsh.nc
gmt xyz2grd grid_in.xyz -R$RB -I0.5 -r -Gin.nc
gmt makecpt -Cjet -T-9.5/-6.5/0.25 -D > rate.cpt
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt
gmt coast -EAF -M -Di > af.txt

gmt begin $OUT/fig_strain_rate pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 9p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray40 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y15.5c
gmt grdimage sh.nc -Crate.cpt -Q -t20
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot af.txt -W1.4p,black
gmt velo cross_sh.txt -Sx0.62c -W3.4p,white -Gwhite -A0.2c+e+a50
gmt velo cross_sh.txt -Sx0.62c -W1.4p,black -Gblack -A0.18c+e+a50
printf "59.15 39.85 (a) Shallow events (<= 70 km), 1 x 1 degree cells\n" | \
  gmt text -F+f9p,Helvetica-Bold+jTL -Gwhite@15 -C0.06c/0.05c
echo "$RB" | awk -F/ '{printf "%s %s\n%s %s\n%s %s\n%s %s\n%s %s\n",$1,$3,$2,$3,$2,$4,$1,$4,$1,$3}' | gmt plot -W1p,black,-
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.1c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
gmt colorbar -DJBC+w12c/0.3c+o0/0.75c+h+ef -Crate.cpt -Bxa0.5f0.25+l"log@-10@- horizontal scalar strain rate (yr@+-1@+)"
gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 7p,Helvetica
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$RB -J$JB -Y-10.8c
gmt grdimage in.nc -Crate.cpt -Q -t20
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa1f0.5g0.5 -Bya1f0.5g0.5 -BWSen
gmt plot af.txt -W1.4p,black
gmt velo cross_in.txt -Sx0.48c -W3.0p,white -Gwhite -A0.17c+e+a50
gmt velo cross_in.txt -Sx0.48c -W1.2p,black -Gblack -A0.15c+e+a50
printf "68.6 38.93 (b) 70-300 km, 0.5 x 0.5 degree cells\n" | \
  gmt text -F+f8p,Helvetica-Bold+jTL -Gwhite@15 -C0.05c/0.04c
gmt histogram hist_sh.txt -R-10.5/-6/0/14 -JX6.6c/5.8c -X10.2c -Y1.4c -T0.25 -Z0 -G220/20/60@30 -W0.6p,220/20/60 \
    -Bxa1f0.25+l"log@-10@- scalar strain rate (yr@+-1@+)" -Bya2f1+l"Cells" -BWSen+t"(c) Cell strain rates"
gmt histogram hist_in.txt -T0.25 -Z0 -G0/70/255@30 -W0.6p,0/70/255
cat > leg.txt <<EOF
S 0.2c r 0.3c 220/20/60@30 0.6p,220/20/60 0.45c Shallow ($CSH cells)
S 0.2c r 0.3c 0/70/255@30 0.6p,0/70/255 0.45c Intermediate ($CIN cells)
EOF
gmt legend leg.txt -DjTL+w3.4c+o0.15c -F+gwhite+p0.3p,gray50 --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-10.2c -Y-4.4c <<EOF
0 1.55 Depicts: Kostrov summation of $NSH shallow and $NIN intermediate-depth Global CMT moment tensors, 1976-2025 (cells with >= $NMIN events);
0 1.15 layer thickness $H_S km (a) and $H_I km (b), shear modulus 30 GPa; crosses: normalised principal horizontal strain rates (inward arrows = shortening).
0 0.75 Data: Global CMT Project (NDK); SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
