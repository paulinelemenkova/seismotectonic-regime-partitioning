#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
W=58; E=76; S=27; N=41; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work
python3 - <<'PY'
import json, os, urllib.request
it = []
for p in (1, 2, 3):
    u = ("https://www.ngdc.noaa.gov/hazel/hazard-service/api/v1/earthquakes?minYear=1900&maxYear=2025"
         f"&minLatitude=27&maxLatitude=41&minLongitude=58&maxLongitude=76&page={p}")
    d = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "afg/1.0"}), timeout=100))
    it += d["items"]
    if p >= d["totalPages"]:
        break
def cls(e):
    n = e.get("deathsTotal") or e.get("deaths")
    o = e.get("deathsAmountOrderTotal") or e.get("deathsAmountOrder")
    if n:
        return 1 if n <= 100 else 2 if n <= 1000 else 3 if n <= 10000 else 4, n
    if o:
        return {1: 1, 2: 1, 3: 2, 4: 3}[o], None
    return 0, None
rows, tl = [], []
for e in it:
    if e.get("latitude") is None:
        continue
    k, n = cls(e)
    m = e.get("eqMagnitude") or e.get("eqMagMs") or 5.0
    rows.append((m, e["longitude"], e["latitude"], k, 0.03 * 1.9 ** (m - 3)))
    if n:
        afg = os.environ.get("TARGET_TAG", "\x00") in (e.get("locationName") or "")
        t = e["year"] + ((e.get("month") or 6) - 0.5) / 12
        tl.append((t, n, k, afg))
rows.sort()
open("eq.txt", "w").write("".join(f"{x[1]} {x[2]} {x[3]} {x[4]:.3f}\n" for x in rows))
open("tl_afg.txt", "w").write("".join(f"{t:.2f} {n} {k}\n" for t, n, k, a in tl if a))
open("tl_oth.txt", "w").write("".join(f"{t:.2f} {n} {k}\n" for t, n, k, a in tl if not a))
print("events", len(rows), "with deaths", len(tl), "afg", sum(a for *_, a in tl),
      "total deaths", sum(n for _, n, _, _ in tl), "afg deaths", sum(n for _, n, _, a in tl if a))
open("n.txt", "w").write(f"{len(rows)} {len(tl)} {sum(a for *_, a in tl)}\n")
PY
read NEV NDE NAF < n.txt
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt coast -EAF -M -Di > af.txt
gmt makecpt -Cgray -T-16000/6000/250 > relief.cpt
cat > fat.cpt <<'EOF'
-0.5 200/200/200 0.5 200/200/200 ;No fatalities reported
0.5 255/221/87 1.5 255/221/87 ;1-100
1.5 244/134/40 2.5 244/134/40 ;101-1000
2.5 205/30/38 3.5 205/30/38 ;1001-10 000
3.5 103/0/31 4.5 103/0/31 ;> 10 000
EOF

gmt begin $OUT/fig_damaging_events pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        MAP_GRID_PEN_PRIMARY 0.25p,gray55 FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y14.0c
gmt coast -Di -W0.3p,gray30 -N1/0.6p,gray25,- -I1/0.35p,steelblue -C+l+gsteelblue -A50 -Bxa2f1g2 -Bya2f1g2 -BWSen
gmt plot af.txt -W1.4p,black
gmt plot eq.txt -Sc -Cfat.cpt -W0.4p,black
cat > lab.txt <<'EOF'
58.364 38.084 58.70 39.00 ML Ashgabat 1948: 110 000
73.649 34.451 74.25 33.55 TC Kashmir 2005: 76 213
66.80 29.50 65.90 28.55 ML Quetta 1935: 60 000
58.311 28.995 58.70 28.40 ML Bam 2003: 31 000
67.90 38.50 67.60 39.55 MR Karatag 1907: 12 000
59.083 34.158 59.40 33.25 ML Dasht-e Bayaz 1968: 10 488
69.315 36.062 68.60 35.40 MR Nahrin 2002: 1000
EOF
awk '{printf ">\n%s %s\n%s %s\n",$1,$2,$3,$4}' lab.txt | gmt plot -W0.6p,gray15
awk '{print $1,$2}' lab.txt | gmt plot -Sc0.08c -Gblack
awk '{printf "%s %s %s", $3, $4, $5; for (i = 6; i <= NF; i++) printf " %s", $i; print ""}' lab.txt | \
  gmt text -F+f8p,Helvetica-Bold,black+j -Gwhite@15 -C0.05c/0.04c -W0.3p,gray30 -N
printf "58.25 40.75 (a)\n" | gmt text -F+f11p,Helvetica-Bold+jTL -Gwhite@20 -C0.06c/0.05c
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.2c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<EOF
H 9p,Helvetica-Bold Reported fatalities
N 5
S 0.25c c 0.3c 200/200/200 0.4p,black 0.5c None reported
S 0.25c c 0.3c 255/221/87 0.4p,black 0.5c 1-100
S 0.25c c 0.3c 244/134/40 0.4p,black 0.5c 101-1000
S 0.25c c 0.3c 205/30/38 0.4p,black 0.5c 1001-10 000
S 0.25c c 0.3c 103/0/31 0.4p,black 0.5c > 10 000
G 0.1c
H 9p,Helvetica-Bold Magnitude
N 4
S 0.45c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^2}')c white 0.4p,black 1.0c M5
S 0.45c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^3}')c white 0.4p,black 1.0c M6
S 0.45c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^4}')c white 0.4p,black 1.0c M7
S 0.45c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^4.6}')c white 0.4p,black 1.0c M7.6
EOF
gmt legend leg.txt -DJBC+w17c+o0/0.8c -F+gwhite+p0.4p,gray40
gmt basemap -R1900/2030/1/300000 -JX17c/4.6cl -Bxa10f5+l"Year" -Bya1pf3+l"Fatalities per event" -BWSen -Y-8.4c
awk '{print $1, $2, $3}' tl_oth.txt | gmt plot -Sb0.12cb1 -Cfat.cpt -W0.2p,gray40 -t45
awk '{print $1, $2, $3}' tl_afg.txt | gmt plot -Sb0.12cb1 -Cfat.cpt -W0.5p,black
gmt plot tl_afg.txt -Sc0.13c -Gblack
cat > tlab.txt <<'EOF'
1998.1 2323 7000 BR Rostaq 1998
2002.2 1000 1800 BL Nahrin 2002
EOF
awk '{printf ">\n%s %s\n%s %s\n",$1,$2,$1,$3}' tlab.txt | gmt plot -W0.5p,gray15
awk '{x = ($4 == "BR") ? $1 - 0.3 : $1 + 0.3; printf "%s %s %s", x, $3, $4; for (i = 5; i <= NF; i++) printf " %s", $i; print ""}' tlab.txt | \
  gmt text -F+f7.5p,Helvetica,black+j -Gwhite@15 -C0.04c/0.03c -N
printf "1902 200000 (b)\n" | gmt text -F+f11p,Helvetica-Bold+jTL -Gwhite@20 -C0.06c/0.05c
cat > leg2.txt <<'EOF'
N 2
S 0.3c r 0.12c/0.35c 244/134/40 0.5p,black 0.7c Target area: opaque bar, black outline, dot at top
S 0.3c r 0.12c/0.35c 250/196/150 0.2p,gray60 0.7c Neighbouring countries: translucent bar
EOF
gmt legend leg2.txt -DJBC+w17c+o0/1.05c -F+gwhite+p0.3p,gray50
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -Y-3.7c <<EOF
0 1.55 Depicts: (a) $NEV damaging earthquakes in 58-76E, 27-41N, 1900-2025, coloured by reported fatalities and sized by magnitude;
0 1.15 (b) fatalities per event for the $NDE events with reported deaths ($NAF in the target area), logarithmic axis.
0 0.75 Data: NOAA NCEI/WDS Global Significant Earthquake Database; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
