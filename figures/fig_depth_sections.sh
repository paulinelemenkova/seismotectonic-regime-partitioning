#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
mkdir -p work && cd work
HW=50

python3 - <<'PY'
import csv, io, urllib.request
u = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=27&maxlatitude=42"
     "&minlongitude=58&maxlongitude=78&minmagnitude=4&starttime=1900-01-01&endtime=2025-12-31T23:59:59"
     "&orderby=magnitude")
r = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "afg-paper/1.0"}), timeout=300)
rows = list(csv.DictReader(io.TextIOWrapper(r, encoding="utf-8")))
n_all = len(rows); keep = []
for e in rows:
    try:
        d = float(e["depth"]); m = float(e["mag"])
    except ValueError:
        continue
    if round(d, 1) in (10.0, 33.0, 35.0):
        continue
    keep.append((float(e["longitude"]), float(e["latitude"]), d, m))
keep.sort(key=lambda x: -x[3])
with open("eq.txt", "w") as f:
    for k in keep:
        f.write(f"{k[0]} {k[1]} {k[2]} {k[3]}\n")
open("n.txt", "w").write(f"{n_all} {len(keep)}\n")
print("events", n_all, "kept (non-default depth)", len(keep))
PY
read NALL NKEEP < n.txt

cat > lines.txt <<'EOF'
69.6 39.2 71.6 34.6 A
75.5 40.0 72.8 36.2 B
65.5 32.6 70.8 30.4 C
62.0 34.3 68.4 37.0 D
EOF

gmt grdcut @earth_relief_30s -R58/78/27/42 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt
gmt makecpt -Cjet -T0/300/25 -D -I > dep.cpt
gmt coast -EAF -M -Di > af.txt

while read x0 y0 x1 y1 L; do
  gmt project eq.txt -C$x0/$y0 -E$x1/$y1 -W-$HW/$HW -Lw -Q -Fpz > sec_$L.txt
  gmt project -C$x0/$y0 -E$x1/$y1 -G1 -Q > pts_$L.txt
  gmt grdtrack pts_$L.txt -Grelief.nc | awk '{print $3, $4/1000}' > topo_$L.txt
  LEN=$(gmt project -C$x0/$y0 -E$x1/$y1 -G1000 -Q | tail -1 | awk '{printf "%.0f", $3}')
  echo "$L $LEN $(wc -l < sec_$L.txt)" >> len.txt
done < lines.txt
cat len.txt

msize() { awk '{s=0.035*1.55^($3-3); print $1, $2, $2, s}'; }

gmt begin $OUT/fig_depth_sections pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7.5p,Helvetica FONT_LABEL 8.5p,Helvetica FONT_TITLE 9p,Helvetica-Bold \
        MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.25p,gray60 MAP_TICK_LENGTH_PRIMARY 3p \
        FORMAT_GEO_MAP ddd:mmF
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R60/77/29/41 -JM9.6c -Y30.0c
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa4f1g2 -Bya2f1g2 -BWSen \
          --MAP_FRAME_TYPE=fancy --MAP_FRAME_WIDTH=2.5p
gmt plot af.txt -W1.2p,black
awk '{print $1, $2, $3, 0.035*1.55^($4-3)}' eq.txt | gmt plot -Sc -Cdep.cpt -W0.15p,gray20
while read x0 y0 x1 y1 L; do
  printf "%s %s\n%s %s\n" $x0 $y0 $x1 $y1 | gmt plot -W2.6p,white
  printf "%s %s\n%s %s\n" $x0 $y0 $x1 $y1 | gmt plot -W1.3p,black
  printf "%s %s %s\n%s %s %s'\n" $x0 $y0 $L $x1 $y1 $L | \
    gmt text -F+f9p,Helvetica-Bold,black+jCM -Gwhite -W0.5p,black -C0.05c/0.04c
done < lines.txt
printf "60.2 40.85 (a)\n" | gmt text -F+f10p,Helvetica-Bold+jTL -Gwhite@15 -C0.05c/0.04c
gmt basemap -LjBR+c35+w200k+f+l"km"+o0.3c/0.4c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt colorbar -Cdep.cpt -DJMR+w6.0c/0.32c+o0.6c/0+e -Bxa50f25+l"Hypocentral depth (km)"
cat > leg.txt <<EOF
H 8.5p,Helvetica-Bold Magnitude
S 0.35c c $(awk 'BEGIN{printf "%.3f", 0.035*1.55^1}')c white 0.4p,black 0.8c M4
S 0.35c c $(awk 'BEGIN{printf "%.3f", 0.035*1.55^2}')c white 0.4p,black 0.8c M5
S 0.35c c $(awk 'BEGIN{printf "%.3f", 0.035*1.55^3}')c white 0.4p,black 0.8c M6
S 0.35c c $(awk 'BEGIN{printf "%.3f", 0.035*1.55^4}')c white 0.4p,black 0.8c M7
G 0.15c
L 7p,Helvetica L Swath: $HW km either side
L 7p,Helvetica L $NKEEP of $NALL events
L 7p,Helvetica L (default depths excluded)
EOF
gmt legend leg.txt -DJMR+w3.7c+o2.7c/0 -F+gwhite+p0.4p,gray50

sec() {
  local L=$1 ZM=$2 H=$3 YS=$4 P=$5
  local LEN=$(awk -v l=$L '$1==l{print $2}' len.txt)
  local N=$(awk -v l=$L '$1==l{print $3}' len.txt)
  local HD=$(awk -v h=$H 'BEGIN{print h + 0.05}')
  gmt basemap -R0/$LEN/-1/7 -JX16.2c/1.0c -Y$YS -Bya3f1+l"km" -BWsen
  awk 'NR==1{print $1, -1} {print} END{print $1, -1}' topo_$L.txt | gmt plot -G150/150/150 -W0.4p,gray20
  gmt basemap -R0/$LEN/0/$ZM -JX16.2c/-${H}c -Y-${HD}c \
    -Bxa100f25+l"Distance along section (km)" -Bya50f25g50+l"Depth (km)" -BWSen
  msize < sec_$L.txt | gmt plot -Sc -Cdep.cpt -W0.15p,gray20
  printf "0 0 $L\n$LEN 0 $L'\n" | gmt text -F+f10p,Helvetica-Bold+jCB -D0/1.15c -N
  echo "($P) $L-$L', $N hypocentres" | gmt text -F+cBR+f8p,Helvetica-Bold -Dj0.12c/0.12c -Gwhite@15 -C0.05c/0.03c
}
sec A 300 4.6 -2.5c b
sec B 300 4.6 -2.9c c
sec C 120 2.3 -2.9c d
sec D 120 2.3 -2.9c e
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-0.0c -Y-3.2c <<EOF
0 1.55 Depicts: (a) $NKEEP USGS ComCat hypocentres (M >= 4, 1900-2025; default depths of 10, 33 and 35 km excluded) and section lines;
0 1.15 (b-e) hypocentres within $HW km of each line coloured by depth and sized by magnitude; SRTM15+ topography above each section.
0 0.75 Data: USGS ComCat; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
