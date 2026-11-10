#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
CELLS=${CELLS:-moment_cells.xyz}; ANN=${ANN:-moment_annual.txt}
R=58/76/27/41; J=M17c
mkdir -p work && cp "$CELLS" work/cells.xyz && cp "$ANN" work/ann.txt && cd work
NC=$(wc -l < cells.xyz)

gmt xyz2grd cells.xyz -R58/76/27/41 -I0.5 -r -Gm0.nc
gmt makecpt -Cjet -T16/21.5/0.25 -D > m0.cpt
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt
gmt coast -EAF -M -Di > af.txt
cat > recent.txt <<'EOF'
EOF

C_HKP=0/70/255; C_NEC=0/175/200; C_PBZ=220/20/60; C_WEST=240/160/0

gmt begin $OUT/fig_moment_release pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        MAP_GRID_PEN_PRIMARY 0.25p,gray40 FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y13.6c
gmt grdimage m0.nc -Cm0.cpt -Q -t20 -nn
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot af.txt -W1.4p,black
gmt plot recent.txt -Sa0.42c -Gwhite -W0.9p,black
awk '{print $1, $2, $3" "$4}' recent.txt | gmt text -F+f8p,Helvetica-Bold+jLM -D0.28c/0 -Gwhite@15 -C0.04c/0.03c
printf "58.2 40.85 (a) Cumulative seismic moment per 0.5-degree cell, 1900-2025\n" | \
  gmt text -F+f9p,Helvetica-Bold+jTL -Gwhite@15 -C0.06c/0.05c
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.1c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
gmt colorbar -DJBC+w12c/0.3c+o0/0.75c+h+e -Cm0.cpt -Bxa0.5f0.25+l"log@-10@- cumulative seismic moment M@-0@- (N m)"
gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 7.5p,Helvetica MAP_GRID_PEN_PRIMARY 0.25p,gray80
gmt basemap -R1900/2026/1e18/1e23 -JX15.6c/6.8cl -X1.0c -Y-9.0c \
    -Bxa10f5g10+l"Year" -Bya1pf3g1p+l"Cumulative seismic moment M@-0@- (N m)" -BWSen
i=2
for d in HKP NEC PBZ WEST; do
  eval col=\$C_$d
  awk -v c=$i '$c>0{print $1, $c}' ann.txt | gmt plot -W1.6p,$col
  i=$((i+1))
done
awk '$1==2022{print 2022.5, $4, "E1"} $1==2023{print 2023.8, $5, "E2"} $1==2025{print 2025.7, $4, "E3"; print 2025.85, $3, "E4"}' ann.txt > rec_b.txt
gmt plot rec_b.txt -Sa0.26c -Gwhite -W0.7p,black
cat > lab_b.txt <<'EOF'
1921.9 HKP 1
1935.4 PBZ 2
2005.8 PBZ 3
2013.3 WEST 4
EOF
while read yr dom k; do
  col=$(awk -v d=$dom 'BEGIN{split("HKP NEC PBZ WEST",a," "); for(i in a) if(a[i]==d) print i+1}')
  v=$(awk -v y=${yr%.*} -v c=$col '$1==y{print $c}' ann.txt)
  echo "$yr $v $k"
done < lab_b.txt > num_b.txt
gmt plot num_b.txt -Sc0.36c -Gwhite -W0.8p,black
gmt text num_b.txt -F+f7p,Helvetica-Bold,black+jCM
printf "1902 1.3e18 (b) Cumulative moment by domain\n" | gmt text -F+f9p,Helvetica-Bold+jBL -Gwhite@15 -C0.05c/0.04c
cat > key.txt <<'EOF'
L 7p,Helvetica L (1) 1921, M@-w@- 7.8
L 7p,Helvetica L (2) Quetta 1935, M@-w@- 7.6
L 7p,Helvetica L (3) Kashmir 2005, M@-w@- 7.6
L 7p,Helvetica L (4) Saravan 2013, M@-w@- 7.7
EOF
gmt legend key.txt -DjTL+w4.2c+o0.15c -F+gwhite+p0.3p,gray50
cat > leg.txt <<EOF
N 1
S 0.25c - 0.5c - 1.6p,$C_HKP 0.7c HKP intermediate-depth domain
S 0.25c - 0.5c - 1.6p,$C_NEC 0.7c NEC north-eastern crust
S 0.25c - 0.5c - 1.6p,$C_PBZ 0.7c PBZ plate-boundary zone
S 0.25c - 0.5c - 1.6p,$C_WEST 0.7c WEST western interior
S 0.25c a 0.3c white 0.8p,black 0.7c 2022-2025 destructive events
EOF
gmt legend leg.txt -DjBR+w6.4c+o0.15c -F+gwhite+p0.3p,gray50 --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-1.0c -Y-3.4c <<EOF
0 1.55 Depicts: (a) cumulative moment of all events in $NC cells, M@-0@- = 10@+1.5 M@-w@- + 9.1@+ N m, 1900-2025, with the 2022-2025 destructive earthquakes;
0 1.15 (b) cumulative moment of each seismotectonic domain through time (logarithmic axis), with the largest events annotated.
0 0.75 Data: merged USGS ComCat and ISC-GEM catalogue (M@-w@-); SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
