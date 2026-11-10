#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
BG=${BG:-bgrid_mb.csv}
RA=62/76/28/41; JA=M17c
RB=68/76/34.5/40; JB=M8.6c
mkdir -p work && cp "$BG" work/bg.csv && cd work

awk -F, 'NR>1 && $3=="crust"{print $1,$2,$4 > "b_cr.xyz"; print $1,$2,$5 > "s_cr.xyz"}
         NR>1 && $3=="deep" {print $1,$2,$4 > "b_dp.xyz"; print $1,$2,$5 > "s_dp.xyz"}
         NR>1 && $3=="crust"{print $4 > "h_cr.txt"} NR>1 && $3=="deep"{print $4 > "h_dp.txt"}' bg.csv
NCR=$(wc -l < b_cr.xyz); NDP=$(wc -l < b_dp.xyz)
MCR=$(sort -n h_cr.txt | awk '{a[NR]=$1} END{printf "%.2f", (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}')
MDP=$(sort -n h_dp.txt | awk '{a[NR]=$1} END{printf "%.2f", (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}')
for L in cr dp; do
  gmt xyz2grd b_$L.xyz -R58.25/75.75/27.25/40.75 -I0.25 -Gb_$L.nc
  gmt xyz2grd s_$L.xyz -R58.25/75.75/27.25/40.75 -I0.25 -Gs_$L.nc
done
gmt makecpt -Cjet -T1.2/2.6/0.1 -D > b.cpt

curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
ogr2ogr -f OGR_GMT qf.gmt "$(find . -name 'faults.shp' | head -1)" -t_srs EPSG:4326
gmt grdcut @earth_relief_30s -R58/77/27/42 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt
gmt coast -EAF -M -Di > af.txt

gmt begin $OUT/fig_bvalue_map pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 9p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray40 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$RA -J$JA -Y15.2c
gmt grdimage b_cr.nc -Cb.cpt -Q -t25 -nn
gmt grdcontour s_cr.nc -C0.1 -A0.1+f6p,Helvetica,black+gwhite@20 -W0.6p,black,4_2:0p -S4
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot qf.gmt -W0.35p,gray10
gmt plot af.txt -W1.4p,black
printf "62.15 40.85 (a) Crustal events (< 70 km), median b = $MCR\n" | \
  gmt text -F+f9p,Helvetica-Bold+jTL -Gwhite@15 -C0.06c/0.05c
echo "$RB" | awk -F/ '{printf "%s %s\n%s %s\n%s %s\n%s %s\n%s %s\n",$1,$3,$2,$3,$2,$4,$1,$4,$1,$3}' | gmt plot -W1p,black,-
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.1c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
gmt colorbar -DJBC+w12c/0.3c+o0/0.75c+h+e -Cb.cpt -Bxa0.2f0.1+l"b-value (body-wave magnitude, declustered, 1973-2025)"
gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 7p,Helvetica
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$RB -J$JB -Y-10.9c
gmt grdimage b_dp.nc -Cb.cpt -Q -t25 -nn
gmt grdcontour s_dp.nc -C0.1 -A0.1+f6p,Helvetica,black+gwhite@20 -W0.6p,black,4_2:0p -S4
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya1f0.5g0.5 -BWSen
gmt plot af.txt -W1.4p,black
printf "68.1 39.9 (b) >= 70 km, median b = $MDP\n" | gmt text -F+f8p,Helvetica-Bold+jTL -Gwhite@15 -C0.05c/0.04c
gmt histogram h_cr.txt -R1.0/3.0/0/14 -JX6.4c/5.6c -X10.4c -Y1.6c -T0.05 -Z1 -G220/20/60@35 -W0.6p,220/20/60 \
    -Bxa0.5f0.1+l"b-value" -Bya2f1+l"Nodes (%)" -BWSen+t"(c) Node b-values"
gmt histogram h_dp.txt -T0.05 -Z1 -G0/70/255@35 -W0.6p,0/70/255
printf "$MCR 0\n$MCR 14\n" | gmt plot -W0.9p,220/20/60,6_3:0p
printf "$MDP 0\n$MDP 14\n" | gmt plot -W0.9p,0/70/255,6_3:0p
cat > leg.txt <<EOF
S 0.2c r 0.3c 220/20/60@35 0.6p,220/20/60 0.45c Crustal ($NCR nodes)
S 0.2c r 0.3c 0/70/255@35 0.6p,0/70/255 0.45c Intermediate ($NDP nodes)
S 0.2c - 0.35c - 0.9p,gray20,6_3:0p 0.45c Median
EOF
gmt legend leg.txt -DjTR+w3.4c+o0.15c -F+gwhite+p0.3p,gray50 --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-10.4c -Y-4.6c <<EOF
0 1.55 Depicts: b-values at 0.25-degree nodes from the nearest 150 background events within 150 km (>= 50 events above Mc, Mc = MAXC + 0.2);
0 1.15 dashed contours: Shi-Bolt standard error sigma@-b@- (0.1 interval); black lines: USGS Quaternary faults (a).
0 0.75 Data: USGS ComCat and ISC-GEM merged catalogue (mb, declustered); USGS OFR 2007-1103; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
