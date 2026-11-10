#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
W=58; E=76; S=27; N=41; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work

B="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=$S&maxlatitude=$N&minlongitude=$W&maxlongitude=$E&minmagnitude=4&orderby=magnitude-asc"
curl -s -m 150 "$B&starttime=1900-01-01&endtime=1999-12-31T23:59:59" >  c1.csv
curl -s -m 150 "$B&starttime=2000-01-01&endtime=2025-12-31T23:59:59" >  c2.csv
( tail -n +2 c1.csv; tail -n +2 c2.csv ) | awk -F, '$4!=""{printf "%s %s %s %s %s\n",$3,$2,$4,$5,substr($1,1,10)}' \
  | sort -k4,4g > eq.txt
NEQ=$(wc -l < eq.txt)
awk '{print $3}' eq.txt > depth.txt
awk '{print $4}' eq.txt > mag.txt
awk '$5=="2022-06-21"&&$4>=5.8 || $5=="2023-10-07"&&$4>=6.2 || $5=="2025-08-31"&&$4>=5.9 || $5=="2025-11-02"&&$4>=6.2' eq.txt | sort -u -k5,5 > key.txt
cat key.txt

curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
ogr2ogr -f OGR_GMT qf.gmt faults.shp -t_srs EPSG:4326
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt coast -EAF -M -Di > af.txt

gmt makecpt -Cgray -T-14000/6000/250 > relief.cpt
gmt makecpt -Cturbo -T0/300/10 -I -D > depth.cpt

gmt begin $OUT/fig_catalogue_overview pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 10p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray55 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J
gmt coast -Di -W0.3p,gray30 -N1/0.6p,gray25,- -I1/0.35p,steelblue -C+l+gsteelblue -A50 \
          -Bxa2f1g2 -Bya2f1g2 -BWSen
gmt plot qf.gmt -W0.3p,gray20
gmt plot af.txt -W1.4p,black
awk '{printf "%s %s %s %.3f\n",$1,$2,$3,0.03*exp(($4-3)*log(1.9))}' eq.txt \
  | gmt plot -Sc -Cdepth.cpt -W0.1p,black@40
awk '{print $1,$2}' key.txt | gmt plot -Sa0.55c -W1.2p,black -Gwhite@100
gmt text -F+f8p,Helvetica-Bold,black+j -Gwhite@20 -C0.05c/0.04c -N <<'EOF'
EOF
gmt text -F+f8p,Helvetica-Oblique,gray15+j -Gwhite@45 -C0.04c/0.03c <<'EOF'
59.2 32.5 CM Iran
72.6 29.4 CM Pakistan
61.5 39.0 CM Turkmenistan
66.0 40.2 CM Uzbekistan
69.4 39.7 CM Tajikistan
75.0 37.6 CM China
EOF
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.2c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<EOF
H 9p,Helvetica-Bold Magnitude
N 4
S 0.3c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^1}')c gray70 0.3p,black 0.6c M4
S 0.3c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^2}')c gray70 0.3p,black 0.6c M5
S 0.3c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^3}')c gray70 0.3p,black 0.6c M6
S 0.3c c $(awk 'BEGIN{printf "%.3f",0.03*1.9^4}')c gray70 0.3p,black 0.6c M7
N 1
S 0.3c a 0.4c white 1.0p,black 0.6c Destructive earthquake, 2022-2025
EOF
gmt legend leg.txt -DjBL+w6.2c+o0.3c/0.3c -F+gwhite@10+p0.4p,gray40+r0.1c

gmt colorbar -Cdepth.cpt -DJBC+w17c/0.32c+o0/1.0c+h+ef -Bxa50f10+l"Hypocentral depth (km)"

gmt histogram depth.txt -R0/300/1/10000 -JX7.9c/3.4cl -T10 -Gsteelblue3 -W0.3p,white \
   -Bxa50f10+l"Depth (km)" -Bya1pf3+l"Events" -BWSne -Y-6.0c
gmt text -F+cTL+f10p,Helvetica-Bold+t"(a)" -D0.15c/-0.15c
gmt histogram mag.txt -R3.9/8/1/10000 -JX7.9c/3.4cl -T0.1 -Gfirebrick3 -W0.3p,white \
   -Bxa0.5f0.1+l"Magnitude" -Bya1pf3+l"Events" -BWSne -X9.1c
gmt text -F+cTL+f10p,Helvetica-Bold+t"(b)" -D0.15c/-0.15c

gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -X-9.1c -Y-2.9c <<EOF
0 1.55 Depicts: earthquakes M@~\263@~4, 1900-2025 (n = $NEQ), coloured by hypocentral depth and sized by magnitude;
0 1.15 Quaternary faults; destructive events of 2022-2025 (stars); (a) depth and (b) magnitude histograms, logarithmic count axis.
0 0.75 Data: USGS ComCat; USGS OFR 2007-1103; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version). Source: authors.
EOF
gmt end
cd ..
