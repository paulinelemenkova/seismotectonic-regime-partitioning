#!/usr/bin/env bash
set -e
HALO='=~1.2p,white'
OUT=${OUT:-$(pwd)}
W=55; E=80; S=25; N=42; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work

curl -s -m 120 -L -o pb_b.txt http://peterbird.name/oldFTP/PB2002/PB2002_boundaries.dig.txt
curl -s -m 120 -L -o pb_o.txt http://peterbird.name/oldFTP/PB2002/PB2002_orogens.dig.txt
awk -F, '/^\*\*\*/{next} /^[^ +-]/{print ">"; next} {gsub(/ /,""); printf "%.4f %.4f\n",$1,$2}' pb_b.txt > pb_bound.gmt
awk -F, 'BEGIN{k=0} /^\*\*\*/{k=0;next} /^Persia-Tibet-Burma/{k=1;print ">";next} /^[^ +-]/{k=0;next} k==1{gsub(/ /,""); printf "%.4f %.4f\n",$1,$2}' pb_o.txt > pb_ptb.gmt

curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
for cls in A B C; do
  ogr2ogr -f OGR_GMT qf_R_$cls.gmt  faults.shp -t_srs EPSG:4326 -where "Type='R' AND ActvCls LIKE '%($cls)%'"
  ogr2ogr -f OGR_GMT qf_SS_$cls.gmt faults.shp -t_srs EPSG:4326 -where "Type IN ('RL','LL','SS') AND ActvCls LIKE '%($cls)%'"
  ogr2ogr -f OGR_GMT qf_N_$cls.gmt  faults.shp -t_srs EPSG:4326 -where "Type='N' AND ActvCls LIKE '%($cls)%'"
done

Q="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=$S&maxlatitude=$N&minlongitude=$W&maxlongitude=$E&minmagnitude=4&mindepth=70&starttime=1900-01-01&endtime=2025-12-31T23:59:59"
curl -s -m 120 "$Q" | awk -F, 'NR>1{printf "%s %s %s %s\n",$3,$2,$4,$5}' > deep.txt
NDEEP=$(wc -l < deep.txt)

gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.7 -Gint.nc
gmt coast -EAF -M -Di > af.txt

cat > cities.txt <<'EOF'
EOF
cat > cities_other.txt <<'EOF'
73.081 33.689 Islamabad
71.533 34.007 Peshawar
67.023 30.222 Quetta
68.774 38.560 Dushanbe
59.568 36.272 Mashhad
74.344 31.549 Lahore
69.240 41.311 Tashkent
60.830 29.500 Zahedan
58.383 37.950 Ashgabat
EOF

gmt makecpt -Cgeo -T-4000/7000/250 -Z > topo.cpt

gmt begin $OUT/fig_tectonic_setting pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 11p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.14p,gray50 \
        MAP_TICK_LENGTH_PRIMARY 3p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel

gmt grdimage relief.nc -Iint.nc -Ctopo.cpt -R$R -J$J -t45
gmt clip af.txt
gmt grdimage relief.nc -Iint.nc -Ctopo.cpt
gmt clip -C
gmt coast -Di -W0.16p,gray25 -N1/0.33p,gray20,- -I1/0.22p,dodgerblue3 -C+l+gdodgerblue3 -A50 \
          -Bxa5f1g5 -Bya5f1g5 -BWSen
gmt plot af.txt -W0.76p,black

gmt plot pb_ptb.gmt -W0.49p,mediumpurple3,- -L
gmt plot pb_bound.gmt -W1.2p,purple4

for t in "R firebrick3" "SS royalblue3" "N forestgreen"; do
  set -- $t
  gmt plot qf_$1_C.gmt -W0.14p,$2
  gmt plot qf_$1_B.gmt -W0.33p,$2
  gmt plot qf_$1_A.gmt -W0.76p,$2
done

gmt plot deep.txt -Sc0.07c -Gmagenta3 -W0.05p,black

printf "70.5 27.0 9 1.65c\n"  | gmt plot -SV0.38c+e+a40+h0.3 -W0.87p,black -Gblack
printf "59.2 25.25 14.6 1.62c\n" | gmt plot -SV0.38c+e+a40+h0.3 -W0.87p,black -Gblack
printf "70.8 27.25 41 mm/yr\n" | gmt text -F+f8p,Helvetica-Bold,black${HALO}+jML
printf "59.65 26.05 40 mm/yr\n" | gmt text -F+f8p,Helvetica-Bold,black${HALO}+jML

gmt plot cities.txt -Ss0.2c -Gred3 -W0.33p,white
gmt plot cities_other.txt -Sc0.14c -Gwhite -W0.27p,black
gmt text cities_other.txt -F+f8p,Helvetica,black${HALO}+jBL -D0.1c/0.06c
printf "66.975 39.655\n" | gmt plot -Sc0.14c -Gwhite -W0.27p,black
printf "66.975 39.655 Samarkand\n" | gmt text -F+f8p,Helvetica,black${HALO}+jBR -D-0.1c/0.06c

gmt text -F+f+a+j -N <<EOF
76.2 28.0 11p,Helvetica-Bold,gray10${HALO} 0 CM Indian Plate
56.6 26.45 10p,Helvetica-Bold,gray10${HALO} 0 CM Arabian Plate
63.6 25.45 8p,Helvetica-Oblique,purple4${HALO} 0 CM Makran subduction zone
72.0 39.55 9p,Helvetica-Oblique,black${HALO} 0 CM Pamir
69.75 29.4 8p,Helvetica-Oblique,black${HALO} 70 CM Sulaiman Range
58.5 37.3 8p,Helvetica-Oblique,black${HALO} 0 CM Kopet Dag
EOF
gmt text -F+f+a+j -N <<'EOF'
EOF

gmt text -F+f+a+j -N <<EOF
61.8 38.5 8p,Helvetica-Oblique,92/51/23 0 CM Karakum Desert
63.4 40.5 8p,Helvetica-Oblique,92/51/23 0 CM Kyzylkum Desert
67.3 38.0 8p,Helvetica-Oblique,black${HALO} 0 CM Tajik Depression
71.6 40.75 8p,Helvetica-Oblique,black${HALO} 0 CM Fergana Valley
75.2 41.2 8p,Helvetica-Oblique,black${HALO} 0 CM Tien Shan
77.7 39.1 8p,Helvetica-Oblique,black${HALO} 0 CM Tarim Basin
77.0 35.6 8p,Helvetica-Oblique,black${HALO} -25 CM Karakoram
77.2 32.2 8p,Helvetica-Oblique,black${HALO} -38 CM Himalaya
73.0 26.2 8p,Helvetica-Oblique,92/51/23 0 CM Thar Desert
64.8 27.8 8p,Helvetica-Oblique,black${HALO} 0 CM Balochistan
58.0 30.5 8p,Helvetica-Oblique,92/51/23 0 CM Dasht-e Lut
71.0 31.6 8p,Helvetica-Oblique,dodgerblue4${HALO} 75 CM Indus R.
EOF
gmt text -F+f9p,Helvetica-Oblique,black${HALO}+j -N <<'EOF'
57.8 32.4 CM IRAN
72.6 29.7 CM PAKISTAN
61.8 39.3 CM TURKMENISTAN
64.8 41.4 CM UZBEKISTAN
68.6 39.95 CM TAJIKISTAN
76.3 29.6 CM INDIA
77.6 37.8 CM CHINA
EOF
gmt basemap -LjBR+c33+w300k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.22p,gray40+r0.1c
gmt basemap -TdjTR+w1.3c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold

gmt inset begin -DjTL+w2.9c+o0.15c/0.15c
  gmt coast -Rg -JG67.5/33.5/2.7c -Gwheat3 -Sazure2 -Dc -A5000 -EAF+gred3 -B0
  printf "$W $S\n$E $S\n$E $N\n$W $N\n$W $S\n" | gmt plot -Rg -JG67.5/33.5/2.7c -W0.44p,black -Am
gmt inset end

gmt colorbar -Ctopo.cpt -DJBC+w17c/0.32c+o0/1.0c+h+e -Bxa1000f250+l"Elevation / bathymetry (m)" -I0.5 \
  --FONT_ANNOT_PRIMARY=8p --FONT_LABEL=9p

cat > leg.txt <<'EOF'
N 3
S 0.25c - 0.5c - 1.2p,purple4 0.6c PB2002 plate boundary
S 0.25c - 0.5c - 0.49p,mediumpurple3,- 0.6c PB2002 orogen
S 0.25c c 0.12c magenta3 0.05p,black 0.6c Earthquake, depth > 70 km
S 0.25c - 0.5c - 0.76p,firebrick3 0.6c Thrust / reverse fault
S 0.25c - 0.5c - 0.76p,royalblue3 0.6c Strike-slip fault
S 0.25c - 0.5c - 0.76p,forestgreen 0.6c Normal fault
S 0.25c - 0.5c - 0.76p,gray20 0.6c Activity A: > 10 mm/yr
S 0.25c - 0.5c - 0.33p,gray20 0.6c Activity B: 1-10 mm/yr
S 0.25c - 0.5c - 0.14p,gray20 0.6c Activity C: indeterminate
S 0.25c c 0.13c white 0.27p,black 0.6c Regional city
S 0.25c v0.3c+e+a40+h0.3 0.45c black 0.65p,black 0.6c Plate velocity vs. Eurasia
EOF
gmt legend leg.txt -DJBC+w17c+o0/2.55c -F+gwhite+p0.22p,gray40 --FONT_ANNOT_PRIMARY=8p

gmt text -R0/17/0/1.5 -JX17c/1.5c -N -F+f8p,Helvetica,gray20+jTL -Y-6.3c <<EOF
0 1.45 Depicts: relief; PB2002 plate boundaries and Persia-Tibet-Burma orogen; Quaternary faults by kinematics and activity class;
0 1.05 earthquakes deeper than 70 km (M@~\263@~4, 1900-2025, n = $NDEEP); plate velocities relative to Eurasia from PB2002 Euler poles.
0 0.65 Data: SRTM15+ (GMT earth_relief, 30 arc-sec); Bird (2003) PB2002; USGS OFR 2007-1103; USGS ComCat; GSHHG; DCW.
0 0.25 Software: GMT $(gmt --version), GhostScript $(gs --version). Source: authors.
EOF
gmt end
cd ..
