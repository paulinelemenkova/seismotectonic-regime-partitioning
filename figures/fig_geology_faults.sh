#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
W=60; E=75.5; S=29; N=38.8; R=$W/$E/$S/$N; J=M17c
GDRIVE_ID=${GDRIVE_ID:-1TLTNVIpkC4p7Mff4VkBo8mUVRfA3Yh6v}
mkdir -p work && cd work

curl -s -f -m 600 -L -o afg.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o afg.zip -x '__MACOSX/*'
for f in glgafg fltafg; do
  P=$(find . -name "$f.shp" | head -1); for x in shp shx dbf prj; do cp "${P%.shp}.$x" $f.$x; done
done
P=$(find . -name 'geoBoundaries-AFG-ADM1_simplified.shp' | head -1)
for x in shp shx dbf prj; do cp "${P%.shp}.$x" adm1.$x; done
curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)

SQL="SELECT GEOMETRY, CASE
 WHEN GEN_AGE LIKE 'Q%' THEN 1
 WHEN GEN_AGE LIKE 'N%' THEN 2
 WHEN GEN_AGE LIKE 'Pg%' THEN 3
 WHEN GEN_AGE IN ('K','K1','K2','KJ') THEN 4
 WHEN GEN_AGE IN ('J','J1-2','J2-3','J3','JTr') THEN 5
 WHEN GEN_AGE IN ('Tr','TrP','TK') THEN 6
 WHEN GEN_AGE IN ('P','PC','C','CD') THEN 7
 WHEN GEN_AGE IN ('D','DS','S','O','Cm') THEN 8
 WHEN GEN_AGE IN ('Pt','A','CmPt') THEN 9
 ELSE 0 END AS CL FROM glgafg"
ogr2ogr -f OGR_GMT geo.gmt glgafg.shp -t_srs EPSG:4326 -dialect SQLite -sql "$SQL"
ogr2ogr -f OGR_GMT flt.gmt fltafg.shp -t_srs EPSG:4326 -where "TYPE NOT LIKE '%trust%'"
ogr2ogr -f OGR_GMT thr.gmt fltafg.shp -t_srs EPSG:4326 -where "TYPE LIKE '%trust%'"
ogr2ogr -f OGR_GMT adm1.gmt adm1.shp -t_srs EPSG:4326
ogr2ogr -f OGR_GMT qR.gmt  "$QSHP" -t_srs EPSG:4326 -where "Type='R'"
ogr2ogr -f OGR_GMT qSS.gmt "$QSHP" -t_srs EPSG:4326 -where "Type IN ('RL','LL','SS')"
ogr2ogr -f OGR_GMT qN.gmt  "$QSHP" -t_srs EPSG:4326 -where "Type='N'"
fc() { ogrinfo -so -al "$@" | awk -F': ' '/Feature Count/{print $2}'; }
NPOLY=$(fc glgafg.shp); NFLT=$(fc fltafg.shp -where "TYPE NOT LIKE '%trust%'"); NTHR=$(fc fltafg.shp -where "TYPE LIKE '%trust%'")
NQR=$(fc "$QSHP" -where "Type='R'"); NQS=$(fc "$QSHP" -where "Type IN ('RL','LL','SS')"); NQN=$(fc "$QSHP" -where "Type='N'")
NQ=$((NQR + NQS + NQN))
echo "polygons $NPOLY faults $NFLT thrusts $NTHR quaternary $NQ"

cat > geo.cpt <<'EOF'
-0.5 200/200/200 0.5 200/200/200 ;Undifferentiated
0.5 249/249/127 1.5 249/249/127 ;Quaternary
1.5 255/230/25 2.5 255/230/25 ;Neogene
2.5 253/154/82 3.5 253/154/82 ;Paleogene
3.5 127/198/78 4.5 127/198/78 ;Cretaceous
4.5 52/178/201 5.5 52/178/201 ;Jurassic
5.5 129/43/146 6.5 129/43/146 ;Triassic
6.5 103/165/153 7.5 103/165/153 ;Carboniferous-Permian
7.5 203/140/55 8.5 203/140/55 ;Cambrian-Devonian
8.5 247/53/99 9.5 247/53/99 ;Precambrian
EOF

gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.9 -Gint.nc
gmt makecpt -Cgray -T-12000/7000/250 > relief.cpt
gmt coast -EAF -M -Di > af.txt

cat > cities.txt <<'EOF'
EOF
cat > blocks.txt <<'EOF'
68.25 32.25 30 Katawaz Basin
73.8 37.95 0 Pamir
EOF
cat > faults_lab.txt <<'EOF'
EOF

gmt begin $OUT/fig_geology_faults pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        MAP_GRID_PEN_PRIMARY 0.25p,gray45 FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y11.5c
gmt plot geo.gmt -Cgeo.cpt -aZ=CL -G+z -t30
gmt coast -Di -W0.3p,gray25 -N1/0.7p,gray15,- -I1/0.4p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot adm1.gmt -W0.35p,gray20,4_2:0p
gmt plot flt.gmt -W0.3p,gray15
gmt plot thr.gmt -W1.0p,gray35
gmt plot qN.gmt  -W1.5p,white
gmt plot qN.gmt  -W0.75p,0/70/255
gmt plot qSS.gmt -W1.5p,white
gmt plot qSS.gmt -W0.75p,black
gmt plot qR.gmt  -W1.5p,white
gmt plot qR.gmt  -W0.75p,220/20/60
gmt plot af.txt -W1.6p,black
gmt text blocks.txt -F+f9p,Helvetica-BoldOblique,black+a+jCM -Gwhite@30 -C0.06c/0.04c
gmt text faults_lab.txt -F+f8p,Helvetica-Oblique,220/20/60+a+jCM -Gwhite@20 -C0.05c/0.03c
gmt plot cities.txt -Ss0.2c -Gwhite -W0.8p,black
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -TdjTL+w1.2c+l,,,N+o0.4c/0.4c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<EOF
H 9p,Helvetica-Bold Chronostratigraphic units (USGS 1:850 000 geology)
N 5
S 0.25c s 0.32c 249/249/127 0.3p,gray30 0.55c Quaternary
S 0.25c s 0.32c 255/230/25 0.3p,gray30 0.55c Neogene
S 0.25c s 0.32c 253/154/82 0.3p,gray30 0.55c Paleogene
S 0.25c s 0.32c 127/198/78 0.3p,gray30 0.55c Cretaceous
S 0.25c s 0.32c 52/178/201 0.3p,gray30 0.55c Jurassic
S 0.25c s 0.32c 129/43/146 0.3p,gray30 0.55c Triassic
S 0.25c s 0.32c 103/165/153 0.3p,gray30 0.55c Carboniferous-Permian
S 0.25c s 0.32c 203/140/55 0.3p,gray30 0.55c Cambrian-Devonian
S 0.25c s 0.32c 247/53/99 0.3p,gray30 0.55c Precambrian
S 0.25c s 0.32c 200/200/200 0.3p,gray30 0.55c Undifferentiated
G 0.1c
H 9p,Helvetica-Bold Faults
N 3
S 0.3c - 0.5c - 1.2p,220/20/60 0.7c Quaternary reverse ($NQR)
S 0.3c - 0.5c - 1.2p,black 0.7c Quaternary strike-slip ($NQS)
S 0.3c - 0.5c - 1.1p,0/70/255 0.7c Quaternary normal ($NQN)
S 0.3c - 0.5c - 1.0p,gray35 0.7c Mapped thrust ($NTHR)
S 0.3c - 0.5c - 0.3p,gray15 0.7c Mapped fault ($NFLT)
S 0.3c - 0.5c - 0.35p,gray20,4_2:0p 0.7c Provincial boundary
EOF
gmt legend leg.txt -DJBC+w17c+o0/0.8c -F+gwhite+p0.4p,gray40
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -Y-5.3c <<EOF
0 1.55 Depicts: $NPOLY geological polygons in 10 chronostratigraphic classes, $NFLT mapped faults and $NTHR thrusts, $NQ Quaternary faults (USGS);
0 1.15 provincial boundaries, principal cities, tectonic blocks and major faults on SRTM15+ shaded relief.
0 0.35 Software: GMT $(gmt --version), GDAL $(gdal-config --version 2>/dev/null || ogr2ogr --version | awk '{print $2}' | tr -d ,), GhostScript $(gs --version). Source: authors.
EOF
gmt end
