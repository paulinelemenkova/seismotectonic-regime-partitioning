#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
IN=${IN:-.}
R=58/77/27/41.5; J=M17c
mkdir -p work && cp $IN/dom.xyz $IN/hkp.xyz $IN/meca.txt $IN/cross.txt work/ && cd work

gmt xyz2grd dom.xyz -R58/77/27/41.5 -I0.5 -r -Gdom.nc
gmt xyz2grd hkp.xyz -R58/77/27/41.5 -I0.5 -r -Ghkp.nc
gmt grdmath hkp.nc 0 DENAN = hkp0.nc
gmt grdsample hkp0.nc -I0.05 -Ghkp_f.nc -nb
cat > dom.cpt <<'EOF'
0.5 0/175/200 1.5 0/175/200 ;NEC
1.5 220/20/60 2.5 220/20/60 ;PBZ
2.5 240/160/0 3.5 240/160/0 ;WEST
EOF
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt makecpt -Cgray -T-14000/6500/250 > relief.cpt
gmt coast -EAF -M -Di > af.txt
cat > recent.txt <<'EOF'
EOF
C_HKP=0/70/255; C_NEC=0/175/200; C_PBZ=220/20/60; C_WEST=240/160/0

gmt begin $OUT/fig_synthesis pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        MAP_GRID_PEN_PRIMARY 0.25p,gray40 FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J -Y11.0c
gmt grdimage dom.nc -Cdom.cpt -Q -t45 -nn
gmt coast -Di -W0.3p,gray25 -N1/0.6p,gray20,- -I1/0.35p,steelblue -A50 -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt grdcontour hkp_f.nc -C0.5 -W3.2p,white
gmt grdcontour hkp_f.nc -C0.5 -W1.8p,$C_HKP
gmt plot af.txt -W1.5p,black
gmt velo cross.txt -Sx0.42c -W2.2p,white -Gwhite -A0.14c+e+a50
gmt velo cross.txt -Sx0.42c -W0.9p,black -Gblack -A0.12c+e+a50
awk 'NR==1' meca.txt | gmt meca -Sm1.25c+m -G$C_HKP -Ewhite -W0.8p,black
awk 'NR==2' meca.txt | gmt meca -Sm1.25c+m -G$C_NEC -Ewhite -W0.8p,black
awk 'NR==3' meca.txt | gmt meca -Sm1.25c+m -G$C_PBZ -Ewhite -W0.8p,black
awk 'NR==4' meca.txt | gmt meca -Sm1.25c+m -G$C_WEST -Ewhite -W0.8p,black
awk '{print $1, $2}' meca.txt | paste - <(printf "HKP\nNEC\nPBZ\nWEST\n") | \
  gmt text -F+f9p,Helvetica-Bold,black+jTC -D0/-0.75c -Gwhite@10 -C0.06c/0.04c
gmt plot recent.txt -Sa0.45c -Gwhite -W1p,black
awk '{print $1, $2, $3" "$4}' recent.txt | gmt text -F+f7.5p,Helvetica-Bold+jLM -D0.3c/0 -Gwhite@15 -C0.04c/0.03c
printf "58.2 41.35 (a) Seismicity regimes\n" | gmt text -F+f9.5p,Helvetica-Bold+jTL -Gwhite@10 -C0.06c/0.05c
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.1c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<EOF
N 3
S 0.25c s 0.32c $C_NEC@45 0.3p,gray30 0.55c NEC crustal domain
S 0.25c s 0.32c $C_PBZ@45 0.3p,gray30 0.55c PBZ crustal domain
S 0.25c s 0.32c $C_WEST@45 0.3p,gray30 0.55c WEST crustal domain
S 0.25c - 0.45c - 1.8p,$C_HKP 0.55c HKP footprint (70-300 km)
S 0.25c c 0.32c white 0.6p,black 0.55c Summed CMT mechanism
S 0.25c a 0.32c white 0.8p,black 0.55c 2022-2025 destructive event
G 0.05c
L 7.5p,Helvetica L Crosses: principal horizontal strain axes from Kostrov summation (inward arrows = shortening)
EOF
gmt legend leg.txt -DJBC+w17c+o0/0.75c -F+gwhite+p0.4p,gray40
gmt basemap -R0/17/0/4.1 -JX17c/4.1c -Y-6.8c -B+n
cat > boxes.txt <<EOF
0.0|$C_HKP|HKP|Intermediate-depth slab|b = 0.97 @~\261@~ 0.07|Thrust 72 %, P N165@.E|Moment share 45 %|Median depth 152 km|T(M@-w@- @~\263@~ 7) = 11 yr
4.3|$C_NEC|NEC|North-eastern crust|b = 1.33 @~\261@~ 0.16|Thrust 43 %, P N151@.E|Moment share 21 %|Median depth 30 km|T(M@-w@- @~\263@~ 7) = 70 yr
8.6|$C_PBZ|PBZ|Plate-boundary zone|b = 1.38 @~\261@~ 0.19|Thrust 54 %, P N173@.E|Moment share 19 %|Median depth 25 km|T(M@-w@- @~\263@~ 7) = 87 yr
12.9|$C_WEST|WEST|Western interior|b = 0.69 @~\261@~ 0.11|Strike-slip 32 %, P N041@.E|Moment share 15 %|Median depth 15 km|T(M@-w@- @~\263@~ 7) = 19 yr
EOF
while IFS='|' read -r x0 col t1 t2 t3 t4 t5 t6 t7; do
  x1=$(awk -v a=$x0 'BEGIN{print a+4.1}')
  xm=$(awk -v a=$x0 'BEGIN{print a+0.15}')
  printf "%s 0\n%s 0\n%s 3.6\n%s 3.6\n" $x0 $x1 $x1 $x0 | gmt plot -Gwhite -W0.8p,$col -L
  printf "%s 3.6\n%s 3.6\n%s 4.1\n%s 4.1\n" $x0 $x1 $x1 $x0 | gmt plot -G$col -W0.8p,$col -L
  echo "$xm 3.85 $t1" | gmt text -F+f9p,Helvetica-Bold,white+jLM -N
  y=3.25
  for t in "$t2" "$t3" "$t4" "$t5" "$t6" "$t7"; do
    echo "$xm $y $t" | gmt text -F+f7.5p,Helvetica,black+jLM -N
    y=$(awk -v a=$y 'BEGIN{print a-0.56}')
  done
done < boxes.txt
printf "0 4.4 (b) Domain statistics (this study)\n" | gmt text -F+f9.5p,Helvetica-Bold+jBL -N
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -Y-2.3c <<EOF
0 1.55 Depicts: (a) majority crustal domain per 0.5-degree cell, intermediate-depth footprint, summed Global CMT moment tensors (1976-2025),
0 1.15 Kostrov strain axes and 2022-2025 destructive events; (b) b-values (M@-w@- >= 5.5, 1964-2025), mechanisms, moment shares, depths, return periods.
0 0.75 Data: this study (merged ComCat and ISC-GEM catalogue; Global CMT); SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python. Source: authors.
EOF
gmt end
