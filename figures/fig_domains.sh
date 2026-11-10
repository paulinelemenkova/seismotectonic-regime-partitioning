#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
CAT=${CAT:-catalogue_final.csv}
W=58; E=76; S=27; N=41; R=$W/$E/$S/$N; J=M17c
mkdir -p work && cd work
case "$CAT" in http*) curl -s -m 120 -o cat.csv "$CAT";; *) cp "$CAT" cat.csv;; esac
python3 - <<'PY'
import pandas as pd
c = pd.read_csv("cat.csv")
c = c.sort_values("mw")
col = {"HKP": 1, "NEC": 2, "PBZ": 3, "WEST": 4}
c["k"] = c.domain.map(col)
c["s"] = 0.03 * 1.9 ** (c.mw - 3)
c[["longitude", "latitude", "k", "s"]].to_csv("eq.txt", sep=" ", header=False, index=False)
bg = c[c.hdb_label >= 0]
med = bg.groupby("hdb_label").agg(lon=("longitude", "median"), lat=("latitude", "median"),
                                  dep=("depth", "median"), n=("mw", "size")).reset_index()
med.to_csv("clusters.txt", sep=" ", header=False, index=False,
           columns=["lon", "lat", "hdb_label"])
k = c[c.time.str[:10].isin(["2022-06-21", "2023-10-07", "2023-10-11", "2023-10-15", "2025-08-31", "2025-11-02"]) & (c.mw >= 5.9)]
k[["longitude", "latitude"]].to_csv("key.txt", sep=" ", header=False, index=False)
print(len(c), med.shape)
PY
NEQ=$(wc -l < eq.txt)
curl -s -m 120 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
ogr2ogr -f OGR_GMT qf.gmt faults.shp -t_srs EPSG:4326
gmt grdcut @earth_relief_30s -R$R -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.8 -Gint.nc
gmt coast -EAF -M -Di > af.txt
gmt makecpt -Cgray -T-16000/6000/250 > relief.cpt
cat > dom.cpt <<'EOF'
0.5 123/50/148 1.5 123/50/148
1.5 27/158/119 2.5 27/158/119
2.5 217/95/2 3.5 217/95/2
3.5 44/127/184 4.5 44/127/184
EOF

gmt begin $OUT/fig_domains pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        MAP_GRID_PEN_PRIMARY 0.25p,gray55 FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel
gmt grdimage relief.nc -Iint.nc -Crelief.cpt -R$R -J$J
gmt coast -Di -W0.3p,gray30 -N1/0.6p,gray25,- -I1/0.35p,steelblue -C+l+gsteelblue -A50 -Bxa2f1g2 -Bya2f1g2 -BWSen
gmt plot qf.gmt -W0.25p,gray35
gmt plot af.txt -W1.4p,black
gmt plot eq.txt -Sc -Cdom.cpt -W0.1p,black@50 -t25
gmt plot key.txt -Sa0.55c -W1.1p,black -Gwhite@100
gmt plot clusters.txt -Sc0.42c -Gwhite -W0.8p,black
gmt text clusters.txt -F+f6.5p,Helvetica-Bold,black+jCM
gmt text -F+f+j -Gwhite@25 -C0.06c/0.05c <<'EOF'
72.9 40.55 10p,Helvetica-Bold,27/158/119 CM North-eastern crust
73.75 36.68 10p,Helvetica-Bold,123/50/148 CM intermediate-depth zone
68.6 28.15 10p,Helvetica-Bold,217/95/2 CM Plate-boundary domain
60.3 30.2 10p,Helvetica-Bold,44/127/184 CM Western interior
EOF
gmt basemap -LjBR+c34+w200k+f+l"km"+o0.4c/0.5c -F+gwhite@20+p0.4p,gray40+r0.1c
gmt basemap -TdjTR+w1.2c+l,,,N+o0.35c/0.35c --FONT_TITLE=8p,Helvetica-Bold
cat > leg.txt <<'EOF'
N 2
S 0.25c c 0.25c 123/50/148 0.2p,black 0.55c Intermediate-depth domain
S 0.25c c 0.25c 27/158/119 0.2p,black 0.55c North-eastern crust
S 0.25c c 0.25c 217/95/2 0.2p,black 0.55c Plate boundary
S 0.25c c 0.25c 44/127/184 0.2p,black 0.55c Western interior
S 0.25c c 0.35c white 0.8p,black 0.55c HDBSCAN cluster (median, label)
S 0.25c a 0.4c white 1.0p,black 0.55c Destructive earthquake, 2022-2025
EOF
gmt legend leg.txt -DJBC+w17c+o0/0.9c -F+gwhite+p0.4p,gray40
gmt text -R0/17/0/1.6 -JX17c/1.6c -N -F+f8p,Helvetica,gray20+jTL -Y-4.0c <<EOF
0 1.55 Depicts: earthquakes Mw@~\263@~4, 1900-2025 (n = $NEQ), coloured by seismotectonic domain and sized by magnitude; medians of the
0 1.15 14 HDBSCAN clusters of declustered hypocentres; Quaternary faults; destructive events of 2022-2025 (stars).
0 0.75 Data: USGS ComCat; ISC-GEM; USGS OFR 2007-1103; SRTM15+ (GMT earth_relief, 30 arc-sec); GSHHG; DCW.
0 0.35 Software: GMT $(gmt --version), GhostScript $(gs --version), Python scikit-learn. Source: authors.
EOF
gmt end
