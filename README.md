# Seismotectonic Regime Partitioning

Reproducible pipeline for the statistical partitioning of a composite hypocentral catalogue into seismotectonic regimes, with frequency-magnitude, moment-tensor and fault-structure diagnostics.

## Methods

- Catalogue fusion with space-time duplicate matching and magnitude homogenisation (general orthogonal regression, Scordilis relations)
- Declustering by nearest-neighbour proximity (Zaliapin-Ben-Zion), thresholded with a two-component Gaussian mixture; Gardner-Knopoff windows as a control
- Density-based regime delineation with HDBSCAN on depth-weighted hypocentral features and k-nearest-neighbour propagation
- Completeness by maximum curvature; Aki-Utsu maximum-likelihood b-values with Shi-Bolt errors, bootstrap intervals and gridded b-value mapping
- Poisson recurrence and exceedance probabilities from Gutenberg-Richter parameters
- Frohlich ternary classification of centroid moment tensors; orientation-tensor P-axis statistics
- Kostrov summation of moment tensors for principal strain-rate fields
- Seismic moment budgets and their temporal release
- Fault-density gridding, length-weighted strike roses and KD-tree epicentre-fault distances tested against a random null (two-sample Kolmogorov-Smirnov)

## Layout

```
analysis/   catalogue processing and statistics (Python)
figures/    map and diagram generation (GMT 6, Bash)
```

## Requirements

Python 3.11+ (`pip install -r requirements.txt`), GMT 6.5+, GDAL/OGR, curl.

## Usage

```
python analysis/analysis.py <comcat.csv> <iscgem.csv> <outdir>
python analysis/bgrid_mb.py <outdir>/catalogue_final.csv <outdir>
python analysis/meca_stats.py <outdir>/catalogue_final.csv gcmt_ndk_subset.txt <outdir>/meca_stats.json
python analysis/moment_release.py <outdir>/catalogue_final.csv <outdir>
python analysis/synthesis_inputs.py <outdir>/catalogue_final.csv gcmt_ndk_subset.txt kostrov_cells_shallow.txt <outdir>
OUT=<figdir> bash figures/<script>.sh
```

All inputs are open data retrieved at run time (USGS ComCat, ISC-GEM, Global CMT, SRTM15+, GSHHG, PB2002, USGS fault databases).
