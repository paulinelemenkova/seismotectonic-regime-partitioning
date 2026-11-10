#!/usr/bin/env python3
import json, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from scipy.stats import gaussian_kde

A, F = sys.argv[1], sys.argv[2]
os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 8, "axes.labelsize": 9,
                     "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.fontsize": 8,
                     "legend.frameon": False, "pdf.fonttype": 42})
DOMS = ["HKP", "NEC", "PBZ", "WEST"]
LAB = {"HKP": "Intermediate-depth domain", "NEC": "North-eastern crust",
       "PBZ": "Plate boundary", "WEST": "Western interior"}
COL = {"HKP": "#7b3294", "NEC": "#1b9e77", "PBZ": "#d95f02", "WEST": "#2c7fb8"}
MK = {"HKP": "o", "NEC": "s", "PBZ": "^", "WEST": "D"}
S = json.load(open(os.path.join(A, "stats.json")))
cat = pd.read_csv(os.path.join(A, "catalogue_final.csv"))
cat["time"] = pd.to_datetime(cat["time"], format="ISO8601", utc=True)
cat["yr"] = cat.time.dt.year + cat.time.dt.dayofyear / 365.25

def style(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", lw=0.3, color="0.85", zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def panel(ax, t):
    ax.text(-0.13, 1.03, t, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(F, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

mt = pd.read_csv(os.path.join(A, "mc_time.csv"))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
ax[0].fill_between(mt.year, mt.lo, mt.hi, color="#9ecae1", lw=0, label="Bootstrap 16-84 %")
ax[0].plot(mt.year, mt.mc, color="#08519c", lw=1.2, label="Maximum curvature + 0.2")
ax[0].axhline(5.5, color="#cb181d", lw=0.9, ls="--", label="Mc of Mw analysis (5.5)")
ax[0].set_xlabel("Year")
ax[0].set_ylabel("Magnitude of completeness (Mw)")
ax[0].set_xlim(1970, 2026)
ax[0].set_ylim(4.4, 5.8)
ax[0].legend(loc="lower left", ncol=1)
style(ax[0]); panel(ax[0], "(a)")
c = cat.sort_values("time")
for lab, sel, colr, ls in (("All events", np.ones(len(c), bool), "0.2", "-"),
                           ("Nearest-neighbour background", c.nnd_bg.values, "#08519c", "-"),
                           ("Gardner-Knopoff mainshocks", c.gk_main.values, "#cb181d", "--")):
    s = c[sel]
    ax[1].plot(s.yr, np.arange(1, len(s) + 1), color=colr, lw=1.1, ls=ls, label=f"{lab} (n = {len(s)})")
ax[1].set_xlabel("Year")
ax[1].set_ylabel("Cumulative number of events")
ax[1].set_xlim(1900, 2026)
ax[1].legend(loc="upper left")
style(ax[1]); panel(ax[1], "(b)")
fig.tight_layout()
save(fig, "fig_completeness")

fig, ax = plt.subplots(figsize=(4.6, 3.6))
for D in DOMS:
    m = pd.read_csv(os.path.join(A, f"fmd_{D}.csv")).m.values
    bins = np.round(np.arange(4.0, 8.01, 0.1), 1)
    cum = np.array([(m >= b - 1e-9).sum() for b in bins]) / (2026 - 1964)
    ok = cum > 0
    ax.semilogy(bins[ok], cum[ok], MK[D], ms=3.6, mfc=COL[D], mec="k", mew=0.3, ls="none")
    s = S["domains"][D]
    mm = np.linspace(5.5, 8.0, 50)
    ax.semilogy(mm, s["rate_mc"] * 10 ** (-s["b"] * (mm - 5.5)), color=COL[D], lw=1.2,
                label=f"{LAB[D]}: b = {s['b']:.2f} $\\pm$ {s['sb']:.2f}")
ax.axvline(5.5, color="0.4", lw=0.7, ls=":")
ax.text(5.53, 0.012, "Mc = 5.5", fontsize=8, color="0.3")
ax.set_xlabel("Moment magnitude, Mw")
ax.set_ylabel("Annual rate of events $\\geq$ Mw")
ax.set_xlim(4.4, 8.0)
ax.set_ylim(0.01, 200)
ax.legend(loc="upper right", fontsize=7.5)
style(ax)
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
save(fig, "fig_fmd_domains")

fig, ax = plt.subplots(figsize=(4.6, 3.4))
zz = np.linspace(0, 300, 601)
for D in DOMS:
    z = cat[(cat.domain == D) & ~cat.fixed_depth & cat.nnd_bg].depth.values
    k = gaussian_kde(z, bw_method=0.08)(zz)
    ax.plot(k / k.max(), zz, color=COL[D], lw=1.3, label=f"{LAB[D]} (n = {len(z)})")
    med = np.median(z)
    ax.plot(1.03, med, MK[D], color=COL[D], mec="k", mew=0.3, ms=4, clip_on=False)
for m_, s_, w_ in S["hkp_gmm"]:
    g = np.exp(-0.5 * ((zz - m_) / s_) ** 2)
    ax.plot(g * w_ / 0.53, zz, color=COL["HKP"], lw=0.8, ls="--")
ax.invert_yaxis()
ax.set_xlabel("Normalised kernel density")
ax.set_ylabel("Hypocentral depth (km)")
ax.set_xlim(0, 1.08)
ax.set_ylim(300, 0)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5)
style(ax)
save(fig, "fig_depth_histograms")

var = [("nnd", "Nearest-neighbour, Mc 5.5"), ("raw", "Undeclustered"), ("gk", "Gardner-Knopoff"),
       ("nnd_nofixed", "No default depths"), ("mc5.3", "Mc 5.3"), ("mc5.8", "Mc 5.8"),
       ("conv_gor", "Converted, regional GOR"), ("conv_scordilis", "Converted, Scordilis")]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.35, 1]})
for k, D in enumerate(DOMS):
    s = S["domains"][D]
    for i, (key, name) in enumerate(var):
        v = [s["mc"], s["b"], s["sb"], s["n"]] if key == "nnd" else s["sensitivity"].get(key)
        if v is None:
            continue
        y = i + (k - 1.5) * 0.18
        ax[0].errorbar(v[1], y, xerr=v[2], fmt=MK[D], color=COL[D], mec="k", mew=0.3, ms=3.8,
                       elinewidth=0.8, capsize=1.5, label=LAB[D] if i == 0 else None)
ax[0].set_yticks(range(len(var)))
ax[0].set_yticklabels([n for _, n in var])
ax[0].invert_yaxis()
ax[0].set_xlabel("b-value")
ax[0].legend(loc="lower center", bbox_to_anchor=(0.45, 1.02), ncol=2, fontsize=7)
style(ax[0]); ax[0].yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax[0].text(-0.62, 1.03, "(a)", transform=ax[0].transAxes, fontsize=11, fontweight="bold", va="bottom")
x = np.arange(len(DOMS))
for j, (key, name, hatch) in enumerate((("ret_m7_yr", "GR model, 1964-2025", None),)):
    vals = [S["domains"][D][key] for D in DOMS]
    ci = np.array([S["domains"][D]["ret_m7_ci"] for D in DOMS]).T
    ax[1].bar(x - 0.18, vals, 0.34, color=[COL[D] for D in DOMS], edgecolor="k", lw=0.4,
              label="Gutenberg-Richter model")
    ax[1].errorbar(x - 0.18, vals, yerr=[np.array(vals) - ci[0], ci[1] - np.array(vals)], fmt="none",
                   ecolor="k", elinewidth=0.7, capsize=2)
obs = [126 / max(S["domains"][D]["n_obs_m7"], 1) for D in DOMS]
ax[1].bar(x + 0.18, obs, 0.34, color="white", edgecolor=[COL[D] for D in DOMS], lw=1.0, hatch="////",
          label="Observed, 1900-2025")
ax[1].set_yscale("log")
ax[1].set_xticks(x)
ax[1].set_xticklabels(["HKP", "NEC", "PBZ", "WEST"])
ax[1].set_ylabel("Return period of Mw $\\geq$ 7 (yr)")
ax[1].legend(loc="upper left", fontsize=7)
style(ax[1]); ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); panel(ax[1], "(b)")
fig.tight_layout()
save(fig, "fig_sensitivity")
print("figures done")
