"""Render the four article figures from released numerical inputs.

Arial 10 pt, maximum width 7.5 in, embedded vector PDF and 600-dpi LZW TIFF.
Existing figures are replaced atomically, preserving historical hard links.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARMS = [(0, 0), (0, 1), (1, 0), (1, 1)]
COLORS = ["#737373", "#E76F51", "#2A9D8F", "#195E56"]
MARKERS = ["o", "s", "^", "D"]
LABELS = [f"E{e}H{h}" for e, h in ARMS]
STYLES = [":", "--", "-.", "-"]


def finish(fig, out, name):
    for ext in ("pdf", "png", "tif"):
        target = out/f"{name}.{ext}"
        fd, temp = tempfile.mkstemp(prefix=".figure_", suffix="."+ext, dir=out)
        os.close(fd)
        kw = {"dpi": 600 if ext == "tif" else 180}
        if ext == "tif":
            kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        try:
            fig.savefig(temp, **kw)
            os.replace(temp, target)
        finally:
            if os.path.exists(temp):
                os.unlink(temp)
    plt.close(fig)


def axes_style(ax, panel):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.15)
    ax.set_title(panel, loc="left", fontweight="bold", fontsize=12)


def dots(ax, table, column):
    for x, ((e, h), color, marker) in enumerate(zip(ARMS, COLORS, MARKERS)):
        v = table.loc[(table.E == e) & (table.H == h), column].to_numpy(float)
        if len(v) != 3:
            raise ValueError("expected three independent seed blocks")
        ax.scatter(x+np.linspace(-.08, .08, 3), v, s=30, facecolor="white",
                   edgecolor=color, marker=marker, zorder=4)
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color=color, marker="_",
                    markersize=17, capsize=3, lw=1.3, zorder=3)
    ax.set_xticks(range(4), LABELS, rotation=35, ha="right")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--font", help="Explicit Matplotlib font family; otherwise select an available font")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    candidates = [a.font] if a.font else ["Arial", "Liberation Sans", "DejaVu Sans"]
    family = None
    for candidate in candidates:
        try:
            font = font_manager.findfont(candidate, fallback_to_default=False)
            family = candidate
            break
        except ValueError:
            continue
    if family is None:
        raise ValueError("Requested figure font is unavailable")
    plt.rcParams.update({"font.family": family, "font.size": 10,
                         "axes.labelsize": 10, "xtick.labelsize": 9,
                         "ytick.labelsize": 9, "legend.fontsize": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    endpoints = pd.read_csv(ROOT/"data/processed/run_endpoints.csv")
    frames = pd.read_csv(ROOT/"data/processed/frame_metrics.csv")
    interface = pd.read_csv(ROOT/"results/interface_diagnostics/interface_run_metrics.csv")
    mobility = pd.read_csv(ROOT/"results/mobility_diagnostic/run_level_mobility.csv")
    metrics = ["late_mean_largest_component_fraction", "late_mean_synapsin_recruited_fraction",
               "late_mean_largest_component_vesicle_fraction"]
    labels = ["Largest component /\nall agents", "Non-isolated Syn /\nall Syn", "SVs in largest component /\nall SVs"]
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 3.25), layout="constrained")
    for k, (ax, metric, label) in enumerate(zip(axs, metrics, labels)):
        dots(ax, endpoints, metric)
        ax.set(ylabel=label, ylim=(0, 1))
        axes_style(ax, chr(65+k))
    finish(fig, a.output, "force_factorial_endpoints")

    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.7), layout="constrained")
    series = {}
    for (e, h), color, marker, label, style in zip(ARMS, COLORS, MARKERS, LABELS, STYLES):
        g = frames[(frames.E == e) & (frames.H == h)]
        s = g.groupby("frame_index").agg(time=("time_ms", "mean"),
            mean=("largest_component_fraction", "mean"), sd=("largest_component_fraction", "std"))
        series[e, h] = s
        x = s.time.to_numpy()
        axs[0].plot(x, s["mean"], color=color, ls=style, lw=1.5, label=label)
        axs[0].fill_between(x, s["mean"]-s.sd, s["mean"]+s.sd, color=color, alpha=.12)
        values = endpoints[(endpoints.E == e) & (endpoints.H == h)]
        xx, yy = values[metrics[1]], values[metrics[2]]
        axs[1].scatter(xx, yy, marker=marker, edgecolor=color, facecolor="white", s=31, zorder=4)
        axs[1].errorbar(xx.mean(), yy.mean(), xerr=xx.std(ddof=1), yerr=yy.std(ddof=1),
                        marker=marker, color=color, markersize=5, capsize=3)
    ref = series[1, 0]["mean"]+series[0, 1]["mean"]-series[0, 0]["mean"]
    axs[0].plot(series[0, 0].time, ref, color="black", ls=(0, (4, 2)), lw=1, label="Additive reference")
    late_start = frames.loc[frames.frame_index == 140, "time_ms"].mean()
    axs[0].axvspan(late_start, 28, color="gray", alpha=.08, zorder=0)
    axs[0].set(xlabel="Time (ms)", ylabel="Largest component / all agents", ylim=(0, .8), xlim=(0, 28))
    axs[0].legend(loc="upper left", frameon=False)
    axs[1].set(xlabel="Non-isolated Syn / all Syn", ylabel="SVs in largest component / all SVs",
               xlim=(.25, .95), ylim=(0, .8))
    for k, ax in enumerate(axs):
        axes_style(ax, chr(65+k))
    finish(fig, a.output, "force_factorial_timecourse")

    fig, axs = plt.subplots(1, 3, figsize=(7.5, 3.25), layout="constrained")
    cols = ["syn_near_sv_fraction", "syn_multivesicle_neighbor_fraction", "edges_SV"]
    labels = ["Syn with at least one SV\nneighbor / all Syn", "Syn with at least two SV\nneighbors / all Syn", "Syn-SV edges\nper snapshot"]
    for k, (ax, col, label) in enumerate(zip(axs, cols, labels)):
        dots(ax, interface, col)
        ax.set(ylabel=label, ylim=(0, 1) if k < 2 else (0, 2700))
        axes_style(ax, chr(65+k))
    finish(fig, a.output, "force_interface_diagnostics")

    full = mobility[mobility.window == "full_28ms"].copy()
    full["E"] = full.arm.str[1].astype(int)
    full["H"] = full.arm.str[3].astype(int)
    fig, ax = plt.subplots(figsize=(6.6, 3.6), layout="constrained")
    dots(ax, full, "D_mean_um2_s")
    reference = json.loads((ROOT/"data/external/hoffmann_mobility_reference.json").read_text())
    ax.axhline(reference["with_SV"]["mean"], color="black", ls="--", lw=1.1,
               label="Experiment with SVs (mean)")
    nominal = 1.380649e-23*310.15/(6*np.pi*.35*3.972e-9)*1e12
    ax.axhline(nominal, color="gray", ls=":", lw=1.1, label="Nominal free Syn input")
    ax.set(ylabel="Apparent short-lag mobility ($\\mu$m$^2$/s)", ylim=(0, .185))
    ax.legend(loc="upper right", frameon=False)
    axes_style(ax, "")
    finish(fig, a.output, "force_mobility_external_comparison")
    (a.output/"figure_provenance.json").write_text(json.dumps({
        "font": family, "font_resolved": Path(font).name, "font_size_pt": 10,
        "width_inches_max": 7.5, "tiff_dpi": 600, "compression": "LZW",
        "palette": dict(zip(LABELS, COLORS)), "independent_replicates_per_arm": 3,
        "errors": "between-seed SD", "nominal_Syn_D_um2_s": nominal}, indent=2))
    print("Four figure sets rendered: PDF, PNG and TIFF")


if __name__ == "__main__":
    main()
