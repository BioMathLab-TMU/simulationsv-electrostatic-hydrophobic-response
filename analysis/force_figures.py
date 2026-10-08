"""Prespecified, print-safe figures for the force-response experiment."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plotting import COLORS, save_figure

COLORS_BY_ARM = {(0, 0): COLORS["neutral"], (0, 1): COLORS["syn_alone"],
                 (1, 0): COLORS["ves_alone"], (1, 1): COLORS["both_inside"]}


def make_figures(run_level: pd.DataFrame, contrasts: pd.DataFrame,
                 shells: pd.DataFrame, output_dir) -> None:
    metrics = (
        ("late_mean_largest_component_fraction", "Largest component / all particles"),
        ("late_mean_synapsin_recruited_fraction", "Recruited synapsin fraction"),
        ("late_mean_largest_component_vesicle_fraction", "Vesicles in largest component / all vesicles"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4), constrained_layout=True)
    for ax, (metric, label) in zip(axes, metrics):
        for x, arm in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            group = run_level[(run_level.E == arm[0]) & (run_level.H == arm[1])]
            values = group[metric].to_numpy(dtype=float)
            ax.scatter(np.full(len(values), x) + np.linspace(-0.08, 0.08, len(values)), values,
                       color=COLORS_BY_ARM[arm], s=22, zorder=3)
            ax.plot([x - 0.22, x + 0.22], [values.mean()] * 2,
                    color=COLORS_BY_ARM[arm], lw=2.0)
        ax.set_xticks(range(4), ["E0H0", "E0H1", "E1H0", "E1H1"])
        ax.set(ylim=(0, 1), ylabel=label)
        ax.grid(axis="y", alpha=0.5)
    save_figure(fig, output_dir, "registered_endpoints")

    primary = contrasts[contrasts.metric == "late_mean_largest_component_fraction"]
    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    labels = primary["contrast"].tolist()
    y = np.arange(len(primary))
    x = primary["mean"].to_numpy(dtype=float)
    errors = np.vstack((x - primary["ci95_low"].to_numpy(dtype=float),
                        primary["ci95_high"].to_numpy(dtype=float) - x))
    ax.errorbar(x, y, xerr=errors, fmt="o", color=COLORS["both_inside"], capsize=3)
    ax.axvline(0, color=COLORS["neutral"], lw=0.8, ls="--")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Paired difference in LCC fraction (95% t CI)")
    fig.subplots_adjust(left=0.28, bottom=0.23, right=0.97, top=0.96)
    fig.text(0.5, 0.04, "n=3 paired seeds; exact sign-flip p≥0.25",
             ha="center", fontsize=7)
    save_figure(fig, output_dir, "planned_paired_contrasts")

    late = shells[shells.frame_index >= 140]  # zero-based final 60 / 200 frames
    summaries = late.groupby(["seed", "E", "H", "analysis_shell_nm"], as_index=False)["largest_component_fraction"].mean()
    fig, ax = plt.subplots(figsize=(5.5, 3.5), constrained_layout=True)
    for arm, group in summaries.groupby(["E", "H"]):
        stat = group.groupby("analysis_shell_nm")["largest_component_fraction"].agg(["mean", "std"])
        ax.plot(stat.index, stat["mean"], marker="o", color=COLORS_BY_ARM[arm],
                label=f"E{arm[0]}H{arm[1]}")
        ax.fill_between(stat.index.to_numpy(dtype=float),
                        (stat["mean"] - stat["std"]).to_numpy(dtype=float),
                        (stat["mean"] + stat["std"]).to_numpy(dtype=float),
                        color=COLORS_BY_ARM[arm], alpha=0.12)
    ax.set(xlabel="Surface-gap analysis shell (nm)",
           ylabel="Late-window largest-component fraction", ylim=(0, 1))
    ax.set_xticks([0, 2, 4, 6])
    ax.legend()
    save_figure(fig, output_dir, "analysis_shell_sensitivity")
