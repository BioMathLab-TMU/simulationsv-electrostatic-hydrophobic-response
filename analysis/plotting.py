"""Shared publication palette and vector/raster figure export."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, Mapping

import matplotlib as mpl

# Batch runs must not depend on an available desktop display. If pyplot was
# already imported by Jupyter, its inline backend is preserved.
if "matplotlib.pyplot" not in sys.modules:
    mpl.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


COLORS: Dict[str, str] = {
    "ves_alone": "#2A6EBB",
    "ves_inside": "#1B4F8A",
    "ves_outside": "#7FB3D9",
    "syn_alone": "#D4541A",
    "syn_inside": "#9B2C0A",
    "syn_outside": "#F0A070",
    "both_inside": "#5B2D8E",
    "both_outside": "#B08DC8",
    "neutral": "#6B7280",
    "highlight": "#F5C518",
    "background": "#F8F9FA",
    "text_dark": "#1A1A2E",
    "grid": "#D9DEE7",
}

CONDITION_COLORS = [
    COLORS["both_inside"], COLORS["ves_alone"], COLORS["syn_alone"],
    COLORS["ves_inside"], COLORS["syn_inside"], COLORS["both_outside"],
    COLORS["ves_outside"], COLORS["syn_outside"],
]
LINE_STYLES = ["-", "--", "-.", ":"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]

# Factor-specific mappings are fixed before looking at article outcomes.  Every
# factor remains distinguishable in greyscale through marker and line style.



CMAP_VESICLE = LinearSegmentedColormap.from_list(
    "vesicle_cmap", ["#EAF3FB", "#7FB3D9", "#2A6EBB", "#1B4F8A", "#0D2B4E"]
)
CMAP_SYNAPSIN = LinearSegmentedColormap.from_list(
    "synapsin_cmap", ["#FDF0E8", "#F0A070", "#D4541A", "#9B2C0A", "#4A1005"]
)
CMAP_BOTH = LinearSegmentedColormap.from_list(
    "both_cmap", ["#F3EEF8", "#B08DC8", "#5B2D8E", "#2E0F55"]
)


def apply_publication_style() -> None:
    """Apply the shared print-safe style without library-dependent defaults."""

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "axes.titleweight": "semibold",
            "axes.linewidth": 1.0,
            "axes.facecolor": COLORS["background"],
            "figure.facecolor": "white",
            "text.color": COLORS["text_dark"],
            "axes.labelcolor": COLORS["text_dark"],
            "axes.edgecolor": COLORS["text_dark"],
            "xtick.color": COLORS["text_dark"],
            "ytick.color": COLORS["text_dark"],
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.6,
            "grid.alpha": 0.65,
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def condition_style(labels: list[str]) -> Mapping[str, tuple[str, str, str]]:
    """Return deterministic colour, line, and marker styles for conditions."""

    return {
        label: (
            CONDITION_COLORS[index % len(CONDITION_COLORS)],
            LINE_STYLES[index % len(LINE_STYLES)],
            MARKERS[index % len(MARKERS)],
        )
        for index, label in enumerate(sorted(set(labels)))
    }


def _safe_stem(stem: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", stem.strip()).strip("._")
    if not cleaned:
        raise ValueError("Figure stem must contain at least one portable character")
    return cleaned


def save_figure(
    figure: plt.Figure,
    output_dir: str | Path,
    stem: str,
    *,
    show: bool = False,
    close: bool = True,
) -> Dict[str, Path]:
    """Export one Figure as vector PDF, 600-dpi TIFF, and 300-dpi PNG."""

    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    base = destination / _safe_stem(stem)
    paths = {
        "pdf": base.with_suffix(".pdf"),
        "tiff": base.with_suffix(".tif"),
        "png": base.with_suffix(".png"),
    }
    common = {"bbox_inches": "tight", "edgecolor": "none", "pad_inches": 0.04}
    figure.savefig(
        paths["pdf"], format="pdf", facecolor="white", transparent=False,
        metadata={"Creator": "SimulationSV Data Analysis Pipeline"}, **common,
    )
    figure.savefig(
        paths["tiff"], format="tiff", dpi=600, facecolor="white",
        transparent=False, pil_kwargs={"compression": "tiff_lzw"}, **common,
    )
    figure.savefig(
        paths["png"], format="png", dpi=300, facecolor="none",
        transparent=True, **common,
    )
    if show:
        plt.show()
    if close:
        plt.close(figure)
    return paths


apply_publication_style()
