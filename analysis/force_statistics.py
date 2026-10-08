"""Paired, seed-level statistics for the registered 2x2 force experiment."""
from __future__ import annotations

import itertools
import math
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

ENDPOINTS = {
    "late_mean_largest_component_fraction": "primary",
    "late_mean_synapsin_recruited_fraction": "secondary",
    "late_mean_largest_component_vesicle_fraction": "secondary",
}
CONTRASTS = {
    "E_at_H0": {(1, 0): 1, (0, 0): -1},
    "E_at_H1": {(1, 1): 1, (0, 1): -1},
    "H_at_E0": {(0, 1): 1, (0, 0): -1},
    "H_at_E1": {(1, 1): 1, (1, 0): -1},
    "E_H_interaction": {(1, 1): 1, (1, 0): -1, (0, 1): -1, (0, 0): 1},
}


def mean_ci(values: Sequence[float]) -> dict[str, float | int]:
    x = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError("Nonfinite values cannot enter inference")
    n = int(len(x))
    if not n:
        raise ValueError("No observations")
    mean = float(x.mean())
    sd = float(x.std(ddof=1)) if n > 1 else math.nan
    half = float(student_t.ppf(0.975, n - 1) * sd / math.sqrt(n)) if n > 1 else math.nan
    return {"n_seeds": n, "mean": mean, "sd": sd,
            "ci95_low": mean - half, "ci95_high": mean + half}


def exact_sign_flip(values: Sequence[float]) -> float:
    x = np.asarray(values, dtype=float)
    if not len(x) or not np.all(np.isfinite(x)):
        raise ValueError("Require finite paired differences")
    observed = abs(float(x.mean()))
    return sum(abs(float((x * signs).mean())) >= observed - 1e-15
               for signs in itertools.product((-1, 1), repeat=len(x))) / 2 ** len(x)


def holm(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(p[order] * (len(p) - np.arange(len(p)))))
    return adjusted


def paired_statistics(run_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["seed", "electrostatic_enabled", "hydrophobic_enabled"]
    if run_level.duplicated(keys).any():
        raise ValueError("Duplicate seed/arm: do not pool repeated extractions")
    summaries, contrasts = [], []
    for metric, role in ENDPOINTS.items():
        for (e, h), group in run_level.groupby(keys[1:]):
            summaries.append({"metric": metric, "role": role, "E": int(e), "H": int(h),
                              **mean_ci(group[metric].tolist())})
        pivot = run_level.pivot(index="seed", columns=keys[1:], values=metric)
        required = {(0, 0), (0, 1), (1, 0), (1, 1)}
        if not required.issubset(set(pivot.columns)) or pivot.isna().any().any():
            raise ValueError("Paired factorial inference requires every arm of every seed")
        for name, weights in CONTRASTS.items():
            differences = sum(weight * pivot[arm] for arm, weight in weights.items())
            contrasts.append({"metric": metric, "role": role, "contrast": name,
                              **mean_ci(differences.tolist()),
                              "exact_sign_flip_p": exact_sign_flip(differences.tolist()),
                              "paired_seed_differences": ";".join(f"{seed}:{value:.17g}" for seed, value in differences.items()),
                              "inference_scope": "model endpoint; sign-flip requires symmetry; n=3 is low resolution"})
    contrast_frame = pd.DataFrame(contrasts)
    for _, indices in contrast_frame.groupby("metric").groups.items():
        contrast_frame.loc[indices, "exact_sign_flip_p_holm"] = holm(contrast_frame.loc[indices, "exact_sign_flip_p"])
    return pd.DataFrame(summaries), contrast_frame
