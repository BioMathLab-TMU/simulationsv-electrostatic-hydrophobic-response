"""Reproduce the post hoc xy-TAMSD diagnostic from released trajectories.

No simulation is changed, no negative fitted slopes are clipped, and the
experimental scale reference is not a matched validation or fitting target.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def fitted_mobility(xy, spacing_s):
    lags = np.arange(int(np.ceil(.002 / spacing_s)),
                     int(np.floor(.004 / spacing_s)) + 1)
    lag_s = lags * spacing_s
    if len(lags) < 2 or lags[-1] >= len(xy):
        raise ValueError("insufficient track support for the declared fit")
    msd = np.array([np.mean(np.sum((xy[k:] - xy[:-k])**2, axis=2), axis=0)
                    for k in lags])
    centered = lag_s - lag_s.mean()
    estimates = (centered @ msd / np.sum(centered**2)) / 4
    np.testing.assert_allclose(estimates.mean(),
                               np.polyfit(lag_s, msd.mean(axis=1), 1)[0]/4,
                               atol=1e-12, rtol=1e-12)
    return estimates, lag_s, msd.mean(axis=1)


def self_test():
    xy = np.zeros((200, 3, 2))
    np.testing.assert_allclose(fitted_mobility(xy, .00014)[0], 0)
    xy[:, :, 0] = np.arange(200)[:, None]*.00014*2
    d, lags, _ = fitted_mobility(xy, .00014)
    np.testing.assert_allclose(d, .5*2**2*lags.mean(), rtol=1e-12)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    self_test()
    root = Path(__file__).resolve().parents[1]
    availability = json.loads((root/"data/availability.json").read_text())
    rows, agents, curves, inputs = [], [], [], []
    for spec in availability["runs"]:
        path = root/"data/raw"/spec["run_id"]/"trajectory.h5"
        with path.open("rb") as handle:
            sha = hashlib.file_digest(handle, "sha256").hexdigest()
        if sha != spec["publication_h5_sha256"]:
            raise ValueError("checksum mismatch: "+spec["run_id"])
        with h5py.File(path, "r") as f:
            frames = sorted((k for k in f if re.fullmatch(r"step_\d+", k)),
                            key=lambda k: int(k[5:]))
            if len(frames) != 200:
                raise ValueError("expected 200 frames")
            ids, types = f[frames[0]+"/ids"][:], f[frames[0]+"/types"][:]
            syn = types == 1
            if syn.sum() != 1807:
                raise ValueError("wrong Syn population")
            for k in frames:
                np.testing.assert_array_equal(f[k+"/ids"][:], ids)
                np.testing.assert_array_equal(f[k+"/types"][:], types)
            time_ms = np.array([float(f[k+"/physical_time_ms"][0]) for k in frames])
            np.testing.assert_allclose(np.diff(time_ms), .14, atol=1e-8)
            spacing_s = float(np.median(np.diff(time_ms)))/1000
            xy = np.array([f[k+"/positions"][:][syn, :2] for k in frames],
                          dtype=np.float64)/1000
        identity = dict(arm=f"E{int(spec['E'])}H{int(spec['H'])}", seed=spec["seed"])
        for window, coordinates in (("full_28ms", xy), ("late_60frames", xy[-60:])):
            ds, lags, msd = fitted_mobility(coordinates, spacing_s)
            if not np.isfinite(ds).all():
                raise ValueError("nonfinite fitted slopes")
            rows.append(dict(**identity, window=window, n_agents=len(ds),
                             D_mean_um2_s=float(ds.mean()),
                             D_median_um2_s=float(np.median(ds)),
                             negative_slope_fraction=float(np.mean(ds < 0)),
                             min_lag_ms=1000*lags.min(), max_lag_ms=1000*lags.max()))
            agents.extend(dict(**identity, window=window, agent_id=int(i), D_um2_s=float(d))
                          for i, d in zip(ids[syn], ds))
            curves.extend(dict(**identity, window=window, lag_ms=1000*l, msd_xy_um2=float(m))
                          for l, m in zip(lags, msd))
        inputs.append(dict(run_id=spec["run_id"], sha256=sha))
    a.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output/"run_level_mobility.csv", index=False)
    pd.DataFrame(agents).to_csv(a.output/"agent_level_mobility.csv", index=False)
    pd.DataFrame(curves).to_csv(a.output/"lag_msd_xy.csv", index=False)
    (a.output/"audit.json").write_text(json.dumps(dict(
        analysis="post_hoc_scale_diagnostic_not_matched_validation",
        tests="stationary, constant velocity, independent ensemble slope, identity and SHA256 checks passed",
        units=dict(position_input="nm", fitted_position="um", fit_time="s", D="um^2/s"),
        reference="data/external/hoffmann_mobility_reference.json",
        negative_slopes="retained", inputs=inputs), indent=2), encoding="utf-8")
    print(pd.DataFrame(rows).groupby(["arm", "window"])["D_mean_um2_s"].agg(["mean", "std"]))


if __name__ == "__main__":
    main()
