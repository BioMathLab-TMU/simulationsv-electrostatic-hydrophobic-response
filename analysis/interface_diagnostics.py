"""Post hoc geometric interface diagnostics for the fixed 12-run force design.

Reads published trajectories; never changes simulations or primary endpoints.
At shell=2 nm, an edge satisfies d_ij <= outer_radius_i+outer_radius_j+shell.
Syn with >=2 adjacent SVs is a geometric bridge candidate, not a measured bond.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def graph_metrics(positions, radii, types, shell=2.0, brute=False):
    positions, radii, types = map(np.asarray, (positions, radii, types))
    n = len(types)
    if positions.shape != (n, 3) or radii.shape != (n,):
        raise ValueError("inconsistent particle arrays")
    if not np.isfinite(positions).all() or not np.isfinite(radii).all():
        raise ValueError("nonfinite particle data")
    if set(np.unique(types)) != {1, 2}:
        raise ValueError("expected Syn=1 and SV=2")
    if brute:
        i, j = np.triu_indices(n, 1)
        pairs = np.column_stack((i, j))
    else:
        pairs = cKDTree(positions).query_pairs(2 * radii.max() + shell,
                                              output_type="ndarray")
    if len(pairs):
        d = np.linalg.norm(positions[pairs[:, 0]] - positions[pairs[:, 1]], axis=1)
        pairs = pairs[d <= radii[pairs[:, 0]] + radii[pairs[:, 1]] + shell]
    left, right = pairs.T if len(pairs) else (np.array([], int), np.array([], int))
    ss = (types[left] == 1) & (types[right] == 1)
    vv = (types[left] == 2) & (types[right] == 2)
    sv = types[left] != types[right]
    syn_vertex = np.where(types[left[sv]] == 1, left[sv], right[sv])
    ves_vertex = np.where(types[left[sv]] == 2, left[sv], right[sv])
    syn_sv_degree = np.bincount(syn_vertex, minlength=n)[types == 1]
    ves_syn_degree = np.bincount(ves_vertex, minlength=n)[types == 2]
    degree = np.bincount(np.concatenate((left, right)), minlength=n)
    adjacency = coo_matrix((np.ones(2*len(left)),
                            (np.concatenate((left, right)),
                             np.concatenate((right, left)))), shape=(n, n)).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    # Same declared tie rule as the spatial pipeline: largest size, then lowest
    # node index encountered by scipy's component labeling.
    largest = np.argmax(np.bincount(labels))
    return {
        "edges_SS": int(ss.sum()), "edges_SV": int(sv.sum()),
        "edges_VV": int(vv.sum()), "edges_all": int(len(left)),
        "syn_near_sv_fraction": float(np.mean(syn_sv_degree >= 1)),
        "syn_multivesicle_neighbor_fraction": float(np.mean(syn_sv_degree >= 2)),
        "sv_near_syn_fraction": float(np.mean(ves_syn_degree >= 1)),
        "syn_sv_degree_mean": float(np.mean(syn_sv_degree)),
        "sv_syn_degree_mean": float(np.mean(ves_syn_degree)),
        "syn_nonisolated_fraction": float(np.mean(degree[types == 1] >= 1)),
        "largest_component_fraction": float(np.mean(labels == largest)),
        "largest_component_sv_fraction": float(np.mean(labels[types == 2] == largest)),
    }


def self_test():
    # One Syn touches two SVs; the SVs do not touch each other.
    p = np.array([[0, 0, 0], [2, 0, 0], [-2, 0, 0]], float)
    r = np.ones(3)
    t = np.array([1, 2, 2])
    got = graph_metrics(p, r, t, shell=0)
    assert got["edges_SV"] == 2 and got["edges_VV"] == 0
    assert got["syn_multivesicle_neighbor_fraction"] == 1
    assert got["largest_component_fraction"] == 1
    assert got == graph_metrics(p, r, t, shell=0, brute=True)
    # Boundary inclusion, translation/permutation invariance, zero-edge case.
    translated = graph_metrics(p+np.array([51, -23, 8]), r, t, shell=0)
    assert translated == got
    order = np.array([2, 0, 1])
    assert graph_metrics(p[order], r[order], t[order], shell=0) == got
    isolated = graph_metrics(10*p, r, t, shell=0)
    assert isolated["edges_all"] == 0
    assert isolated["syn_near_sv_fraction"] == 0
    assert isolated["largest_component_fraction"] == 1/3
    rng = np.random.default_rng(4821)
    for _ in range(8):
        p = rng.uniform(-15, 15, size=(40, 3))
        r = rng.uniform(.5, 3, size=40)
        t = np.tile([1, 2], 20)
        assert graph_metrics(p, r, t) == graph_metrics(p, r, t, brute=True)


def render(run_table, figures):
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42, "ps.fonttype": 42})
    arms = [(0, 0), (0, 1), (1, 0), (1, 1)]
    colors = ["#737373", "#E76F51", "#2A9D8F", "#195E56"]
    markers = ["o", "s", "^", "D"]
    columns = ["syn_near_sv_fraction", "syn_multivesicle_neighbor_fraction", "edges_SV"]
    labels = ["Syn adjacent to at least one SV", "Syn adjacent to at least two SVs", "Syn-SV edges"]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.6), layout="constrained")
    for k, (ax, column, label) in enumerate(zip(axes, columns, labels)):
        for x, ((e, h), color, marker) in enumerate(zip(arms, colors, markers)):
            values = run_table[(run_table.E == e) & (run_table.H == h)][column].to_numpy()
            ax.scatter(x + np.linspace(-.09, .09, len(values)), values,
                       edgecolor=color, facecolor="white", marker=marker, s=37, zorder=3)
            ax.errorbar(x, values.mean(), yerr=values.std(ddof=1), color=color,
                        marker="_", markersize=18, capsize=4, lw=1.5)
        ax.set(xticks=range(4), xticklabels=["E0H0", "E0H1", "E1H0", "E1H1"], ylabel=label)
        ax.set_title(chr(65+k), loc="left", fontweight="bold")
        ax.set_ylim(bottom=0)
        if k < 2:
            ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.15)
    figures.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "tif"):
        kwargs = {"dpi": 600} if ext == "tif" else {"dpi": 180}
        if ext == "tif":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(figures / f"force_interface_diagnostics.{ext}", **kwargs)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path)
    args = parser.parse_args()
    self_test()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "configs/manifest.json").read_text())
    availability = json.loads((root / "data/availability.json").read_text())
    expected = {r["run_id"]: r for r in availability["runs"]}
    old = pd.read_csv(root / "data/processed/run_endpoints.csv").set_index("run_id")
    frames, runs, provenance = [], [], []
    for spec in manifest["runs"]:
        run_id = spec["run_id"]
        path = root / "data/raw" / run_id / "trajectory.h5"
        sha = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
        if sha != expected[run_id]["publication_h5_sha256"]:
            raise ValueError("trajectory checksum mismatch: " + run_id)
        local = []
        with h5py.File(path, "r") as h5:
            names = sorted((k for k in h5 if k.startswith("step_")), key=lambda k: int(k[5:]))
            assert len(names) == 200
            for index, name in enumerate(names):
                frame = h5[name]
                p, r, t = (frame[k][:] for k in ("positions", "radii", "types"))
                # Use float64 for distance arithmetic, without changing stored values.
                p, r = p.astype(float), r.astype(float)
                metrics = graph_metrics(p, r, t)
                if index in (0, 199):
                    assert metrics == graph_metrics(p, r, t, brute=True)
                row = {"run_id": run_id, "seed": spec["seed"],
                       "E": int(spec["electrostatic_enabled"]),
                       "H": int(spec["hydrophobic_enabled"]),
                       "step": int(name[5:]),
                       "time_ms": float(frame["physical_time_ms"][()].item()), **metrics}
                assert metrics["edges_all"] == sum(metrics[k] for k in ("edges_SS", "edges_SV", "edges_VV"))
                assert metrics["syn_multivesicle_neighbor_fraction"] <= metrics["syn_near_sv_fraction"] <= metrics["syn_nonisolated_fraction"]
                local.append(row)
        table = pd.DataFrame(local)
        late = table.iloc[-60:]
        endpoints = late[list(metrics)].mean().to_dict()
        for new, previous in [("largest_component_fraction", "late_mean_largest_component_fraction"),
                              ("largest_component_sv_fraction", "late_mean_largest_component_vesicle_fraction"),
                              ("syn_nonisolated_fraction", "late_mean_synapsin_recruited_fraction")]:
            np.testing.assert_allclose(endpoints[new], old.loc[run_id, previous], rtol=0, atol=1e-12)
        runs.append({k: row[k] for k in ("run_id", "seed", "E", "H")} | endpoints)
        frames.extend(local)
        provenance.append({"run_id": run_id, "sha256": sha, "frames": 200})
        print("INTERFACE_DONE", run_id, flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    run_table = pd.DataFrame(runs)
    run_table.to_csv(args.output / "interface_run_metrics.csv", index=False)
    pd.DataFrame(frames).to_csv(args.output / "interface_frame_metrics.csv", index=False)
    summary = run_table.groupby(["E", "H"]).agg({k: ["mean", "std"] for k in metrics})
    summary.to_csv(args.output / "interface_condition_summary.csv")
    # Independently recompute the registered primary contrast; no statistical helper reused.
    primary = run_table.pivot(index="seed", columns=["E", "H"], values="largest_component_fraction")
    interaction = primary[1, 1] - primary[1, 0] - primary[0, 1] + primary[0, 0]
    audit = {"schema_version": 1, "shell_nm": 2, "late_frames": 60,
             "analysis_status": "post_hoc_descriptive", "n_runs": 12, "n_seed_blocks": 3,
             "primary_interaction_by_seed": interaction.to_dict(),
             "primary_interaction_mean": float(interaction.mean()),
             "tests": {"toy_and_random_bruteforce": "passed", "first_last_frame_bruteforce_each_run": "passed",
                       "all_12_primary_endpoint_reconciliations": "passed", "trajectory_hashes": "passed"},
             "interpretation": "Snapshot proximity is not a biochemical bond, bridge lifetime or pair-class causal ablation.",
             "inputs": provenance}
    (args.output / "interface_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    render(run_table, args.figures or args.output / "figures")
    print(summary.to_string())


if __name__ == "__main__":
    main()
