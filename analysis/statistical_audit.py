"""Independent small-sample audit of the fixed 2x2, three-block design."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
METRICS = ["late_mean_largest_component_fraction",
           "late_mean_synapsin_recruited_fraction",
           "late_mean_largest_component_vesicle_fraction"]


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.minimum(1, np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
    result = np.empty_like(p)
    result[order] = adjusted
    return result


def main():
    assert np.allclose(holm([.01, .04, .03]), [.03, .06, .06])
    data = pd.read_csv(ROOT/"data/processed/run_endpoints.csv").sort_values(["seed", "E", "H"])
    paired = pd.read_csv(ROOT/"data/processed/paired_contrasts.csv")
    if len(data) != 12 or data.groupby(["E", "H"]).size().tolist() != [3]*4:
        raise ValueError("incomplete design")
    seeds = sorted(data.seed.unique())
    X = np.column_stack([np.ones(len(data)), data.E, data.H, data.E*data.H,
                         data.seed.eq(seeds[1]), data.seed.eq(seeds[2])]).astype(float)
    if np.linalg.matrix_rank(X) != 6:
        raise ValueError("rank-deficient design")
    inverse = np.linalg.inv(X.T@X)
    leverage = np.sum((X@inverse)*X, axis=1)
    fits, contrasts = [], []
    names = ["intercept", "E_at_H0", "H_at_E0", "E_by_H", "seed_2", "seed_3"]
    for metric in METRICS:
        y = data[metric].to_numpy(float)
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        residual = y-X@beta
        weighted = residual**2/(1-leverage)**2
        covariance = inverse@(X.T@(X*weighted[:, None]))@inverse
        se = np.sqrt(np.maximum(0, np.diag(covariance)))
        probabilities = 2*t.sf(np.abs(beta/se), df=6)
        adjusted = holm(probabilities)
        for j, name in enumerate(names):
            fits.append(dict(metric=metric, term=name, estimate=beta[j], se_HC3=se[j],
                             ci95_low=beta[j]-t.ppf(.975, 6)*se[j],
                             ci95_high=beta[j]+t.ppf(.975, 6)*se[j],
                             p_raw=probabilities[j], p_holm=adjusted[j], residual_df=6))
        table = data.pivot(index="seed", columns=["E", "H"], values=metric)
        effects = {"E_at_H0": table[1, 0]-table[0, 0],
                   "E_at_H1": table[1, 1]-table[0, 1],
                   "H_at_E0": table[0, 1]-table[0, 0],
                   "H_at_E1": table[1, 1]-table[1, 0],
                   "E_by_H": table[1, 1]-table[1, 0]-table[0, 1]+table[0, 0]}
        for name, values in effects.items():
            old = paired[(paired.metric == metric) & (paired.contrast == name)]
            # The existing pipeline uses interaction_DID for this term.
            if old.empty and name == "E_by_H":
                old = paired[(paired.metric == metric) & paired.contrast.str.contains("interaction|DID", case=False)]
            mean, sd = float(values.mean()), float(values.std(ddof=1))
            se_d = sd/np.sqrt(3)
            lo, hi = mean-t.ppf(.975, 2)*se_d, mean+t.ppf(.975, 2)*se_d
            if len(old) != 1:
                raise ValueError("missing published contrast: "+name)
            np.testing.assert_allclose([mean, lo, hi], old[["mean", "ci95_low", "ci95_high"]].iloc[0], atol=1e-12)
            contrasts.append(dict(metric=metric, contrast=name, mean=mean, sd=sd,
                                  ci95_low=lo, ci95_high=hi,
                                  paired_t_p_raw=2*t.sf(abs(mean/se_d), 2),
                                  standardized_dz=mean/sd))
    out = ROOT/"results/statistical_audit"
    out.mkdir(parents=True, exist_ok=True)
    c = pd.DataFrame(contrasts)
    c["paired_t_p_holm"] = c.groupby("metric")["paired_t_p_raw"].transform(lambda v: holm(v))
    c.to_csv(out/"paired_parametric.csv", index=False)
    pd.DataFrame(fits).to_csv(out/"factorial_HC3.csv", index=False)
    (out/"audit.json").write_text(json.dumps({"published_contrast_checks": 15,
        "absolute_tolerance": 1e-12, "replicate": "complete run / seed block",
        "caution": "Three blocks; HC3 is approximate, not cluster-robust; no nonparametric p<0.25 is possible here.",
        "multiplicity": "Holm within each outcome and contrast family; factorial terms within outcome"}, indent=2))
    print("15 published paired contrasts independently reconciled; HC3 outputs saved")


if __name__ == "__main__":
    main()
