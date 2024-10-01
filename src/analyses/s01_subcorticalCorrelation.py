import numpy as np
import pandas as pd
import utilities as util
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM
from scipy.stats import pearsonr


# Helper functions
def load_abcd():
    """
    Load ABCD data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, site, PC10, threshold, prs, and sv information.
    """
    return util.load_data(
        "../../data/processed/abcd_data.npz",
        ["age", "sex", "site", "pc10", "thresh", "prs_all", "sv", "icv"],
    )


# Main analysis
def main():
    # Load data
    age, sex, site, pc10, thresh, prs_all, sv, icv = load_abcd()

    thr = "Pt_0.1"
    prs = prs_all[:, thresh == thr].flatten()

    print("------------------------------")
    print("Subcortical effects of PRS-TLE")
    print("------------------------------")

    print()
    print("Correlate PRS with each region")
    print("------------------------------")
    # Regress out effects of age, sex, and term_pc10
    term_age = FixedEffect(age, "Age")
    term_sex = FixedEffect(sex, "Sex")
    term_pc10 = FixedEffect(pc10, [f"PC{i}" for i in range(10)])
    term_icv = FixedEffect(icv, "ICV")
    model = term_age + term_sex + term_pc10 + term_icv
    contrast = np.ones(len(prs))
    slm = SLM(model, contrast)
    slm.fit(sv)
    residual = sv - np.dot(slm.X, slm.coef)

    n_sctx = sv.shape[1]
    # Lobe indices
    subcortical = {
        "lh_whole": np.array([True] * int(n_sctx / 2) + [False] * int(n_sctx / 2)),
        "lh_accumbens": 0,
        "lh_amygdala": 1,
        "lh_caudate": 2,
        "lh_hippocampus": 3,
        "lh_pallidum": 4,
        "lh_putamen": 5,
        "lh_thalamus": 6,
        "rh_whole": np.array([False] * int(n_sctx / 2) + [True] * int(n_sctx / 2)),
        "rh_accumbens": 7,
        "rh_amygdala": 8,
        "rh_caudate": 9,
        "rh_hippocampus": 10,
        "rh_pallidum": 11,
        "rh_putamen": 12,
        "rh_thalamus": 13,
    }

    # Correlate with PRS-TLE
    mean_sv = np.zeros((len(prs), len(subcortical)))
    r = np.zeros(len(subcortical))
    r2 = np.zeros(len(subcortical))
    p = np.zeros(len(subcortical))

    for i, sctx in enumerate(subcortical):
        idx = subcortical[sctx]
        mean_sv[:, i] = (
            np.sum(residual[:, idx], axis=1) if "whole" in sctx else residual[:, idx]
        )
        corr = pearsonr(prs, mean_sv[:, i])
        r[i], r2[i], p[i] = corr.statistic, corr.statistic**2, corr.pvalue
        print(f"{sctx}: r = {r[i]}, r-squared = {r2[i]}, p = {p[i]}")

    global_association = pd.DataFrame(
        {"subcortical": subcortical, "r": r, "r2": r2, "p": p}
    )

    print()
    print("Save results")
    print("-----------------------------")
    util.save_to_pickle(
        "../../data/results/s01_subcorticalCorrelation/global_association.pkl",
        {"global_association": global_association, "mean_sv": mean_sv, "prs": prs},
    )
