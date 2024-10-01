import numpy as np
import pandas as pd
import utilities as util
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM
from enigmatoolbox.utils import parcel_to_surface
from scipy.stats import pearsonr


# Helper functions
def load_abcd():
    """
    Load ABCD data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, site, PC10, threshold, prs, and ct information.
    """
    return util.load_data(
        "../../data/processed/abcd_data.npz",
        ["age", "sex", "site", "pc10", "thresh", "prs_all", "ct_vertex", "ct_aparc"],
    )


# Main analysis
def main():
    # Load data
    age, sex, site, pc10, thresh, prs_all, ct_vertex, ct_aparc = load_abcd()

    thr = "Pt_0.1"
    prs = prs_all[:, thresh == thr].flatten()

    print("-----------------------------")
    print("Brain-wide effects of PRS-TLE")
    print("-----------------------------")

    print()
    print("Correlate PRS with each lobe")
    print("-----------------------------")
    # Regress out effects of age, sex, and term_pc10
    term_age = FixedEffect(age, "Age")
    term_sex = FixedEffect(sex, "Sex")
    term_pc10 = FixedEffect(pc10, [f"PC{i}" for i in range(10)])
    model = term_age + term_sex + term_pc10
    contrast = np.ones(len(prs))
    slm = SLM(model, contrast)
    slm.fit(ct_vertex)
    residual = ct_vertex - np.dot(slm.X, slm.coef)

    # Lobe indices
    lobe_names = [
        "lh_whole",
        "lh_frontal",
        "lh_limbic",
        "lh_occipital",
        "lh_parietal",
        "lh_temporal",
        "rh_whole",
        "rh_frontal",
        "rh_limbic",
        "rh_occipital",
        "rh_parietal",
        "rh_temporal",
    ]
    lh_frontal = [17, 18, 25, 16, 22, 2, 26, 10, 12, 30, 15]
    lh_limbic = [8, 21, 1, 24, 33]
    lh_occipital = [11, 19, 9, 3]
    lh_parietal = [20, 27, 29, 6, 23]
    lh_temporal = [4, 5, 7, 13, 14, 28, 31, 32, 0]
    rh_frontal = [46, 64, 60, 49, 44, 59, 52, 50, 56, 36, 51]
    rh_limbic = [58, 35, 55, 42, 67]
    rh_occipital = [45, 53, 37, 43]
    rh_parietal = [54, 61, 40, 63, 57]
    rh_temporal = [41, 47, 62, 38, 65, 48, 39, 66, 34]
    lobes = np.zeros(68)
    lobes[lh_frontal] = 1
    lobes[lh_limbic] = 2
    lobes[lh_occipital] = 3
    lobes[lh_parietal] = 4
    lobes[lh_temporal] = 5
    lobes[rh_frontal] = 7
    lobes[rh_limbic] = 8
    lobes[rh_occipital] = 9
    lobes[rh_parietal] = 10
    lobes[rh_temporal] = 11
    lobes = parcel_to_surface(lobes, "aparc_conte69")

    # Correlate with PRS-TLE
    mean_ct = np.zeros((len(prs), len(lobe_names)))
    r = np.zeros(len(lobe_names))
    r2 = np.zeros(len(lobe_names))
    p = np.zeros(len(lobe_names))

    for i, lobe in enumerate(lobe_names):
        if "whole" in lobe:
            idx = (
                np.array([True] * int(len(lobes) / 2) + [False] * int(len(lobes) / 2))
                if "lh" in lobe
                else np.array(
                    [False] * int(len(lobes) / 2) + [True] * int(len(lobes) / 2)
                )
            )
        else:
            idx = lobes == i

        mean_ct[:, i] = np.mean(residual[:, idx], axis=1)
        corr = pearsonr(prs, mean_ct[:, i])
        r[i], r2[i], p[i] = corr.statistic, corr.statistic**2, corr.pvalue
        print(f"{lobe}: r = {r[i]}, r-squared = {r2[i]}, p = {p[i]}")

    global_association = pd.DataFrame({"lobe": lobe_names, "r": r, "r2": r2, "p": p})

    print()
    print("Save results")
    print("-----------------------------")
    util.save_to_pickle(
        "../../data/results/01_geneticCorrelation/global_association.pkl",
        {"global_association": global_association, "mean_ct": mean_ct, "prs": prs},
    )

    print()
    print("-------------------------")
    print("Local effects of PRS-TLE")
    print("------------------------")

    print()
    print("Correlate across regions")
    print("------------------------")
    term_prs = FixedEffect(prs)
    model = term_age + term_sex + term_pc10 + term_prs
    contrast = prs
    slm = SLM(model, contrast, correction="fdr", two_tailed=True)
    slm.fit(ct_aparc)

    print()
    print("Save results")
    print("-----------------------------")
    util.save_to_pickle(
        "../../data/results/01_geneticCorrelation/regional_association.pkl",
        {"slm": slm},
    )


if __name__ == "__main__":
    main()
