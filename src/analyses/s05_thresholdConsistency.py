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
    Load abcd data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, site, PC10, threshold, prs, and ct information.
    """
    return util.load_data(
        "../../data/processed/abcd_data.npz",
        ["age", "sex", "site", "pc10", "thresh", "prs_all", "ct_vertex", "ct_aparc"],
    )


def load_multi_tle():
    """
    Load multisite TLE (ENIGMA) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, group, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/multi_tle_data.npz",
        ["age", "sex", "focus", "ct"],
    )
    
    
def load_multi_ige():
    """
    Load multisite IGE (ENIGMA) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, group, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/multi_ige_data.npz",
        ["age", "sex", "group", "ct"],
    )



# Main analysis
def main():
    # Load data
    age, sex, site, pc10, thresh, prs_all, ct_vertex, ct_aparc = load_abcd()

    print("------------------------------------------")
    print("Effects of PRS across different thresholds")
    print("------------------------------------------")

    thresholds = [
        "Pt_0.00100005",
        "Pt_0.0500001",
        "Pt_0.1",
        "Pt_0.2",
        "Pt_0.3",
        "Pt_0.4",
        "Pt_0.5",
    ]

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
    contrast = np.ones(prs_all.shape[0])
    slm = SLM(model, contrast)
    slm.fit(ct_vertex)
    residual = ct_vertex - np.dot(slm.X, slm.coef)

    # Lobe indices
    lobe_names = ["whole", "frontal", "limbic", "occipital", "parietal", "temporal"]
    frontal = [17, 18, 25, 16, 22, 2, 26, 10, 12, 30, 15, 46, 64, 60, 49, 44, 59, 52, 50, 56, 36, 51,]
    limbic = [8, 21, 1, 24, 33, 58, 35, 55, 42, 67]
    occipital = [11, 19, 9, 3, 45, 53, 37, 43]
    parietal = [20, 27, 29, 6, 23, 54, 61, 40, 63, 57]
    temporal = [4, 5, 7, 13, 14, 28, 31, 32, 0, 41, 47, 62, 38, 65, 48, 39, 66, 34]
    lobes = np.zeros(68)
    lobes[frontal] = 1
    lobes[limbic] = 2
    lobes[occipital] = 3
    lobes[parietal] = 4
    lobes[temporal] = 5
    lobes = parcel_to_surface(lobes, "aparc_conte69")

    # Correlate with PRS-TLE
    global_association = {}
    for i, lobe in enumerate(lobe_names):
        if "whole" in lobe:
            idx = np.array([True] * len(lobes))
        else:
            idx = lobes == i

        r = np.zeros(len(thresholds))
        r2 = np.zeros(len(thresholds))
        p = np.zeros(len(thresholds))
        mean_ct = np.mean(residual[:, idx], axis=1)
        for thr in thresholds:
            prs = prs_all[:, thresh == thr].flatten()
            corr = pearsonr(prs, mean_ct)
            r[i], r2[i], p[i] = corr.statistic, corr.statistic**2, corr.pvalue

        global_association[lobe] = pd.DataFrame(
            {"thresholds": thresholds, "r": r, "r2": r2, "p": p}
        )

        print(f"{lobe} lobe done")

    print()
    print("Regional effects across thresholds")
    print("------------------------------------------")
    regional_association = {}
    for thr in thresholds:
        prs = prs_all[:, thresh == thr]
        term_prs = FixedEffect(prs)
        model = term_age + term_sex + term_pc10 + term_prs
        contrast = prs
        slm = SLM(model, contrast, correction="fdr", two_tailed=True)
        slm.fit(ct_aparc)
        regional_association[thr] = slm

        print(f"{thr} done")

    print()
    print("Save results")
    print("------------------------------------------")
    util.save_to_pickle(
        "../../data/results/s05_thresholdConsistency/threshold_global_association.pkl",
        {"global_association": global_association},
    )

    util.save_to_pickle(
        "../../data/results/s05_thresholdConsistency/threshold_regional_association.pkl",
        {"regional_association": regional_association},
    )

    print()
    print("----------------------------------")
    print("Imaging-genetic pattern similarity")
    print("----------------------------------")

    print()
    print("Psnp thresholds - regional")
    print("----------------------------------")
    thr_similarity_r = np.zeros((len(thresholds), len(thresholds)))
    thr_similarity_p = np.zeros((len(thresholds), len(thresholds)))
    for i in range(len(thresholds)):
        thr1 = regional_association[thresholds[i]]
        for j in range(i, len(thresholds)):
            thr2 = regional_association[thresholds[j]]
            thr_similarity_r[i, j], thr_similarity_p[i, j], _ = (
                util.spatial_correlation(thr1.t, thr2.t, n_rot=5000)
            )
            print(f"{thresholds[i]} x {thresholds[j]} done")

    print()
    print("Comparison to case-control atrophy")
    print("----------------------------------")
    # Generate atrophy for each site:
    epilepsy_atrophy = {}
    
    multi_tle_atrophy = util.load_result(
        "../../data/results/03_atrophyAssociation/multi_atrophy.pkl", ["atrophy"]
    )[0]
    epilepsy_atrophy["ltle"] = multi_tle_atrophy["ltle"]
    epilepsy_atrophy["rtle"] = multi_tle_atrophy["rtle"]
    
    multi_ige_atrophy = util.load_result(
        "../../data/results/s02_igeSpecificity/ige_atrophy.pkl", ["atrophy"]
    )[0]
    epilepsy_atrophy["ige"] = multi_ige_atrophy["ige"]
    
    atrophy_similarity_r = np.zeros((len(thresholds), len(epilepsy_atrophy)))
    atrophy_similarity_p = np.zeros((len(thresholds), len(epilepsy_atrophy)))
    for i in range(len(thresholds)):
        map1 = regional_association[thresholds[i]]
        for j, subtype in enumerate(epilepsy_atrophy):
            map2 = epilepsy_atrophy[subtype]
            atrophy_similarity_r[i, j], atrophy_similarity_p[i, j], _ = util.spatial_correlation(map1.t, map2.t, n_rot=5000)


    print()
    print("Save results")
    print("----------------------------------")
    util.save_to_pickle(
        "../../data/results/s05_thresholdConsistency/threshold_similarity.pkl",
        {"thr_similarity_r": thr_similarity_r, "thr_similarity_p": thr_similarity_p},
    )

    util.save_to_pickle(
        "../../data/results/s05_thresholdConsistency/threshold_atrophy_similarity.pkl",
        {
            "atrophy_similarity_r": atrophy_similarity_r,
            "atrophy_similarity_p": atrophy_similarity_p,
        },
    )


if __name__ == "__main__":
    main()
