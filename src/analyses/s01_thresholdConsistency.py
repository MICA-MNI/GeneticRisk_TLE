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
    numpy.ndarray: The loaded data containing age, sex, dataset, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/abcd_data.npz",
        ["age", "sex", "site", "pc10", "thresh", "prs_all", "ct_vertex", "ct_aparc"],
    )


def load_local_tle():
    """
    Load local TLE (EpiC x MICs x NKG) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/local_tle_data.npz",
        ["age", "sex", "dataset", "focus", "ct"],
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
    frontal = [
        17,
        18,
        25,
        16,
        22,
        2,
        26,
        10,
        12,
        30,
        15,
        46,
        64,
        60,
        49,
        44,
        59,
        52,
        50,
        56,
        36,
        51,
    ]
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
        "../../data/results/s01_thresholdConsistency/threshold_global_association.pkl",
        {"global_association": global_association},
    )

    util.save_to_pickle(
        "../../data/results/s01_thresholdConsistency/threshold_regional_association.pkl",
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
    print("Psnp thresholds - epicentres")
    print("----------------------------------")
    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()

    thr_fc_epicentre = np.zeros((len(thresholds), 82))
    thr_sc_epicentre = np.zeros((len(thresholds), 82))
    fc_epicentre_similarity_r = np.zeros((len(thresholds), len(thresholds)))
    fc_epicentre_similarity_p = np.zeros((len(thresholds), len(thresholds)))
    sc_epicentre_similarity_r = np.zeros((len(thresholds), len(thresholds)))
    sc_epicentre_similarity_p = np.zeros((len(thresholds), len(thresholds)))
    for i in range(len(thresholds)):
        thr1 = regional_association[thresholds[i]].t

        fc_ctx_r1, fc_ctx_p1 = util.epicenter_mapping(thr1, fc_ctx)
        fc_sctx_r1, fc_sctx_p1 = util.epicenter_mapping(thr1, fc_sctx)
        sc_ctx_r1, sc_ctx_p1 = util.epicenter_mapping(thr1, sc_ctx)
        sc_sctx_r1, sc_sctx_p1 = util.epicenter_mapping(thr1, sc_sctx)

        thr_fc_epicentre[i, :] = np.concatenate((fc_ctx_r1, fc_sctx_r1))
        thr_sc_epicentre[i, :] = np.concatenate((sc_ctx_r1, sc_sctx_r1))

        for j in range(i, len(thresholds)):
            thr2 = regional_association[thresholds[j]].t

            fc_ctx_r2, fc_ctx_p2 = util.epicenter_mapping(thr2, fc_ctx)
            fc_sctx_r2, fc_sctx_p2 = util.epicenter_mapping(thr2, fc_sctx)
            sc_ctx_r2, sc_ctx_p2 = util.epicenter_mapping(thr2, sc_ctx)
            sc_sctx_r2, sc_sctx_p2 = util.epicenter_mapping(thr2, sc_sctx)

            fc_epicentre_similarity_r[i, j], fc_epicentre_similarity_p[i, j], _ = (
                util.spatial_correlation(
                    np.concatenate((fc_ctx_r1, fc_sctx_r1)),
                    np.concatenate((fc_ctx_r2, fc_sctx_r2)),
                    surface_name="fsa5_with_sctx",
                    parcellation_name="aparc_aseg",
                    n_rot=5000,
                )
            )

            sc_epicentre_similarity_r[i, j], sc_epicentre_similarity_p[i, j], _ = (
                util.spatial_correlation(
                    np.concatenate((sc_ctx_r1, sc_sctx_r1)),
                    np.concatenate((sc_ctx_r2, sc_sctx_r2)),
                    surface_name="fsa5_with_sctx",
                    parcellation_name="aparc_aseg",
                    n_rot=5000,
                )
            )

            print(f"{thresholds[i]} x {thresholds[j]} done")

    print()
    print("Comparison to case-control atrophy")
    print("----------------------------------")
    # Generate atrophy for each site:
    age, sex, dataset, focus, ct = load_local_tle()
    local_tle_atrophy = util.load_data(
        "../../data/results/03_atrophyAssociation/local_atrophy.pkl", ["atrophy"]
    )
    # Combined
    tle_atrophy = {}
    tle_atrophy["lh_all"] = local_tle_atrophy["ltle"]["slm"].t
    tle_atrophy["rh_all"] = local_tle_atrophy["rtle"]["slm"].t

    # Across different datasets
    for d in ["EpiC", "MICs", "NKG"]:
        dataset_idx = dataset == d
        dataset_ct = ct[dataset_idx, :]
        dataset_age = age[dataset_idx]
        dataset_sex = sex[dataset_idx]
        dataset_focus = focus[dataset_idx]

        x = util.zscore_flip(dataset_ct, dataset_focus, "C", "R")
        for hemi in ["L", "R"]:
            covar = pd.DataFrame({"age": dataset_age, "sex": dataset_sex})
            tle_atrophy[f"{hemi.lower()}h_{d.lower()}"] = util.casecontrol_difference(
                x, covar, dataset_focus, "C", hemi
            ).t

        print(f"{d} atrophy done")

    site_labels = [
        "lh_all",
        "lh_epic",
        "lh_mics",
        "lh_nkg",
        "rh_all",
        "rh_epic",
        "rh_mics",
        "rh_nkg",
    ]
    atrophy_similarity_r = np.zeros((len(thresholds), len(site_labels)))
    atrophy_similarity_p = np.zeros((len(thresholds), len(site_labels)))
    for i in range(len(thresholds)):
        map1 = regional_association[thresholds[i]]
        for site in site_labels:
            map2 = tle_atrophy[site]
            atrophy_similarity_r[i, j], atrophy_similarity_p[i, j], _ = (
                spatial_correlation(map1.t, map2.t, n_rot=5000)
            )

    print()
    print("Save results")
    print("----------------------------------")
    util.save_to_pickle(
        "../../data/results/s01_thresholdConsistency/threshold_similarity.pkl",
        {"thr_similarity_r": thr_similarity_r, "thr_similarity_p": thr_similarity_p},
    )

    util.save_to_pickle(
        "../../data/results/s01_thresholdConsistency/epicentre_similarity.pkl",
        {
            "fc_epicentre_similarity_r": fc_epicentre_similarity_r,
            "fc_epicentre_similarity_p": fc_epicentre_similarity_p,
            "sc_epicentre_similarity_r": sc_epicentre_similarity_r,
            "sc_epicentre_similarity_p": sc_epicentre_similarity_p,
        },
    )


if __name__ == "__main__":
    main()
