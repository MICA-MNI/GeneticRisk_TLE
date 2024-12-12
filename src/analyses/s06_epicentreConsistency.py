import numpy as np
import utilities as util


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


# Main analysis
def main():
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

    regional_association = util.load_result(
        "../../data/results/s01_thresholdConsistency/threshold_regional_association.pkl",
        ["regional_association"],
    )

    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()

    thr_fc_epicentre = np.zeros((len(thresholds), 82))
    thr_sc_epicentre = np.zeros((len(thresholds), 82))
    fc_epicentre_similarity_r = np.zeros((len(thresholds), len(thresholds)))
    fc_epicentre_similarity_p = np.zeros((len(thresholds), len(thresholds)))
    sc_epicentre_similarity_r = np.zeros((len(thresholds), len(thresholds)))
    sc_epicentre_similarity_p = np.zeros((len(thresholds), len(thresholds)))

    percent_fc_ctx = np.zeros((len(thresholds), 68))
    percent_fc_sctx = np.zeros((len(thresholds), 14))
    percent_sc_ctx = np.zeros((len(thresholds), 68))
    percent_sc_sctx = np.zeros((len(thresholds), 14))

    for i in range(len(thresholds)):
        thr1 = regional_association[thresholds[i]].t

        fc_ctx_r1, fc_ctx_p1 = util.epicenter_mapping(thr1, fc_ctx)
        fc_sctx_r1, fc_sctx_p1 = util.epicenter_mapping(thr1, fc_sctx)
        sc_ctx_r1, sc_ctx_p1 = util.epicenter_mapping(thr1, sc_ctx)
        sc_sctx_r1, sc_sctx_p1 = util.epicenter_mapping(thr1, sc_sctx)

        thr_fc_epicentre[i, :] = np.concatenate((fc_ctx_r1, fc_sctx_r1))
        thr_sc_epicentre[i, :] = np.concatenate((sc_ctx_r1, sc_sctx_r1))

        percent_fc_ctx[i, :] = fc_ctx_p1 < 0.05
        percent_fc_sctx[i, :] = fc_sctx_p1 < 0.05
        percent_sc_ctx[i, :] = sc_ctx_p1 < 0.05
        percent_sc_sctx[i, :] = sc_sctx_p1 < 0.05

        for j in range(i, len(thresholds)):
            thr2 = regional_association[thresholds[j]].t

            fc_ctx_r2, _ = util.epicenter_mapping(thr2, fc_ctx)
            fc_sctx_r2, _ = util.epicenter_mapping(thr2, fc_sctx)
            sc_ctx_r2, _ = util.epicenter_mapping(thr2, sc_ctx)
            sc_sctx_r2, _ = util.epicenter_mapping(thr2, sc_sctx)

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

    percent_fc_ctx = np.mean(percent_fc_ctx, axis=0)
    percent_fc_sctx = np.mean(percent_fc_sctx, axis=0)
    percent_sc_ctx = np.mean(percent_sc_ctx, axis=0)
    percent_sc_sctx = np.mean(percent_sc_sctx, axis=0)

    print()
    print("Save results")
    print("------------------------------------------")
    util.save_to_pickle(
        "../../data/results/s06_epicentreConsistency/percent_epicentre.pkl",
        {
            "percent_fc_ctx": percent_fc_ctx,
            "percent_fc_sctx": percent_fc_sctx,
            "percent_sc_ctx": percent_sc_ctx,
            "percent_sc_sctx": percent_sc_sctx,
        },
    )


if __name__ == "__main__":
    main()
