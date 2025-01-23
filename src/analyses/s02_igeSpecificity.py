import numpy as np
import pandas as pd
import utilities as util


# Helper functions
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
    print()
    print("--------------------")
    print("Case-control atrophy")
    print("--------------------")

    age, sex, group, ct = load_multi_ige()

    # Load imaging genetic result
    imaging_genetics = util.load_imaging_genetic("regional")

    atrophy, atrophy_association = {}, {}

    print()
    print("Derive atrophy map")
    print("--------------------")
    x = util.zscore_flip(ct, group, "HC", "_")
    covar = pd.DataFrame({"age": age, "sex": sex})

    slm = util.casecontrol_difference(x, covar, group, "C", "X")
    atrophy["ige"] = slm

    print()
    print("Correlate with PRS map")
    print("--------------------")
    r, p, null = util.spatial_correlation(slm.t, imaging_genetics.t)

    atrophy_association["ige"] = {"r": r, "p": p, "null": null}

    print(f"Correlation: r = {r}, p = {p}")

    print()
    print("Save results")
    print("--------------------")
    util.save_to_pickle(
        "../../data/results/s02_igeSpecificity/ige_atrophy.pkl",
        {"atrophy": atrophy},
    )

    util.save_to_pickle(
        "../../data/results/s02_igeSpecificity/ige_atrophyAssociation.pkl",
        {"association": atrophy_association},
    )

    print()
    print("---------------------------")
    print("Network epilepsy epicentres")
    print("---------------------------")
    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()

    prs_fc_ctx, prs_fc_sctx, prs_sc_ctx, prs_sc_sctx = util.load_imaging_genetic(
        "network"
    )
    prs_fc_epi = np.concatenate((prs_fc_ctx, prs_fc_sctx))
    prs_sc_epi = np.concatenate((prs_sc_ctx, prs_sc_sctx))

    epicentre, network_association = {}, {}

    print()
    print("Derive disease epicentre")
    print("---------------------------")
    atrophy_map = atrophy["ige"].t
    fc_ctx_r, fc_ctx_p = util.epicenter_mapping(atrophy_map, fc_ctx)
    fc_sctx_r, fc_sctx_p = util.epicenter_mapping(atrophy_map, fc_sctx)
    sc_ctx_r, sc_ctx_p = util.epicenter_mapping(atrophy_map, sc_ctx)
    sc_sctx_r, sc_sctx_p = util.epicenter_mapping(atrophy_map, sc_sctx)
    epicentre["ige"] = {
        "fc_ctx": {"r": fc_ctx_r, "p": fc_ctx_p},
        "fc_sctx": {"r": fc_ctx_r, "p": fc_ctx_p},
        "sc_ctx": {"r": sc_ctx_r, "p": sc_ctx_p},
        "sc_sctx": {"r": sc_sctx_r, "p": sc_sctx_p},
    }

    print()
    print("Correlate with PRS epicentres")
    print("---------------------------")
    r_fc_epi, p_fc_epi, null_fc_epi = util.spatial_correlation(
        np.concatenate((fc_ctx_r, fc_sctx_r)),
        prs_fc_epi,
        surface_name="fsa5_with_sctx",
        parcellation_name="aparc_aseg",
        n_rot=5000,
    )

    r_sc_epi, p_sc_epi, null_sc_epi = util.spatial_correlation(
        np.concatenate((sc_ctx_r, sc_sctx_r)),
        prs_sc_epi,
        surface_name="fsa5_with_sctx",
        parcellation_name="aparc_aseg",
        n_rot=5000,
    )

    network_association["ige"] = {
        "fc_epi": {
            "r_fc_epi": r_fc_epi,
            "p_fc_epi": p_fc_epi,
            "null_fc_epi": null_fc_epi,
        },
        "sc_epi": {
            "r_sc_epi": r_sc_epi,
            "p_sc_epi": p_sc_epi,
            "null_sc_epi": null_sc_epi,
        },
    }

    print()
    print("Save results")
    print("---------------------------")
    util.save_to_pickle(
        "../../data/results/s02_igeSpecificity/ige_epicentre.pkl",
        {"epicentre": epicentre},
    )

    util.save_to_pickle(
        "../../data/results/s02_igeSpecificity/ige_epicentreAssociation.pkl",
        {"association": network_association},
    )


if __name__ == "__main__":
    main()
