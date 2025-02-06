import numpy as np
import utilities as util


# Main analysis
def main():
    # Load atrophy from ENIGMA
    atrophy = util.load_result(
        "../../data/results/s03_psychiatryAtrophySpecificity/psychiatric_atrophy.pkl",
        ["atrophy"],
    )[0]

    print("------------------------------")
    print("Network psychiatric epicentres")
    print("------------------------------")
    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()

    prs_fc_ctx, prs_fc_sctx, prs_sc_ctx, prs_sc_sctx = util.load_imaging_genetic(
        "network"
    )
    prs_fc_epi = np.concatenate((prs_fc_ctx, prs_fc_sctx))
    prs_sc_epi = np.concatenate((prs_sc_ctx, prs_sc_sctx))

    epicentre, association = {}, {}
    for psy in atrophy:
        print("------------------------------")
        print(psy.upper())
        print("------------------------------")

        print()
        print("Derive epicentre map")
        print("------------------------------")
        atrophy_map = atrophy[psy]
        fc_ctx_r, fc_ctx_p = util.epicenter_mapping(atrophy_map, fc_ctx)
        fc_sctx_r, fc_sctx_p = util.epicenter_mapping(atrophy_map, fc_sctx)
        sc_ctx_r, sc_ctx_p = util.epicenter_mapping(atrophy_map, sc_ctx)
        sc_sctx_r, sc_sctx_p = util.epicenter_mapping(atrophy_map, sc_sctx)
        epicentre[psy] = {
            "fc_ctx": {"r": fc_ctx_r, "p": fc_ctx_p},
            "fc_sctx": {"r": fc_ctx_r, "p": fc_sctx_p},
            "sc_ctx": {"r": sc_ctx_r, "p": sc_ctx_p},
            "sc_sctx": {"r": sc_sctx_r, "p": sc_sctx_p},
        }

        print()
        print("Correlate with PRS epicentres")
        print("------------------------------")
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

        association[psy] = {
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

        print(f"FC epicentre: {r_fc_epi} ({p_fc_epi})")
        print(f"SC epicentre: {r_sc_epi} ({p_sc_epi})")

    print()
    print("------------")
    print("Save results")
    print("------------")
    util.save_to_pickle(
        "../../data/results/s04_psychiatryEpicentreSpecificity/psychiatric_epicentre.pkl",
        {"epicentre": epicentre},
    )

    util.save_to_pickle(
        "../../data/results/s04_psychiatryEpicentreSpecificity/psychiatric_association.pkl",
        {"association": association},
    )


if __name__ == "__main__":
    main()
