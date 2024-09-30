import numpy as np
import utilities as util


# Main analysis
def main():
    atrophy = util.load_result(
        "../../data/results/03_atrophyAssociation/local_atrophy.pkl", ["atrophy"]
    )[0]

    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()
    prs_fc_ctx, prs_fc_sctx, prs_sc_ctx, prs_sc_sctx = util.load_imaging_genetic(
        "network"
    )

    epicentre = {}
    epicentre_association = {}
    for hemi in ["L", "R"]:
        print("------------------------")
        print(f"{hemi} temporal lobe epilepsy")
        print("------------------------")

        print()
        print("Derive disease epicentre map")
        print("------------------------")
        atrophy_map = atrophy[f"{hemi.lower()}tle"].t
        fc_ctx_r, fc_ctx_p = util.epicenter_mapping(atrophy_map, fc_ctx)
        fc_sctx_r, fc_sctx_p = util.epicenter_mapping(atrophy_map, fc_sctx)
        sc_ctx_r, sc_ctx_p = util.epicenter_mapping(atrophy_map, sc_ctx)
        sc_sctx_r, sc_sctx_p = util.epicenter_mapping(atrophy_map, sc_sctx)
        epicentre[f"{hemi.lower()}tle"] = {
            "fc_ctx": {"r": fc_ctx_r, "p": fc_ctx_p},
            "fc_sctx": {"r": fc_ctx_r, "p": fc_sctx_p},
            "sc_ctx": {"r": sc_ctx_r, "p": sc_ctx_p},
            "sc_sctx": {"r": sc_sctx_r, "p": sc_sctx_p},
        }

        print()
        print("Correlate with PRS epicentre map")
        print("------------------------")
        r_fc_epi, p_fc_epi, null_fc_epi = util.spatial_correlation(
            np.concatenate((fc_ctx_r, fc_sctx_r)),
            np.concatenate((prs_fc_ctx, prs_fc_sctx)),
            surface_name="fsa5_with_sctx",
            parcellation_name="aparc_aseg",
        )
        r_sc_epi, p_sc_epi, null_sc_epi = util.spatial_correlation(
            np.concatenate((sc_ctx_r, sc_sctx_r)),
            np.concatenate((prs_sc_ctx, prs_sc_sctx)),
            surface_name="fsa5_with_sctx",
            parcellation_name="aparc_aseg",
        )

        epicentre_association[f"{hemi.lower()}tle"] = {
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
    print("------------")
    print("Save results")
    print("------------")
    np.savez(
        "../../data/results/04_epicentreAssociation/epicentre_association.npz",
        epicentre=epicentre,
        epicentre_association=epicentre_association,
    )


if __name__ == "__main__":
    main()
