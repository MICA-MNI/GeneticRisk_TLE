import numpy as np
import utilities as util


# Main analysis
def main():
    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()
    prs_fc_ctx, prs_fc_sctx, prs_sc_ctx, prs_sc_sctx = util.load_imaging_genetic(
        "network"
    )

    for site in ["multi"]: #["local", "multi"]:
        print("------------------------")
        print(f"Dataset: {site}")
        print("------------------------")

        # Load atrophy maps
        atrophy = util.load_result(
            f"../../data/results/03_atrophyAssociation/{site}_atrophy.pkl", ["atrophy"]
        )[0]

        epicentre = {}
        association = {}
        for hemi in ["L", "R"]:
            print()
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

            association[f"{hemi.lower()}tle"] = {
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

            print(f"FC epicentre correlation: r = {r_fc_epi}, p = {p_fc_epi}")
            print(f"SC epicentre correlation: r = {r_sc_epi}, p = {p_sc_epi}")

        print()
        print("------------")
        print("Save results")
        print("------------")
        util.save_to_pickle(
            f"../../data/results/04_epicentreAssociation/{site}_epicentre.pkl",
            {"epicentre": epicentre},
        )

        util.save_to_pickle(
            f"../../data/results/04_epicentreAssociation/{site}_association.pkl",
            {"association": association},
        )


if __name__ == "__main__":
    main()
