import numpy as np
import utilities as util
from enigmatoolbox.datasets import load_summary_stats


# Helper function
def load_psychiatry():
    """
    Loads atrophy data for various psychiatric conditions:
    - 'adhd': Attention Deficit Hyperactivity Disorder
    - 'asd': Autism Spectrum Disorder
    - 'bd': Bipolar Disorder
    - 'mdd': Major Depressive Disorder
    - 'ocd': Obsessive-Compulsive Disorder
    - 'scz': Schizophrenia

    Returns:
    dict: A dictionary where each key is a psychiatric condition and each value is
          the cohen's D effect size from the summary statistics for that condition.
    """
    atrophy = {
        "adhd": load_summary_stats("adhd")["CortThick_case_vs_controls_adult"]["d_icv"],
        "asd": load_summary_stats("asd")["CortThick_case_vs_controls_meta_analysis"][
            "d_icv"
        ],
        "bd": load_summary_stats("bipolar")["CortThick_case_vs_controls_adult"][
            "d_icv"
        ],
        "mdd": load_summary_stats("depression")["CortThick_case_vs_controls_adult"][
            "d_icv"
        ],
        "ocd": load_summary_stats("ocd")["CortThick_case_vs_controls_adult"]["d_icv"],
        "scz": load_summary_stats("schizophrenia")["CortThick_case_vs_controls"][
            "d_icv"
        ],
    }

    return atrophy


# Main analysis
def main():
    # Load atrophy from ENIGMA
    atrophy = load_psychiatry()

    # Load imaging genetic result
    imaging_genetics = util.load_imaging_genetic("regional")

    association = {}
    print("----------------------")
    print("Correlate with PRS map")
    print("----------------------")
    for psy in atrophy:
        print()
        print(f"{psy}")
        print("----------------------")
        r, p, null = util.spatial_correlation(atrophy[psy], imaging_genetics.t)

        association[psy] = {"r": r, "p": p, "null": null}

        print(f"Correlation: r = {r}, p = {p}")

    print()
    print("------------")
    print("Save results")
    print("------------")
    util.save_to_pickle(
        "../../data/results/s06_psychiatryAtrophySpecificity/psychiatric_atrophy.pkl",
        {"atrophy": atrophy},
    )

    util.save_to_pickle(
        "../../data/results/s06_psychiatryAtrophySpecificity/psychiatric_association.pkl",
        {"association": association},
    )

    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()

    prs_fc_ctx, prs_fc_sctx, prs_sc_ctx, prs_sc_sctx = util.load_imaging_genetic(
        "network"
    )
    prs_fc_epi = np.concatenate((prs_fc_ctx, prs_fc_sctx))
    prs_sc_epi = np.concatenate((prs_sc_ctx, prs_sc_sctx))

    epicentre, association = {}, {}
    for hemi in ["L", "R"]:
        print("-------------------------")
        print(f"{hemi} temporal lobe epilepsy")
        print("-------------------------")

        print()
        print("Derive epicentre map")
        print("-------------------------")
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
        print("Correlate with PRS epicentres")
        print("-------------------------")
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

    print("--------------------------------")
    print(f"Idiopathic generalized epilepsy")
    print("--------------------------------")

    print()
    print("Derive disease epicentre")
    print("--------------------------------")
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
    print("--------------------------------")
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

    association["ige"] = {
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
    util.save_to_pickle(
        "../../data/results/s05_epicentreAssociation/multi_epicentre.pkl",
        {"epicentre": epicentre},
    )

    util.save_to_pickle(
        "../../data/results/s05_epicentreAssociation/multi_association.pkl",
        {"association": association},
    )
