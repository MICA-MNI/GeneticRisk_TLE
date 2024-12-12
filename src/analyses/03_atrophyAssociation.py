import numpy as np
import pandas as pd
import os
import pickle
import utilities as util


# Helper functions
def load_local_tle():
    """
    Load local TLE (MICs x NKG) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/local_tle_data.npz",
        ["age", "sex", "focus", "ct"],
    )


def load_multi_tle():
    """
    Load multisite TLE (ENIGMA) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, group, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/multi_tle_data.npz",
        ["age", "sex", "site", "group", "focus", "ct"],
    )


# Main analysis
def main():
    load_dataset = {"local": load_local_tle, "multi": load_multi_tle}
    for key in load_dataset:
        # Load LTE data
        age, sex, focus, ct = load_dataset[key]()

        # Load imaging genetic result
        imaging_genetics = util.load_imaging_genetic()

        x = util.zscore_flip(ct, focus, "C", "_")
        covar = pd.DataFrame({"age": age, "sex": sex})
        atrophy, association = {}, {}
        for hemi in ["L", "R"]:
            print()
            print("------------------------")
            print(f"{hemi} temporal lobe epilepsy")
            print("------------------------")

            print()
            print("Derive atrophy map")
            print("------------------------")
            slm = util.casecontrol_difference(x, covar, focus, "C", hemi)
            atrophy[f"{hemi.lower()}tle"] = slm

            print()
            print("Correlate with PRS map")
            print("------------------------")
            r, p, null = util.spatial_correlation(slm.t, imaging_genetics.t)

            association[f"{hemi.lower()}tle"] = {"r": r, "p": p, "null": null}

            print(f"Correlation: r = {r}, p = {p}")

        print()
        print("------------")
        print("Save results")
        print("------------")
        util.save_to_pickle(
            f"../../data/results/03_atrophyAssociation/{key}_atrophy.pkl",
            {"atrophy": atrophy},
        )

        util.save_to_pickle(
            f"/../../data/results/03_atrophyAssociation/{key}_association.pkl",
            {"association": association},
        )


if __name__ == "__main__":
    main()
