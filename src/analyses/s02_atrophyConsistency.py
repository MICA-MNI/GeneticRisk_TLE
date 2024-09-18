import os
import pickle
import pandas as pd
import utilities as util


# Helper functions
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


def load_multi_ige():
    """
    Load multisite IGE (ENIGMA) data

    Returns:
    numpy.ndarray: The loaded data containing age, sex, dataset, group, focus, and ct information.
    """
    return util.load_data(
        "../../data/processed/multi_ige_data.npz",
        ["age", "sex", "site", "group", "ct"],
    )


# Main analysis
def main():
    # Load multi TLE data
    age, sex, site, group, focus, ct = load_multi_tle()

    # Load imaging genetic result
    imaging_genetics = util.load_imaging_genetic()

    x = util.zscore_flip(ct, focus, "C", "_")
    covar = pd.DataFrame({"age": age, "sex": sex})
    atrophy, association = {}, {}
    for hemi in ["L", "R"]:
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

    print("-------------------------------")
    print("Idiopathic generalized epilepsy")
    print("-------------------------------")
    # Load multi IGE data
    age, sex, site, group, ct = load_multi_ige()

    x = util.zscore_flip(ct, group, "HC", "_")
    covar = pd.DataFrame({"age": age, "sex": sex})

    print()
    print("Derive atrophy map")
    print("-------------------------------")
    slm = util.casecontrol_difference(x, covar, group, "HC", "IGE")
    atrophy["ige"] = slm

    print()
    print("Correlate with PRS map")
    print("-------------------------------")
    r, p, null = util.spatial_correlation(slm.t, imaging_genetics.t)
    association["ige"] = {"r": r, "p": p, "null": null}

    print(f"Correlation: r = {r}, p = {p}")

    print()
    print("-----------")
    print("Save resuts")
    print("-----------")
    util.save_to_pickle(
        "../../data/results/s02_atrophyConsistency/multi_atrophy.pkl",
        {"atrophy": atrophy},
    )

    util.save_to_pickle(
        "../../data/results/s02_atrophyConsistency/multi_association.pkl",
        {"association": association},
    )


if __name__ == "__main__":
    main()
