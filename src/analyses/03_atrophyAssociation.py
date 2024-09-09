import numpy as np
import pandas as pd
import pickle
import utilities as util
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM


# Helper functions
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
    # Load LTE data
    age, sex, dataset, focus, ct = load_local_tle()

    # Load imaging genetic result
    imaging_genetics = util.load_imaging_genetic()
    x = zscore_flip(ct, focus, "C", "R")
    covar = pd.concat([age, sex])
    atrophy, atrophy_association = {}, {}
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

        atrophy_association[f"{hemi.lower()}tle"] = {"r": r, "p": p, "null": null}

    print()
    print("------------")
    print("Save results")
    print("------------")
    np.savez(
        "../../data/results/03_atrophyAssociation/atrophy_association.npz",
        atrophy=atrophy,
        atrophy_association=atrophy_association,
    )


if __name__ == "__main__":
    main()
