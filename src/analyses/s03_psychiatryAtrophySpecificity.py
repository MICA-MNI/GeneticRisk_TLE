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
        "../../data/results/s03_psychiatryAtrophySpecificity/psychiatric_atrophy.pkl",
        {"atrophy": atrophy},
    )

    util.save_to_pickle(
        "../../data/results/s03_psychiatryAtrophySpecificity/psychiatric_association.pkl",
        {"association": association},
    )


if __name__ == "__main__":
    main()
