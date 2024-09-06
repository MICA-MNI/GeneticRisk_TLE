import pandas as pd
import numpy as np
from neuroCombat import neuroCombat


def main():
    print("--------------------")
    print("Loading ABCD dataset")
    print("--------------------")

    print()
    print("Demographic data")
    print("--------------------")
    abcd_demographics = pd.read_csv(
        "../data/raw/abcd_demographics.csv", index_col="participant"
    )
    age = abcd_demographics["age"]
    sex = abcd_demographics["sex"]
    site = abcd_demographics["site"]
    print(f"# of participants: {len(abcd_demographics)}")
    print(
        f"mean age: {np.mean(age)}; std age: {np.std(age)}; range age: {np.min(age)} - {np.max(age)}"
    )
    print(f'# of males: {np.sum(sex=="M")}; # of F: {np.sum(sex=="F")}')

    print()
    print("Genetic data")
    print("--------------------")
    abcd_genetics = pd.read_csv(
        "../data/raw/abcd_genetics.csv", index_col="participant"
    )
    pc10 = abcd_genetics.filter(like="PC")
    prs_all = abcd_genetics.filter(like="Pt")
    thresh = prs_all.columns.to_list()

    print()
    print("Thickness data")
    print("--------------------")
    # Vertex-wise
    abcd_ct = pd.read_csv("../data/raw/abcd_ct_vertex.csv", index_col="participant")
    ct_vertex = np.transpose(abcd_ct.to_numpy())
    covars = pd.DataFrame({"age": age, "sex": sex == "M", "site": site})
    categorical_cols = ["sex"]
    batch_col = "site"
    ct_vertex = np.transpose(
        neuroCombat(
            dat=ct_vertex,
            covars=covars,
            batch_col=batch_col,
            categorical_cols=categorical_cols,
        )["data"]
    )

    # Parcellated
    abcd_ct = pd.read_csv("../data/raw/abcd_ct_aparc.csv", index_col="participant")
    ct_aparc = np.transpose(abcd_ct.to_numpy())
    covars = pd.DataFrame({"age": age, "sex": sex == "M", "site": site})
    categorical_cols = ["sex"]
    batch_col = "site"
    ct_aparc = np.transpose(
        neuroCombat(
            dat=ct_aparc,
            covars=covars,
            batch_col=batch_col,
            categorical_cols=categorical_cols,
        )["data"]
    )

    print()
    print("Save data")
    print("--------------------")
    np.savez(
        "../data/processed/abcd_data",
        age=age,
        sex=sex,
        site=site,
        pc10=pc10,
        thresh=thresh,
        prs_all=prs_all,
        ct_vertex=ct_vertex,
        ct_aparc=ct_aparc,
    )
    print()

    print()
    print("-------------------------")
    print("Loading local TLE dataset")
    print("-------------------------")

    print()
    print("Demographic data")
    print("-------------------------")
    local_tle_demographics = pd.read_csv(
        "../data/raw/local_tle_demographics.csv", index_col="participant"
    )
    age = local_tle_demographics["age"]
    sex = local_tle_demographics["sex"]
    group = local_tle_demographics["group"]
    focus = local_tle_demographics["focus"]
    dataset = local_tle_demographics["dataset"]
    print(
        f'# of HC: {np.sum(focus=="C")}; # left TLE: {np.sum(focus=="L")}; # of right TLE: {np.sum(focus=="R")}'
    )
    print(
        f'mean age HC: {np.nanmean(age[group=="HC"])}; std age HC: {np.nanstd(age[group=="HC"])}; range age HC: {np.nanmin(age[group=="HC"])} - {np.nanmax(age[group=="HC"])}'
    )
    print(
        f'mean age TLE: {np.nanmean(age[group=="TLE"])}; std age TLE: {np.nanstd(age[group=="TLE"])}; range age TLE: {np.nanmin(age[group=="TLE"])} - {np.nanmax(age[group=="TLE"])}'
    )
    print(
        f'# of M HC: {np.sum((sex=="M") & (group=="HC"))}; # of F HC: {np.sum((sex=="F") & (group=="HC"))}'
    )
    print(
        f'# of M TLE: {np.sum((sex=="M") & (group=="TLE"))}; # of F TLE: {np.sum((sex=="F") & (group=="TLE"))}'
    )

    print()
    print("Thickness data")
    print("-------------------------")
    local_tle_ct = pd.read_csv("../data/raw/local_tle_ct.csv", index_col="participant")
    ct = np.transpose(local_tle_ct.to_numpy())
    covars = pd.DataFrame(
        {
            "age": age,
            "sex": sex == "M",
            "group": group == "TLE",
            "dataset": dataset,
        }
    )
    categorical_cols = ["sex", "group"]
    batch_col = "dataset"
    ct = neuroCombat(
        dat=ct, covars=covars, batch_col=batch_col, categorical_cols=categorical_cols
    )["data"]

    print()
    print("Save data")
    print("-------------------------")
    np.savez(
        "../data/processed/local_tle_data",
        age=age,
        sex=sex,
        dataset=dataset,
        group=group,
        focus=focus,
        ct=ct,
    )

    print()
    print("-------------------------")
    print("Loading multi TLE dataset")
    print("-------------------------")

    print()
    print("Demographic data")
    print("-------------------------")
    multi_tle_demographics = pd.read_csv(
        "../data/raw/multi_tle_demographics.csv", index_col="participant"
    )
    age = multi_tle_demographics["age"]
    sex = multi_tle_demographics["sex"]
    group = multi_tle_demographics["group"]
    focus = multi_tle_demographics["focus"]
    site = multi_tle_demographics["site"]
    print(
        f"# of HC: {np.sum(focus=='C')}; # left TLE: {np.sum(focus=='L')}; # of right TLE: {np.sum(focus=='R')}"
    )
    print(
        f"mean age HC: {np.nanmean(age[group=='HC'])}; std age HC: {np.nanstd(age[group=='HC'])}; range age HC: {np.nanmin(age[group=='HC'])} - {np.nanmax(age[group=='HC'])}"
    )
    print(
        f"mean age TLE: {np.nanmean(age[group=='TLE'])}; std age TLE: {np.nanstd(age[group=='TLE'])}; range age TLE: {np.nanmin(age[group=='TLE'])} - {np.nanmax(age[group=='TLE'])}"
    )
    print(
        f"# of M HC: {np.sum((sex=='M') & (group=='HC'))}; # of F HC: {np.sum((sex=='F') & (group=='HC'))}"
    )
    print(
        f"# of M TLE: {np.sum((sex=='M') & (group=='TLE'))}; # of F TLE: {np.sum((sex=='F') & (group=='TLE'))}"
    )

    print()
    print("Thickness data")
    print("-------------------------")
    multi_tle_ct = pd.read_csv("../data/raw/multi_tle_ct.csv", index_col="participant")
    ct = np.transpose(multi_tle_ct.to_numpy())

    print()
    print("Save data")
    print("-------------------------")
    np.savez(
        "../data/processed/multi_tle_data",
        age=age,
        sex=sex,
        dataset=dataset,
        group=group,
        focus=focus,
        ct=ct,
    )


if __name__ == "__main__":
    main()
