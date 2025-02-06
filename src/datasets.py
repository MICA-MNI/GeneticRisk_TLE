import glob
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
    print("Morphological data")
    print("--------------------")
    # Vertex-wise cortical thickness
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

    # Parcellated corticla thickness
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

    # Subcortical/intracranial volume
    abcd_sv = pd.read_csv("../data/raw/abcd_sv.csv", index_col="participant")
    sv = np.transpose(abcd_sv.to_numpy())
    covars = pd.DataFrame({"age": age, "sex": sex == "M", "site": site})
    categorical_cols = ["sex"]
    batch_col = "site"
    sv = np.transpose(
        neuroCombat(
            dat=sv,
            covars=covars,
            batch_col=batch_col,
            categorical_cols=categorical_cols,
        )["data"]
    )
    icv, sv = sv[:, sv.shape[1] - 1], sv[:, : sv.shape[1] - 1]

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
        sv=sv,
        icv=icv,
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
    ct = np.transpose(
        neuroCombat(
            dat=ct,
            covars=covars,
            batch_col=batch_col,
            categorical_cols=categorical_cols,
        )["data"]
    )

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
    ct = multi_tle_ct.to_numpy()

    print()
    print("Save data")
    print("-------------------------")
    np.savez(
        "../data/processed/multi_tle_data",
        age=age,
        sex=sex,
        site=site,
        group=group,
        focus=focus,
        ct=ct,
    )

    print()
    print("-------------------------")
    print("Loading multi IGE dataset")
    print("-------------------------")

    print()
    print("Demographic data")
    print("-------------------------")
    multi_ige_demographics = pd.read_csv(
        "../data/raw/multi_ige_demographics.csv", index_col="participant"
    )
    age = multi_ige_demographics["age"]
    sex = multi_ige_demographics["sex"]
    group = multi_ige_demographics["group"]
    site = multi_ige_demographics["site"]
    print(f"# of HC: {np.sum(group=='HC')}; # IGE: {np.sum(group=='IGE')}")
    print(
        f"mean age HC: {np.nanmean(age[group=='HC'])}; std age HC: {np.nanstd(age[group=='HC'])}; range age HC: {np.nanmin(age[group=='HC'])} - {np.nanmax(age[group=='HC'])}"
    )
    print(
        f"mean age IGE: {np.nanmean(age[group=='IGE'])}; std age IGE: {np.nanstd(age[group=='IGE'])}; range age IGE: {np.nanmin(age[group=='IGE'])} - {np.nanmax(age[group=='IGE'])}"
    )
    print(
        f"# of M HC: {np.sum((sex=='M') & (group=='HC'))}; # of F HC: {np.sum((sex=='F') & (group=='HC'))}"
    )
    print(
        f"# of M IGE: {np.sum((sex=='M') & (group=='IGE'))}; # of F IGE: {np.sum((sex=='F') & (group=='IGE'))}"
    )

    print()
    print("Thickness data")
    print("-------------------------")
    multi_ige_ct = pd.read_csv("../data/raw/multi_ige_ct.csv", index_col="participant")
    ct = multi_ige_ct.to_numpy()

    print()
    print("Save data")
    print("-------------------------")
    np.savez(
        "../data/processed/multi_ige_data",
        age=age,
        sex=sex,
        site=site,
        group=group,
        ct=ct,
    )

    print()
    print("-----------------------")
    print("Loading HCP connectomes")
    print("-----------------------")

    print()
    print("Functional connectomes")
    print("-----------------------")
    # Load files
    files = glob.glob("../data/raw/connectomes/*FC*")

    # Compute average connectome
    for i, f in enumerate(files):
        mtx = np.loadtxt(f, delimiter=",")
        mtx_z = np.arctanh(mtx)
        mtx_z[(mtx_z < 0) | np.isinf(mtx_z)] = 0
        mtx_z = np.triu(mtx_z, 1) + mtx_z.T

        if i == 0:
            fc = np.empty((mtx_z.shape[0], mtx_z.shape[1], len(files)))

        fc[:, :, i] = mtx_z

    mean_fc = np.mean(fc, axis=2)

    # Isolate cortico-cortical and cortico-subcortical functional connectivity
    ctx_idx = np.concatenate(
        [np.arange(49, 52), np.arange(53, 84), np.arange(84, 87), np.arange(88, 119)]
    )
    sctx_idx = np.array([6, 5, 1, 4, 3, 2, 0, 13, 12, 8, 11, 10, 9, 7])
    fc_ctx = mean_fc[np.ix_(ctx_idx, ctx_idx)]
    fc_sctx = mean_fc[np.ix_(sctx_idx, ctx_idx)]

    print()
    print("Structural connectomes")
    print("-----------------------")
    # Load files
    files = glob.glob("../data/raw/connectomes/*SC*")

    # Compute average connectome
    for i, f in enumerate(files):
        mtx = np.loadtxt(f, delimiter=",")
        mtx_log = np.log(np.triu(mtx, 1) + mtx.T)
        mtx_log[np.isinf(mtx_log)] = 0

        if i == 0:
            sc = np.empty((mtx_log.shape[0], mtx_log.shape[1], len(files)))

        sc[:, :, i] = mtx_log

    mean_sc = np.mean(sc, axis=2)

    # Isolate cortico-cortical and cortico-subcortical structural connectivity
    ctx_idx = np.concatenate(
        [np.arange(49, 52), np.arange(53, 84), np.arange(85, 88), np.arange(89, 120)]
    )
    sctx_idx = np.array([6, 5, 1, 4, 3, 2, 0, 13, 12, 8, 11, 10, 9, 7])
    sc_ctx = mean_sc[np.ix_(ctx_idx, ctx_idx)]
    sc_sctx = mean_sc[np.ix_(sctx_idx, ctx_idx)]

    print()
    print("Save data")
    print("-----------------------")
    np.savez(
        "../data/processed/hcp_data",
        fc_ctx=fc_ctx,
        fc_sctx=fc_sctx,
        sc_ctx=sc_ctx,
        sc_sctx=sc_sctx,
    )


if __name__ == "__main__":
    main()
