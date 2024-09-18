import utilities as util


# Main analysis
def main():
    print("---------------------------")
    print("Polygenic epicenter mapping")
    print("---------------------------")

    fc_ctx, fc_sctx, sc_ctx, sc_sctx = util.load_connectomes()
    prs_map = util.load_imaging_genetic("regional")

    print()
    print("Functional networks")
    print("---------------------------")
    # Cortical epicenters
    fc_ctx_r, fc_ctx_p = util.epicenter_mapping(prs_map, fc_ctx)
    # Subcortical epicenters
    fc_sctx_r, fc_sctx_p = util.epicenter_mapping(prs_map, fc_sctx)

    print()
    print("Structural networks")
    print("---------------------------")
    # Cortical epicenters
    sc_ctx_r, sc_ctx_p = util.epicenter_mapping(prs_map, sc_ctx)
    # Subcortical epicenters
    sc_sctx_r, sc_sctx_p = util.epicenter_mapping(prs_map, sc_sctx)

    print()
    print("Save results")
    print("---------------------------")
    util.save_to_pickle(
        "../../data/results/02_epicentreMapping/epicentre.pkl",
        {
            "fc_ctx_r": fc_ctx_r,
            "fc_ctx_p": fc_ctx_p,
            "fc_sctx_r": fc_sctx_r,
            "fc_sctx_p": fc_sctx_p,
            "sc_ctx_r": sc_ctx_r,
            "sc_ctx_p": sc_ctx_p,
            "sc_sctx_r": sc_sctx_r,
            "sc_sctx_p": sc_sctx_p,
        },
    )


if __name__ == "__main__":
    main()
