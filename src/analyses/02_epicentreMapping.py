import numpy as np
from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.permutation_testing import spin_test


# Helper function
def epicenter_mapping(map, connectome):
    epi_r = []
    epi_p = []
    for seed in range(connectome.shape[0]):
        seed_con = connectome[:, seed]
        epi_r = np.append(epi_r, np.corrcoef(seed_con, map)[0, 1])
        epi_p = np.append(epi_p, spin_test(seed_con, map, surface_name='fsa5', parcellation_name='aparc', type='pearson', n_rot=5000, null_dist=False))

    return epi_r, epi_p


def main():

    print('---------------------------')
    print('Polygenic epicenter mapping')
    print('---------------------------')

    fc_ctx, _, fc_sctx, _ = load_fc()
    sc_ctx, _, sc_sctx, _ = load_sc()
    prs_map =

    print()
    print('Functional networks')
    print('---------------------------')
    # Cortical epicenters
    fc_ctx_r, fc_ctx_p = epicenter_mapping(prs_map, fc_ctx)
    # Subcortical epicenters
    fc_sctx_r, fc_sctx_p = epicenter_mapping(prs_map, fc_sctx)

    print()
    print('Structural networks')
    print('---------------------------')
    # Cortical epicenters
    sc_ctx_r, sc_ctx_p = epicenter_mapping(prs_map, sc_ctx)
    # Subcortical epicenters
    sc_sctx_r, sc_sctx_p = epicenter_mapping(prs_map, sc_sctx)

    print()
    print('Save results')
    print('---------------------------')
    np.savez('../../data/results/02_epicentreMapping/epicentre.npz',
             fc_ctx_r=fc_ctx_r,  fc_ctx_p=fc_ctx_p,
             fc_sctx_r=fc_sctx_r, fc_sctx_p=fc_sctx_p,
             sc_ctx_r=sc_ctx_r, sc_ctx_p=sc_ctx_p,
             sc_sctx_r=sc_sctx_r, sc_sctx_p=sc_sctx_p)


if __name__ == "__main__":
    main()
