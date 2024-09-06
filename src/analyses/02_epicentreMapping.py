import numpy as np
from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.permutation_testing import spin_test


# Helper functions
def spatial_correlation(map1, map2, n_rot=5000, surface_name='fsa5', parcellation_name='aparc'):
    """
    Calculate the spatial correlation between two maps and perform a spin permuation spin_test

    Parameters:
    map1 (array-like): First spatial map
    map2 (array-like): Second spatial map
    n_rot (int, optional): Number of rotations. Default is 5000
    surface_name (str, optional): Name of surface. Default is 'fsa5'
    parcellation_name (str, optional): Name of parcellation. Default is 'aparc'

    Returns:
    r (float): Pearson's correlation coefficient between map1 and map2
    p (float): Permutation-based p-pvalue
    null (array-like): Null distribution from the permuation spin_test
    """
    r = np.corrcoef(map1,map2)[0,1]
    p, null = spin_test(map1, map2, surface_name=surface_name, parcellation_name=parcellation_name, n_rot=n_rot, null_distr=True)

    return r, p, null

def epicenter_mapping(map, connectome):
    """
    Map epicentres of a given map to a connectome

    Parameters:
    map (array-like): Spatial map
    connectome (array-like): Connectome matrix

    Returns:
    epi_r: Correlation coefficients of each seed regions
    epi_p: Permuation-based p-values of each seed regions
    """
    epi_r = []
    epi_p = []
    for seed in range(connectome.shape[0]):
        seed_con = connectome[:, seed]
        r, p, _ = spatial_correlation(seed_con, map)
        epi_r = np.append(epi_r, r)
        epi_p = np.append(epi_p, p)

    return epi_r, epi_p


# Main analysis
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
