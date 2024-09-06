import numpy as np
import pandas as pd
import pickle
from 02_epicentreMapping import spatial_correlation, epicentre_mapping


# Main analysis
def main():

    atrophy = np.load('../../data/results/03_atrophyAssociation/atrophy_association.npz')
    atrophy = atrophy['atrophy']

    fc_ctx, _, fc_sctx, _ = load_fc()
    sc_ctx, _, sc_sctx, _ = load_sc()

    for hemi in ['L', 'R']:
        print('------------------------')
        print(f'{hemi} temporal lobe epilepsy')
        print('------------------------')

        print()
        print('Derive disease epicentre map')
        print('------------------------')
        atrophy_map = atrophy[f'{hemi.lower()}tle'].t
        fc_ctx_r, fc_ctx_p = epicenter_mapping(atrophy_map, fc_ctx)
        fc_sctx_r, fc_sctx_p = epicenter_mapping(atrophy_map, fc_sctx)
        sc_ctx_r, sc_ctx_p = epicenter_mapping(atrophy_map, sc_ctx)
        sc_sctx_r, sc_sctx_p = epicenter_mapping(atrophy_map, sc_sctx)
        epicentre[f'{hemi.lower()}tle'] = {
            'fc_ctx': {'r':fc_ctx_r, 'p':fc_ctx_p},
            'fc_sctx': {'r':fc_ctx_r, 'p':fc_ctx_p},
            'sc_ctx': {'r':sc_ctx_r, 'p':sc_ctx_p},
            'sc_sctx': {'r':sc_sctx_r, 'p':sc_sctx_p}
        }

        print()
        print('Correlate with PRS epicentre map')
        print('------------------------')
        r_fc_epi, p_fc_epi, null_fc_epi = spatial_correlation(np.concatenate(fc_ctx_r, fc_sctx_r), imaging_genetic_epicentre, parcellation_name='aparc_aseg', n_rot=5000)
        r_sc_epi, p_sc_epi, null_sc_epi = spatial_correlation(np.concatenate(sc_ctx_r, sc_sctx_r), imaging_genetic_epicentre, parcellation_name='aparc_aseg', n_rot=5000)
        epicentre_association[f'{hemi.lower()}tle'] = {
            'fc_epi': {'r_fc_epi':r_fc_epi, 'p_fc_epi':p_fc_epi, 'null_fc_epi':null_fc_epi},
            'sc_epi': {'r_sc_epi':r_sc_epi, 'p_sc_epi':p_sc_epi, 'null_sc_epi':null_sc_epi}
        }


    print()
    print('------------')
    print('Save results')
    print('------------')
    np.savez('../../data/results/04_epicentreAssociation/epicentre_association.npz',
            epicentre=epicentre, epicentre_association=epicentre_association)


if __name__ == "__main__":
    main()
