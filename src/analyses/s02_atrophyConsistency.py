import numpy as np
import pandas as pd
from 02_epicentreMapping import spatial_correlation
from 03_casecontrolAssociation import casecontrol_difference, zscore_flip


# Main analysis
def main():

    x = zscore_flip(ct, covar, focus, 'C', 'R')
    covar = pd.concat([age, sex])
    for hemi in ['L', 'R']:
        print('------------------------')
        print(f'{hemi} temporal lobe epilepsy')
        print('------------------------')

        print()
        print('Derive atrophy map')
        print('------------------------')
        slm = casecontrol_difference(x, covar, focus, 'C', hemi)
        atrophy[f'{hemi.lower()}tle'] = slm

        print()
        print('Correlate with PRS map')
        print('------------------------')
        r, p, null = spatial_correlation(slm.t, imaging_genetics.t)

        association[f'{hemi.lower()}tle'] = {'r':r, 'p':p, 'null':null}


    x = zscore_flip(ct, covar, focus, 'C', 'X')
    covar = pd.concat([age, sex])
    print('-------------------------------')
    print('Idiopathic generalized epilepsy')
    print('-------------------------------')

    print()
    print('Derive atrophy map')
    print('-------------------------------')
    slm = casecontrol_difference(x, covar, group, 'HC', 'IGE')
    atrophy['ige'] = slm

    print()
    print('Correlate with PRS map')
    print('-------------------------------')
    r, p, null = spatial_correlation(slm.t, imaging_genetics.t)
    association['ige'] = {'r':r, 'p':p, 'null':null}

    print()
    print('-----------')
    print('Save resuts')
    print('-----------')
    np.savez('../../data/results/s02_atrophyConsistency/atrophy_association.npz',
            atrophy=atrophy, association=association)


if __name__ == "__main__":
    main()
