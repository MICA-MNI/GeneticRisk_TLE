import numpy as np
from 03_casecontrolAssociation import casecontrol_difference, spatial_correlation

def main()

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

        print()
        print('Correlate with PRS map')
        print('------------------------')
        r, p, null = spatial_correlation(slm.t, imaging_genetics.t)

        epilepsy_atrophy[f'{hemi.lower()}tle'] = {'slm':slm, 'r':r, 'p':p, 'null':null}

    x = zscore_flip(ct, covar, focus, 'C', 'X')
    covar = pd.concat([age, sex])
    print('-------------------------------')
    print('Idiopathic generalized epilepsy')
    print('-------------------------------')

    print()
    print('Derive atrophy map')
    print('-------------------------------')
    slm = casecontrol_difference(x, covar, group, 'HC', 'IGE')

    print()
    print('Correlate with PRS map')
    print('-------------------------------')
    r, p, null = spatial_correlation(slm.t, imaging_genetics.t)

    epilepsy_atrophy['ige'] = {'slm':slm, 'r':r, 'p':p, 'null'null}

    print('-----------')
    print('Save resuts')
    print('-----------')


if __name__ == "__main__":
    main()
