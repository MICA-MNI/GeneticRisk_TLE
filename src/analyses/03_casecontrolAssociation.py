import numpy as np
from enigmatoolbox.permutation_testing import spin_test

# Helper functions
def zscore_flip(data, group, control, flip):
    # z-score normalization
    hc = group==control
    hc_mean = np.mean(data[hc,:], axis=0)
    hc_std = np.mean(data[hc,:], axis=0)
    z = (data - hc_mean) / hc_std

    # Sort based on ipsilateral/contralateral
    rtle = group==flip
    n_region = data.shape[1]/2
    f = z.copy()
    f[rtle, :n_region] = z[rtle, n_region:]
    f[rtle, n_region:] = z[rtle, :n_region]

    return f

def casecontrol_difference(data, covar, group, control, patient):
    model = FixedEffect(group)
    for c in covar.keys():
        model += FixedEffect(covar[c])
    contrast = (group==patient).astype(int) - (group==control).astype(int)
    slm = SLM(model, contrast, correction='fdr', two_tailed=True)
    slm.fit(data)
    return slm

def spatial_correlation(map1, map2, n_rot):
    r = np.corrcoef(map1,map2)[0,1]
    p, null = spin_test(map1, map2, n_rot=n_rot, null_distr=True)

    return r, p, null


def main():

    x = zscore_flip(ct,covar, focus, 'C', 'R')
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

        tle_atrophy[f'{hemi.lower()}tle'] = {'slm':slm, 'r':r, 'p':p, 'null'null}


    print()
    print('Save results')
    print('------------------------')



if __name__ == "__main__":
    main()
