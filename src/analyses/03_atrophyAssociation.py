import numpy as np
import pandas as pd
import pickle
from 02_epicentreMapping import spatial_correlation
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM


# Helper functions
def zscore_flip(data, group, control, flip):
    """
    Perform z-score normalization (relative to healthy controls) and sort patients based on ipsilateral/contralateral

    Parameters:
    data (array-like): Case-control data
    group (array-like): Group labels
    control (str): Control-identifying label
    flip (str): Right-identifying label

    Returns:
    f (array-like): z-score normalized and sorted (ipsi/contra) data
    """
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
    """
    Compute between-group differents using a fixed effect model

    Parameters:
    data (array-like): DataFram
    covar (dict): Covariates in the model
    group: Group labels
    control: Control-identifying label
    patient: Patient-identifying labels

    Returns:
    slm (SLM): Statistical linear model
    """
    model = FixedEffect(group)
    for c in covar.keys():
        model += FixedEffect(covar[c])
    contrast = (group==patient).astype(int) - (group==control).astype(int)
    slm = SLM(model, contrast, correction='fdr', two_tailed=True)
    slm.fit(data)
    return slm


# Main analysis
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
        atrophy[f'{hemi.lower()}tle'] = slm

        print()
        print('Correlate with PRS map')
        print('------------------------')
        r, p, null = spatial_correlation(slm.t, imaging_genetics.t)

        atrophy_association[f'{hemi.lower()}tle'] = {'r':r, 'p':p, 'null':null}


    print()
    print('------------')
    print('Save results')
    print('------------')
    np.savez('../../data/results/03_atrophyAssociation/atrophy_association.npz',
            atrophy=atrophy, atrophy_association=atrophy_association)


if __name__ == "__main__":
    main()
