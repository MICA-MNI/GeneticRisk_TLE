import numpy as np
import pandas as pd
import pickle
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM
from scipy.stats.stats import pearsonr

def correlate_lobes(lobes, hemisphere_lobes, ct_vertex, prs, hemisphere):
    mean_ct = np.zeros(len(hemisphere_lobes))
    r = np.zeros(len(hemisphere_lobes))
    r2 = np.zeros(len(hemisphere_lobes))
    p = np.zeros(len(hemisphere_lobes))


def main():

    # Load data
    data = np.load('../data/processed/abcd_data.npz')
    thr = 'Pt_0.1'; prs = prs_all[thr]

    print('-----------------------------')
    print('Brain-wide effects of PRS-TLE')
    print('-----------------------------')

    print()
    print('Correlate PRS with each lobe')
    print('-----------------------------')
    # Regress out effects of age, sex, and term_pc10
    term_age = FixedEffect(age)
    term_sex = FixedEffect(sex)
    term_pc10 = FixedEffect(pc10)
    model = term_age + term_sex + term_pc10
    contrast = np.ones(len(prs))
    slm = SLM(model, contrast, surf='fsLR32k')
    slm.fit(ct_vertex)
    residual = ct_vertex - np.dot(slm.X, slm.coef)

    # Lobe indices
    lobe_names = ['lh_whole', 'lh_frontal', 'lh_limbic', 'lh_occipital', 'lh_parietal', 'lh_temporal',
                  'rh_whole', 'rh_frontal', 'rh_limbic', 'rh_occipital', 'rh_parietal', 'rh_temporal']
    fl = [18, 19, 26, 17, 23, 3, 27, 11, 13, 31, 16]; fr = [47, 65, 61, 50, 45, 60, 53, 51, 57, 37, 52]
    ll = [9, 22, 2, 25, 34]; lr = [59, 36, 56, 43, 68];
    ol = [12, 20, 10, 4]; or_ = [46, 54, 38, 44]
    pl = [21, 28, 30, 7, 24]; pr = [55, 62, 41, 64, 58]
    tl = [5, 6, 8, 14, 15, 29, 32, 33, 1]; tr = [42, 48, 63, 39, 66, 49, 40, 67, 35]
    lobes = np.zeros(68)
    lobes[fl] = 1; lobes[ll] = 2; lobes[ol] = 3; lobes[pl] = 4; lobes[tl] = 5
    lobes[fr] = 7; lobes[lr] = 8; lobes[or] = 9; lobes[pr] = 10; lobes[tr] = 11
    lobes = parcel_to_surface(lobes, 'aparc_conte69')

    lobe_names = ['lh_whole', 'lh_frontal', 'lh_limbic', 'lh_occipital', 'lh_parietal', 'lh_temporal',
                  'rh_whole', 'rh_frontal', 'rh_limbic', 'rh_occipital', 'rh_parietal', 'rh_temporal']

    # Correlate with PRS-TLE
    mean_ct = np.zeros(len(prs),len(lobe_names))
    r = np.zeros(len(lobe_names))
    r2 = np.zeros(len(lobe_names))
    p = np.zeros(len(lobe_names))

    for i, lobe in enumerate(lobe_names):
        if 'whole' in lobe:
            idx = np.array([True]*(len(lobes)/2) + [False]*(len(lobes)/2)) if 'lh' in lobe else np.array([False]*(len(lobes)/2) + [True]*(len(lobes)/2))
        else:
            idx = lobes == i if 'lh' in lobe else lobes == i+6

        mean_ct[i] = np.mean(ct_vertex[:, idx], axis=1)
        corr = pearsonr(prs, mean_ct)
        r[i], r2[i], p[i] = corr.statistic, corr.statistic**2, corr.pvalue
        print(f'{lobe}: r = {r[i]}, r-squared = {r2[i]}, p = {p[i]}')

    global_association = pd.DataFrame('lobe':lobe_names, 'r':r, 'r2':r2, 'p':p)

    print()
    print('Save results')
    print('-----------------------------')
    global_association.to_pickle('../../data/results/01_geneticCorrelation/global_association.pkl')

    print()
    print('-------------------------')
    print('Local effects of PRS-TLE')
    print('------------------------')

    print()
    print('Correlate across regions')
    print('------------------------')
    term_prs = FixedEffect(prs)
    model = term_age + term_sex + term_pc10 + term_prs
    contrast = prs
    slm = SLM(model,contrast,correction='fdr',two_tailed=True)
    slm.fit(ct_aparc)

    print()
    print('Save results')
    print('-----------------------------')
    pickle.dump(slm, '../../data/resuts/01_geneticCorrelation/regional_association.pkl')


if __name__ == "__main__":
    main()
