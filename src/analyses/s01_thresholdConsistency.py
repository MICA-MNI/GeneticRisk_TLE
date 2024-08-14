import numpy as np
from 03_casecontrolAssociation import casecontrol_difference, spatial_correlation
from 01_geneticCorrelation import zscore_flip, imaging_genetic_association

def main():

    thresholds = ['Pt_0.00100005', 'Pt_0.0500001', 'Pt_0.1', 'Pt_0.2', 'Pt_0.3', 'Pt_0.4', 'Pt_0.5']

    print('------------------------------------------')
    print('Effects of PRS across different thresholds')
    print('------------------------------------------')

    print()
    print('Brain-wide associations')
    print('------------------------------------------')
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
    lobe_names = ['whole', 'frontal', 'limbic', 'occipital', 'parietal', 'temporal']
    f = [18, 19, 26, 17, 23, 3, 27, 11, 13, 31, 16, 47, 65, 61, 50, 45, 60, 53, 51, 57, 37, 52]
    l = [9, 22, 2, 25, 34, 59, 36, 56, 43, 68];
    o = [12, 20, 10, 4, 46, 54, 38, 44]
    p = [21, 28, 30, 7, 24, 55, 62, 41, 64, 58]
    t = [5, 6, 8, 14, 15, 29, 32, 33, 1, 42, 48, 63, 39, 66, 49, 40, 67, 35]
    lobes = np.zeros(68)
    lobes[f] = 1; lobes[l] = 2; lobes[o] = 3; lobes[p] = 4; lobes[t] = 5
    lobes = parcel_to_surface(lobes, 'aparc_conte69')

    # Correlate with PRS-TLE
    for i, lobe in enumerate(lobe_names):
        if 'whole' in lobe:
            idx = np.array([True]*int(len(lobes)/2) + [False]*int(len(lobes)/2)) if 'lh' in lobe else np.array([False]*int(len(lobes)/2) + [True]*int(len(lobes)/2))
        else:
            idx = lobes == i if 'lh' in lobe else lobes == i+6

        r = np.zeros(len(thresholds))
        r2 = np.zeros(len(thresholds))
        p = np.zeros(len(thresholds))
        for thr in thresholds:
            prs = prs_all[thr]
            corr = pearsonr(prs, mean_ct)
            r[i], r2[i], p[i] = corr.statistic, corr.statistic**2, corr.pvalue

        global_association[lobe] = pd.DataFrame('thresholds':thresholds, 'r':r, 'r2':r2, 'p':p)

        print(f'{lobe} lobe done')

    print()
    print('Regional effects across d')
    print('------------------------------------------')
    for thr in thresholds:
        prs = prs_all[thr]
        term_prs = FixedEffect(prs)
        model = term_age + term_sex + term_pc10 + term_prs
        contrast = prs
        slm = SLM(model,contrast,correction='fdr',two_tailed=True)
        slm.fit(ct_aparc)
        regional_association[thr] = slm

        print(f'{thr} done')

    print()
    print('Save results')
    print('------------------------------------------')

    print()
    print('----------------------------------')
    print('Imaging-genetic pattern similarity')
    print('----------------------------------')

    print()
    print('Psnp thresholds')
    print('----------------------------------')
    thr_similarity_r = np.zeros((len(thresholds),len(thresholds)))
    thr_similarity_p = np.zeros((len(thresholds),len(thresholds)))
    for i in range(len(thresholds)):
        thr1 = regional_association[thresholds[i]]
        for j in range(i, len(thresholds)):
            thr2 = regional_association[threhsolds[j]]
            thr_similarity_r[i,j], thr_similarity_p[i,j], _ = spatial_correlation(thr1.t, thr2.t, n_rot=5000)
            print(f'{thresholds[i]} x {thresholds[j]} done')

    print()
    print('Comparison to case-control atrophy')
    print('----------------------------------')
    # Generate atrophy for each site:
    local_tle_ct = pd.read_csv('../data/raw/local_tle_ct.csv', index_col='participant').to_numpy()

    # Combined
    tle_atrophy['lh_all'] = local_tle_atrophy['ltle']['slm']
    tle_atrophy['rh_all'] = local_tle_atrophy['rtle']['slm']

    # Across different datasets
    for d in ['EpiC', 'MICs', 'NKG']:
        dataset_idx = dataset==d
        dataset_ct = ct[dataset_idx,:]
        dataset_age = age[dataset_idx]
        dataset_sex = sex[dataset_idx]
        dataset_focus = focus[dataset_idx]

        x = zscore_flip(dataset_ct, dataset_focus, 'C', 'R')
        for hemi in ['L', 'R']:
            covar = pd.concate(dataset_age, dataset_sex)
            tle_atrophy[f'{hemi.lower()}h_{d.lower()}'] = casecontrol_difference(dataset_ct, dataset_focus, 'C', hemi)

        print(f'{d} done')

    site_labels = ['lh_all', 'lh_epic', 'lh_mics', 'lh_nkg',
                   'rh_all', 'rh_epic', 'rh_mics', 'rh_nkg']
    tle_similarity_r = np.zeros((len(thresholds),len(site_labels)))
    tle_similarity_p = np.zeros((len(thresholds),len(site_labels)))
    for i in range(len(thresholds)):
        map1 = regional_association[thresholds[i]]
        for site in site_labels:
            map2 = tle_atrophy[site]
            tle_similarity_r[i,j], thr_similarity_p[i,j], _ = spatial_correlation(map1.t, map2.t, n_rot=5000)

    print()
    print('Save results')
    print('----------------------------------')
    np.savez('../../data/results')

    
if __name__ == "__main__":
    main()
