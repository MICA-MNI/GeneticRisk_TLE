import numpy as np
from enigmatoolbox.permutation_testing import spin_test


def load_imaging_genetic():
    """
    Load imaging genetic results from 01_geneticCorrelation.

    Returns:
    slm: Statistical linear model object
    """
    file_path = f"{os.path.dirname(os.path.abspath(__file__))}/../../data/results/01_geneticCorrelation/regional_association.pkl"
    with open(file_path, "rb") as f:
        data = pickle.load(f)
    return data["slm"]


def load_data(file_path, variables):
    """
    Load specified variables from a .npz file.

    Parameters:
    file_path (str): The path to the .npz file.
    variables (list of str): A list of variable names to extract from the .npz file.

    Returns:
    tuple: A tuple containing the data for each specified variable in the order they are listed in the variables parameter.
    """
    data = np.load(file_path, allow_pickle=True)
    return tuple(data[var] for var in variables)


def spatial_correlation(
    map1, map2, n_rot=5000, surface_name="fsa5", parcellation_name="aparc"
):
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
    r = np.corrcoef(map1, map2)[0, 1]
    p, null = spin_test(
        map1,
        map2,
        surface_name=surface_name,
        parcellation_name=parcellation_name,
        n_rot=n_rot,
        null_distr=True,
    )

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
    hc = group == control
    hc_mean = np.mean(data[hc, :], axis=0)
    hc_std = np.mean(data[hc, :], axis=0)
    z = (data - hc_mean) / hc_std

    # Sort based on ipsilateral/contralateral
    rtle = group == flip
    n_region = data.shape[1] / 2
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
    contrast = (group == patient).astype(int) - (group == control).astype(int)
    slm = SLM(model, contrast, correction="fdr", two_tailed=True)
    slm.fit(data)
    return slm
