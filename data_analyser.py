import numpy as np


def find_nu_dnu(m_val, n_val, frequency_df):
    """
    Finds nu and dnu for a given m and n in a given dataset

    Parameters
    ----------
    m_val : int
        The azimuthal degree
    n_val : int
        The radial order
    frequency_df : pandas.DataFrame
        The dataset as a pandas dataframe

    Returns
    -------
    nu, dnu : float, float or (np.nan, np.nan)
        The nu and dnu values for the given m_val and n_val. If no valid rows are found, returns np.nan, np.nan
    """

    index = (m_val, n_val)
    if index in frequency_df.index:
        #print('m, n, l')
        return frequency_df.loc[index, ["nu", "dnu"]]
    else:
        return np.nan, np.nan



def find_all_frequencies(m_val, n_val, all_frequency_datasets):
    """
    Finds nu and dnu for a given m and n in all datasets in a given directory
    for a given m and n

    Parameters
    ----------
    m_val : int
        The azimuthal degree
    n_val : int
        The radial order
    all_frequency_datasets : dict
        A dictionary of all datasets in a given directory

    Returns
    -------
    array_nu, array_dnu : array, array
        The nu and dnu values for the given m_val and n_val in all datasets in a given directory
    """
    
    array_nu = []
    array_dnu = []
    for dataset in all_frequency_datasets.values():
        nu, dnu = find_nu_dnu(m_val, n_val, dataset)
        array_nu.append(nu)
        array_dnu.append(dnu)
    return array_nu, array_dnu
