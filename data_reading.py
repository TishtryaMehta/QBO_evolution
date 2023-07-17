import os
import glob
import pandas as pd
import numpy as np


###########################################################
#  OBTAINS A LIST OF ALL THE FILES IN THE GIVEN DIRECTORY
#  AND SORTS THEM BY CHRONOLOGICAL ORDER
#  RETUNRS: ARRAY OF FILENAMES IN STR FORMAT
###########################################################
def filename_finder(data_path):
    filenames = []
    for num in [9, 0, 1, 2]:
        filenames.extend(sorted(glob.iglob(os.path.join(data_path, f"mrv1f{num}"+"*.txt"),recursive=True)))
    #print(len(filenames), "length of filenames")
    return filenames

def filename_finder_ell(data_path, ell):
    filenames = []
    for num in [9, 0, 1, 2]:
        filenames.extend(sorted(glob.iglob(os.path.join(data_path, f"mrv1f{num}"+"*d"+f"{ell:03d}"+".txt"),recursive=True)))
    #print(len(filenames), "length of filenames")
    return filenames

def read_frequency_dataset(directory_path, folder_path, dataset_name, dataset_header, skip_rows):
    """
    Reads a frequency dataset at the given path and returns a pandas dataframe

    Parameters
    ----------
    directory_path : str
        The path to the directory containing the dataset
    dataset_name : str
        The path to the dataset relative to the directory
    dataset_header : list of str
        The header of the dataset
    Returns
    -------
    pandas.DataFrame
        The dataset as a pandas dataframe
    """
    datapath = os.path.join(directory_path, folder_path, dataset_name)
    frequency_df = pd.read_csv(datapath, delim_whitespace=True, skiprows=skip_rows, names=dataset_header)
    return frequency_df.loc[(frequency_df["bad"] == 0) & (frequency_df["ierr"] == 0), :].groupby(["m", "n"]).first()


###########################################################
#  SAVES ARRAY OF NU AND DNU FOUND BY find_all_frequencies
#  FOR A GIVEN M AND N IN A GIVEN DIRECTORY. SAVES TO EITHER
#  GOOD OR BAD DATA DEPENDING ON WHETHER NAN IS FOUND
#  RETURNS: NONE
###########################################################
def savedata(array_nu, array_dnu, m_val, n_val, datapath_save_good, datapath_save_bad, ell):
    df_array = pd.DataFrame({"array_nu": array_nu, "array_dnu": array_dnu})
    filename = "v1f_m_" + str(m_val) + "_n_" + str(n_val)
    if np.count_nonzero(np.isnan(array_nu)) != 0:
        # print(np.count_nonzero(np.isnan(array_nu)), ' NANs were found')
        df_array.to_csv(datapath_save_bad + filename + "_"+f"{ell:03d}"+".csv", index=False)
    else:
        # print('No nans were found')
        df_array.to_csv(datapath_save_good + filename + "_"+f"{ell:03d}"+".csv", index=False)


###########################################################
#  SAVES ARRAY OF NU AND DNU FOUND BY find_all_frequencies
#  FOR A GIVEN M AND N IN A GIVEN DIRECTORY. SAVES TO EITHER
#  GOOD OR BAD DATA DEPENDING ON WHETHER NAN IS FOUND
#  RETURNS: NONE
###########################################################
def savedata_as_np(array_nu, array_dnu, m_val, n_val, ell, datapath_save_np):
    df_array = pd.DataFrame({"array_nu": array_nu, "array_dnu": array_dnu})
    filename = "v1f_m_" + str(m_val) + "_n_" + str(n_val)
    if np.count_nonzero(np.isnan(array_nu)) != 0:
        # print(np.count_nonzero(np.isnan(array_nu)), ' NANs were found')
        df_array.to_csv(datapath_save_bad + filename + "_"+f"{ell:03d}"+".csv", index=False)
    else:
        # print('No nans were found')
        df_array.to_csv(datapath_save_good + filename + "_"+f"{ell:03d}"+".csv", index=False)
