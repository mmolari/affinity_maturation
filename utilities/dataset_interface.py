import numpy as np
import os

import am_sim as ams


data_dir = 'data'


def load_all_datasets():
    '''
    This function loads all the datasets contained in the data_dir
    '''
    # list all the .txt files in the data directory
    files = os.listdir(data_dir)
    files = [x for x in files if x.endswith('.txt')]
    # try to load them as datasets and append them to the dataset list
    ds_list = []
    failed_loading = []
    raise_exception = False
    for fil in files:
        try:
            ds = ams.dataset.load(os.path.join(data_dir, fil))
            ds_list.append(ds)
        except:
            # if loading fails raise an exception for the file
            print(f'WARNING: could not load file {fil}')
            raise_exception = True  # raise an exception later
            failed_loading.append(fil)

    # if at least one file could not be loaded raise an exception
    if raise_exception:
        raise Exception(f'WARNING: could not load file(s): {failed_loading}\n'
                        + 'these files should not be in the data directory.')

    return np.array(ds_list)


def extract_attribute(obj_list, attr, sub_idx=None):
    '''
    Utility function to extract an attribute from a list of objects, and return
    a numpy-array of attributes. The argument 'sub_idx', if specified, selects
    a sub-element with the given index of the returned attribute.
    '''
    attr_list = []  # list containing the specified attribute, to be filled
    for ob in obj_list:
        # get the attribute from the object
        ob_attr = getattr(ob, attr)
        # if specified then extract the sub-element from the attribute
        if sub_idx is not None:
            ob_attr = ob_attr[sub_idx]
        # add it to the list
        attr_list.append(ob_attr)
    # return the list as a numpy array
    return np.array(attr_list)


def extract_scheme_info(dsets):
    '''
    Given a list of datasets extract and returns as numpy arrays the
    measurement protocol, the injected Ag dosage and the injection time delay
    for each dataset.

    Args:
    - dsets (list of dataset objects): list of datasets whose arguments should
        be extracted

    Returns:
    - mp (array of string): list of measurement protocols corresponding to the
        datasets
    - D (array of float): list of Ag dosages of the first injection
    - T (array of int): list of time delays between the firs and second
        injecitons
    '''
    # extract the list of measurement protocols
    mp = extract_attribute(dsets, 'meas_prot')
    # of Ag dosage for the first injection
    D = extract_attribute(dsets, 'D_inj', sub_idx=0)
    # of time delay between the first two injections
    T = extract_attribute(dsets, 'T_delay', sub_idx=0)
    # return the three lists
    return mp, D, T


def dataset_idxs_by_scheme(dsets, scheme_n):
    '''
    Given a list of datasets it returns a list of indices. These indices
    correspond to the dataset objects that belong to the specified scheme.
    They are ordered either based on Ag dosage or injection delay, depending on
    the scheme specified.
    '''
    # extract measurement protocol, Ag dosage and delay between injections for
    # all the datasets in the list. These are returned as np.arrays
    mp, D, T = extract_scheme_info(dsets)

    # create a mask and an order for the datasets, depending on the
    # immunization scheme selected
    if scheme_n == 1:
        mask = mp == '4d'
        order = np.argsort(D)
    elif scheme_n == 2:
        mask = (mp == '1d') & (T == 28)
        order = np.argsort(D)
    elif scheme_n == 3:
        mask = (mp == '1d') & (D == 10)
        order = np.argsort(T)
    else:
        raise Exception('the scheme number must be either 1, 2 or 3.')

    # return a list of selected indices, ordered as specified
    ordered_mask = mask[order]
    selected_idx = order[ordered_mask]
    return selected_idx
