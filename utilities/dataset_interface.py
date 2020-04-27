import am_sim as ams

import os


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

    return ds_list


def select_scheme_1_ds(dsets):
    '''Given the list of datasets it selects only the ones belonging to
    immunization scheme 1.'''

    # select all datasets with measurement protocol '4d'
    s1ds = [ds for ds in dsets if ds.meas_prot == '4d']
    # order by dosage
    order = np.argsort([ds.D_inj[0] for ds in s1ds])
    s1ds = np.array(s1ds)[order]
    return s1ds


def select_scheme_2_ds(dsets):
    '''Given the list of datasets it selects only the ones belonging to
    immunization scheme 2.'''
    # select all datasets with measurement protocol '1d'
    s2ds = [ds for ds in dsets if ds.meas_prot == '1d']
    # select all datasets with 28 days delay between two injections
    s2ds = [ds for ds in s2ds if ds.T_delay[0] == 28]
    # order by dosage
    order = np.argsort([ds.D_inj[0] for ds in s2ds])
    s2ds = np.array(s2ds)[order]
    return s2ds


def select_scheme_3_ds(dsets):
    '''Given the list of datasets it selects only the ones belonging to
    immunization scheme 3.'''
    # select all datasets with measurement protocol '1d'
    s3ds = [ds for ds in dsets if ds.meas_prot == '1d']
    # select all datasets with 10 microgram Ag injection in the first injection
    s3ds = [ds for ds in s3ds if ds.D_inj[0] == 10]
    # order by injection delay
    order = np.argsort([ds.T_delay[0] for ds in s3ds])
    s3ds = np.array(s3ds)[order]
    return s3ds


def select_scheme_ds(dsets, scheme_n):
    '''Given the list of datasets it selects only the ones belonging to
    immunization scheme 'scheme_n', where 'scheme_n' is either 1,2 or 3.'''
    if scheme_n == 1:
        return select_scheme_1_ds(dsets)
    elif scheme_n == 2:
        return select_scheme_2_ds(dsets)
    elif scheme_n == 3:
        return select_scheme_3_ds(dsets)
    else:
        raise Exception('scheme number must be 1,2 or 3')
