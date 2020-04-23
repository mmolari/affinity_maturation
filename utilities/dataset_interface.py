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


# TODO: load all datasets of scheme 1, scheme 2, scheme 3?
# Reuse parts of the previous function in this case
