import numpy as np
import pandas as pd
import ast

from .model_parameters import high_en_exp_cutoff, low_en_exp_cutoff


class dataset:
    '''
    This class implements a container that holds simultaneously:
    - details of the immunization protocol used
    - experimental affinity measurements

    Attributes:
    - D_inj (list of floats): list of dosages for every injection performed
    - T_delay (list of int): list of delays between injections
    - meas_prot (string): wether the measurement is performed 1 day after boost
        or 4 days after recall injection
    - aff_nM (list of list): each sub-list contains affinity measurements
        related to a single mouse. Affinity measurements are in nanoMolar units.
    - en (list of list): same as aff_nM, but affinities are reported as binding
        energies.

    Member functions:
    - save: saves the dataset into a csv-like file
    - load (static): loads the dataset from a save file
    - all_en: returns the list of binding-energy measurements for all mice
    '''

    def __init__(self, nM_affinity_table, D_inj, T_delay, meas_prot):
        '''
        Dataset class initializer.

        Args:
        - nM_affinity_table (list of list): it contains the affinity
            measurements in nanoMolar units. This argument must be a list of
            lists, each sub-list containing measurements for a single mouse.
        - D_inj (list of float): list of injected Ag dosages in micrograms, one
            per injection.
        - T_delay (list of int): list of delays between injections in days. The
            last element must be the delay between the last injection and the
            experimental measurement.
        - meas_prot (str): either '1d' if measurement is performed 1 day after
            boost or '4d' for 4 days after recall injection.
        '''
        # immunization scheme information
        self.D_inj = D_inj
        self.T_delay = T_delay
        self.meas_prot = meas_prot
        # affinity measurements, in nanoMolar and binding energy units
        self.aff_nM = [np.array(mouse_aff) for mouse_aff in nM_affinity_table]
        # order mice by number of measurement per mouse
        N_meas_mice = [mouse_aff.size for mouse_aff in self.aff_nM]
        order = np.argsort(N_meas_mice)[::-1]
        self.aff_nM = list(np.array(self.aff_nM)[order])
        # save affinities also as binding energies
        self.en = [np.log(mouse_aff * 1e-9) for mouse_aff in self.aff_nM]

        # check that the experimental sensitivity limits are respected:
        for mouse_en in self.en:
            above = np.any(mouse_en < low_en_exp_cutoff)
            below = np.any(mouse_en > high_en_exp_cutoff)
            if above or below:
                raise Exception('warning: experimental measurements exceed\
                specified instrumental sensitivity.')

    def __str__(self):
        '''
        Creates a string representation of the dataset, reporting details of
        the immunization scheme and number of measurements.
        '''
        str_repr = f'inj_dosage    : {self.D_inj}\n'
        str_repr += f'T_delay       : {self.T_delay}\n'
        str_repr += f'meas_protocol : {self.meas_prot}\n'
        str_repr += f'N_tot_meas    : {self.all_en().size}\n'
        for n_mouse, mouse_en in enumerate(self.en):
            str_repr += f'N_meas_mouse_{n_mouse}: {mouse_en.size}\n'
        return str_repr

    def __repr__(self):
        '''
        Creates a short string representation for the dataset object, reporting
        details of the immunization scheme and number of measurements.
        '''
        str_repr = f'D : {str(self.D_inj):<15} '
        str_repr += f'T : {str(self.T_delay):<15} '
        str_repr += f'mp : {self.meas_prot:<4} '
        str_repr += f'N : {self.all_en().size}\n'
        return str_repr

    def save(self, filename):
        '''
        Saves the dataset into a csv-like file.

        Args:
        - filename (str): full file path of the file into which the dataset
            must be saved
        '''
        # write header lines to the file specifying the immunization scheme
        # parameters
        with open(filename, 'w') as f:
            f.write(f'inj_dosage    : {self.D_inj}\n')
            f.write(f'T_delay       : {self.T_delay}\n')
            f.write(f'meas_protocol : {self.meas_prot}\n')
            f.close()

        # save the affinity measurement as a pandas dataframe in csv format
        df = pd.DataFrame(self.aff_nM).T
        df.columns = [f'mouse_{n+1}' for n in range(df.columns.size)]
        df.to_csv(filename, mode='a', index=False)

    @staticmethod
    def load(filename):
        '''
        Loads the dataset from a file.

        Args:
        - filename (str): the full file name and path from which the dataset
            should be loaded

        Returns:
        - ds (dataset object): the loaded dataset
        '''
        # load immunization protocol details from first three lines
        with open(filename, 'r') as f:
            D_inj = ast.literal_eval(f.readline()[15:].strip())
            T_del = ast.literal_eval(f.readline()[15:].strip())
            me_pr = f.readline()[15:].strip()
            f.close()

        # load affinity measurments from the rest of the file
        df = pd.read_csv(filename, index_col=False, skiprows=3)
        nM_aff = [x[~np.isnan(x)] for x in df.T.to_numpy()]

        # create and return dataset
        ds = dataset(nM_affinity_table=nM_aff, D_inj=D_inj,
                     T_delay=T_del, meas_prot=me_pr)
        return ds

    def all_en(self):
        '''
        Function that returns the cumulative list of all binding energy
        measurements for all the mice in the dataset.

        Returns:
        - all_en (list of floats): cumulative list of all the experimental
            affinity measurements in the protocol, in binding energy units.
        '''
        # return cumulative measurements
        return np.concatenate(self.en)
