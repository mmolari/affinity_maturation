import numpy as np
import pandas as pd


class dataset:

    def __init__(self, mice_measurements, D_inj, T_delay, meas_protocol):
        self.header_dic = {}
        self.meas_df = pd.DataFrame()
        pass

    def save(self, filename):
        # save a pandas dataframe as csv plus header
        pass

    @staticmethod
    def load(filename):
        # load a pandas dataframe without header
        pass

    def T_delay(self):
        pass

    def D_inj(self):
        pass

    def meas(self):
        pass


def load_all_csv_datasets():
    # with try
    # this should go in the utility folder
    # initialize the utility folder
    pass
