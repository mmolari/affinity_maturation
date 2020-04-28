import numpy as np

data_color = 'C1'


def plot_dset_mean_quantities(ax, dsets, x_variable, qty):
    '''
    This function plots on a given axis object the value of a mean quantity for
    the given dataset list provided. The quantity can either be the mean
    binding energy (qty = 'en') or the high-affinity fraction (qty = 'r_haff').
    This mean quantity is plotted either agains the injected Ag dosage for the
    scheme (x_variable = 'D') or against the injection time delay (x_variable =
    'T'). The average quantity for each single mice is also reported.

    Args:
    - ax (axis object): axis on which the quatities must be plotted
    - dsets (list of datasets): list of datasets for which the average quantity
        should be plotted
    - x_variable (str): either 'D' for dosage or 'T' for injection delay. The
        corresponding quantity will be on the x-axis.
    - qty (str): either 'en' or 'r_haff'. This variable controls which quantity
        is plotted on the y-axis.
    '''
    # x-axis and quantity containers
    x = []
    mean_qty = []

    for ds in dsets:
        # the x-variable corresponding to the dataset considered
        if x_variable == 'D':
            # either dosage
            x_ds = ds.D_inj[0]
        elif x_variable == 'T':
            # or time delay
            x_ds = ds.T_delay[0]
        else:
            raise Exception(
                f' x_variable must be either \'D\' or \'T\', not {x_variable}')

        # plot the average per single mouse
        plot_dset_mice_mean_quantities(ax, ds, x_ds, qty)

        # append the x-variable in the container
        x.append(x_ds)

        # append the mean quantity for the full dataset to the list
        if qty == 'en':
            # either average binding energy
            mean_qty.append(ds.mean_en())
        elif qty == 'r_haff':
            # or high-affinity fraction
            mean_qty.append(ds.r_haff())

    # plot the mean quantity for all the datasets
    ax.plot(x, mean_qty, c=data_color, ls='--',
            marker='o', fillstyle='none', label='data')


def plot_dset_mice_mean_quantities(ax, dset, x_coord, qty):
    '''
    For a single dataset, it plots the value of an average quantity for every
    single mice on the dataset. The quantity can either be the mean binding
    energy (qty = 'en') or the high-affinity fraction (qty = 'r_haff'). The
    result is plotted as a vertical dashed line with crosses as markers. The
    value of the x-coordinate must be specified as input.

    Args:
    - ax (axis object): axis on which the quatities must be plotted
    - dset (dataset object): dataset for which the average quantities should be
        plotted
    - x_coord (float): x-coordinate of the vertical line to be plotted.
    - qty (str): either 'en' or 'r_haff'. This variable controls which quantity
        is plotted on the y-axis.
    '''
    # quantity to be evaluated for every mouse in the dataset.
    if qty == 'en':
        # either mean binding energy
        mean_qty = np.sort(dset.mean_en_per_mouse())
    elif qty == 'r_haff':
        # or high-affinity fraction
        mean_qty = np.sort(dset.r_haff_per_mouse())

    # plot as a vertical dashed lines with crosses as markers.
    x = np.ones_like(mean_qty) * x_coord
    ax.plot(x, mean_qty, c=data_color, ls='--', marker='x')


def plot_dset_hist(ax, ds, bins):
    '''
    This function plots the dataset histogram on the specified axis and with
    the specified binning.

    Args:
    - ax: axis object on which the histogram is plotted.
    - ds (dataset object): dataset object whose energies are to be displayed.
    - bins (list of floats): desired binning of the histogram.
    '''
    ax.hist(ds.all_en(), density=True, bins=bins, alpha=0.4,
            color=data_color, label='data')
