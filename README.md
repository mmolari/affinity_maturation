[![Documentation Status](https://readthedocs.org/projects/affinity-maturation/badge/?version=latest)](https://affinity-maturation.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/github/license/mmolari/affinity_maturation)](https://opensource.org/licenses/MIT)
[![Python version](https://img.shields.io/badge/python-3.7-blue)](https://www.python.org/downloads/)

# Affinity maturation simulation

This repository contains the code for the simulations and figures featured in the paper [*Quantitative modeling of the effect of antigen dosage on the distribution of B-cell affinities in maturating germinal centers*](https://elifesciences.org/articles/55678), by Marco Molari, Klaus Eyer, Jean Baudry, Simona Cocco, and Rémi Monasson, currently available on eLife.
Please don't forget to cite this paper if you use the code for your project.
For questions concerning the code please write to <marco.molari@phys.ens.fr>.

All of the code is written in Python3 (version 3.7.4).

Dependencies required to run the code include `numpy`, `scipy`, `pandas`, `matplotlib`. Moreover an installation of [jupyter](https://jupyter.org) is needed to run the `*.ipynb` notebooks.

## Repository structure:

- The folder `am_sim` contains the library with the definition of most of the core classes and functions.
- The folder `data` contains experimental affinity measurements for each immunization scheme. See [format of the data](#format-of-the-data) for more details.
- The folder `utilities` contains the definition of other functions used to create the figures of the paper.
- The `parallel_tempering.py` file contains the script used to perform the maximum-likelihood inference.
- The `inference_results` folder contains the results of the inference procedure. These results can be reproduced by running the `parallel_tempering.py` file, see [running the inference procedure](#running-the-inference-procedure).
- The `figure_*.ipynb` are interactive [jupyter](https://jupyter.org) notebooks that contain code to create the figures of the main paper.
- The `figures` folder contains sub-folders, one for each figure, with the `.pdf` version of the figures. For simplicity we also include some pre-generated results of the stochastic simulations, but these results can also be recreated with the code included in the notebooks.
- The `docs` folder contains a short tutorial on how to use the code, see [a short tutorial](#how-to-implement-your-own-simulations-a-short-tutorial)

## Running the inference procedure

The results of our inference procedure, needed for the figures of the paper, is already present in the folder `inferece_results`. However one can reproduce these results by launching the python script `parallel_tempering.py`. This script will automatically create a new folder, named `reproduce_inference_results` in which the results will be saved. The new folder will contain the following files:

- `search_setup.txt` : a text file containing details on the initial setup of the search procedure, including a list of the initial model parameters value.
- `t_*_search_history.csv` : a pandas dictionary containing the results of the search up to time t. It is updated every 100 parallel tempering steps.
- `par_i.pkl` : a pickle file containing a dictionary with the initial values of all model parameters in the search, produced as a safety check.
- `dataset_list.pkl` : a pickle file containing the list of the datasets of which the likelihood is evaluated. These are all the datasets contained in the `data` folder, and this file is produced only as a safety check.

## Reproducing the figures of the paper

To produce the figures of the paper please run the corresponding `figure_*.ipynb` jupyter notebook. Some of the figures require running time-intensive stochastic simulations. In this case the user has the choice of skipping the notebook cell containing the simulations and load pre-generated results from the next cell.

For a quick preview of the figures and notebooks:

- [figure 2](https://nbviewer.jupyter.org/github/mmolari/affinity_maturation/blob/master/figure_2.ipynb)
- [figure 3](https://nbviewer.jupyter.org/github/mmolari/affinity_maturation/blob/master/figure_3.ipynb)
- [figure 4](https://nbviewer.jupyter.org/github/mmolari/affinity_maturation/blob/master/figure_4.ipynb)
- [figure 5](https://nbviewer.jupyter.org/github/mmolari/affinity_maturation/blob/master/figure_5.ipynb)
- [figure 6](https://nbviewer.jupyter.org/github/mmolari/affinity_maturation/blob/master/figure_6.ipynb)

## Format of the data

The measurements are saved in `.csv` files, one per immunization scheme and condition considered. Each `.csv` file has a header containing details of the immunization scheme. This includes:
- `inj_dosage` : a list containing the amount of injected antigen in micrograms for each injection.
- `T_delay` : a list of time delays in days. These represent delays between injections, or delay between injection and measurement in the case of the last number.
- `meas_protocol` : wether the measurement is performed 1 day after boost or 4 days after injection.

The rest of the `.csv` file contains a table. Each column corresponds to a different mouse. Entries consist in single-cell affinity measurements. The measured quantity is the antibody dissociation constant Kd, measured in nano-Molar units. For details on the experimental technique used to obtain the measurements please refer to [*Eyer, Klaus, et al. "Single-cell deep phenotyping of IgG-secreting cells for high-resolution immune monitoring." Nature biotechnology 35.10 (2017): 977.*](https://www.nature.com/articles/nbt.3964)

## How to implement your own simulations: a short tutorial

We include in the repository a small tutorial on how to implement custom-made immunization scheme simulations, and how to import new experimental datasets and run on them the inference procedure. This cna be accessed online at [ReadTheDocs](https://affinity-maturation.readthedocs.io/), or alternatively a local version can be generated by running `make html` from the `docs` directory.
