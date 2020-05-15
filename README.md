# Affinity maturation simulation

This repository contains the code for the simulations and figures featured in the paper *Quantitative modeling of the effect of antigen dosage on the distribution of B-cell affinities in maturating germinal centers Supplementary Information*, by Marco Molari, Klaus Eyer, Jean Baudry, Simona Cocco, and Rémi Monasson.
For any question please write to <marco.molari@phys.ens.fr>.

All of the code contained has been written in Python3 (version 3.7.4).

Dependencies required to run the code include `numpy`, `scipy`, `pandas`, `matplotlib`. Moreover an installation of [jupyter](https://jupyter.org) is needed to run the `*.ipynb` notebooks.

## Repository structure:

- The folder `am_sim` contains the library with the definition of most of the core classes and functions.
- The folder `data` contains experimental affinity measurements for each immunization scheme. See "format of the data" below for more details.
- The folder `utilities` contains the definition of other functions used to create the figures of the paper.
- The `parallel_tempering.py` file contains the script used to perform the maximum-likelihood inference.
- The `inference_results` folder contains the results of the inference procedure. These results can be reproduced by running the `parallel_tempering.py` file.
- The `figure_*.ipynb` are interactive [jupyter](https://jupyter.org) notebooks that contain code to create the figures of the main paper.
- The `figures` folder is used to save the produced figures the simulation results.

## Running the inference procedure

The results of our inference procedure, needed for the figures of the paper, is already present in the folder `inferece_results`. However one can reproduce these results by launching the python script `parallel_tempering.py`. This script will automatically create a new folder, named `reproduce_inference_results` in which the results will be saved. The new folder will contain the following file:

- `search_setup.txt` : a text file containing details on the initial setup of the search procedure.
- `t_*_search_history.csv` : a pandas dictionary containing the results of the search up to time t. It is updated every 100 parallel tempering steps.
- `par_i.pkl` : a pickle file containing a dictionary with the initial values of all model parameters in the search.
- `dataset_list.pkl` : a pickle file containing the list of the datasets with which the likelihood is evaluated. These are all the datasets contained in `data`, and this file is produced only as a safety check.

## Producing the figures of the paper

To produce the figures of the paper

## Format of the data

The measurements are saved in a `.csv` file, one per immunization scheme and condition considered. The `.csv` file has an header containing details of the immunization scheme (injected antigen dosage in micrograms per injection, delays between injections and between injection and measurements, measurement protocol used) and each column contains data for a different mouse.

## code documentation

If you are interested in using the code for germinal center simulations or applying the inference procedure to your data, here is a documentation page with more details on how to use the code.

Please if you do so remember to cite our paper.
