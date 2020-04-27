# Affinity maturation simulation

This repository contains the code for the simulations and figures featured in the paper *Quantitative modeling of the effect of antigen dosage on the distribution of B-cell affinities in maturating germinal centers Supplementary Information*, by Marco Molari, Klaus Eyer, Jean Baudry, Simona Cocco, and Rémi Monasson.
All of the code contained has been written in Python3 (version 3.7.4).
Dependencies required to run the code include `numpy`, `scipy`, `pandas`, `matplotlib`. Moreover an installation of [jupyter](https://jupyter.org) is needed to run the interactive notebooks.

## Repository structure:

- The folder `am_sim` contains the library with most of the code core functions
- The folder `data` contains experimental affinity measurements for each immunization scheme. The measurements are saved in a `.csv` file, one per immunization scheme considered. The `.csv` file has an header containing details of the immunization scheme.
- The folder `utilities` contains the definition of other functions used to create the figures of the paper.
- The `parallel_tempering.py` file contains the script used to perform the maximum-likelihood inference.
- The `inference_results` folder contains the results of the inference procedure. These results can be reproduced by running the `parallel_tempering.py` file.
- The `figure_*.ipynb` are interactive [jupyter](https://jupyter.org) notebooks that contain code to create the figures of the main paper.
- The `figures` folder is used to save the produced figures and some of the simulation results needed to generate them.

## `am_sim` library structure

Here I give a short overview of the main classes and functions contained in the `am_sim` library. For a more detailed description one should refer directly to the code, which is well-commented.


- `st_par` :
- `antigen` :
- `stoch_pop` :
- `det_pop` :
- `germinal_center` :
- `simulate_immscheme` :
- `dataset` :
- `responders_from_dset` :
- `parallel_tempering` :

## Running the inference procedure

## Producing the figures of the paper

## TODO:

- [ ] Clean up README
- [ ] Check code documentation
- [ ] Describe how to run the inference and what are the results
- [ ] Describe
