Setting model parameters
========================

The model presented in [Molari2020]_ was based on experimental data collected on mice immunization against Tetanus Toxoid. However, depending on the organism / setup one considers one might want to adapt the value of many of the model parameters. To do so, here is a description of the parameters, their standard value in our model and instructions on how to change them.

.. _par_dict:

The parameters dictionary
_________________________

The value of all the model parameters is included in a single python dictionary. In this way values of single parameters are easy to modify and pass to the library functions. Whenever a function accepts a :obj:`par` argument, this consist of the model parameters dictionary.

Here is a table with a non-exhaustive list of the most important model parameters. For a more complete list one can refer to the :file:`model_parameters.py` `file <param_file_>`_. In the table each parameter is accompanied by the corresponding label in the parameters dictionary, its standard [1]_ value and a small description.

.. list-table:: List of main model parameters, their value and description.
   :widths: auto
   :header-rows: 1

   * - Parameter
     - Label
     - St. Value
     - Description
   * - :math:`a`
     - ``a_selection``
     - ``0.12``
     - probability for a B-cell to survive selection for T-cell help irrespective of its affinity.
   * - :math:`b`
     - ``b_selection``
     - ``0.66``
     - probability for a B-cell to fail selection for T-cell help irrespective of its affinity.
   * - :math:`k_B^-`
     - ``k_consumption_per_day``
     - ``2.06e-5``
     - antigen consumption rate by B-cells, in units of 1/days.
   * - :math:`\alpha`
     - ``alpha_C``
     - ``0.023``
     - Antigen concentration to dosage conversion faction, in microgram units.
   * - :math:`\mu_i`
     - ``mu_i``
     - ``-14.59``
     - Mean of the naive population binding energy initial distribution.
   * - :math:`\sigma_i`
     - ``sigma_i``
     - ``1.66``
     - Standard deviation of the naive population binding energy initial. distribution
   * - :math:`N_\text{found}`
     - ``N_founders``
     - ``100``
     - Number of founder clones that seed the Germinal Center.
   * - :math:`N_i`
     - ``N_i``
     - ``2500``
     - Size of the Germinal Center B-cell population at the beginning of the simulation.
   * - :math:`N_{\max}`
     - ``GC_carrying_capacity``
     - ``2500``
     - Germinal Center carrying capacity, i.e. maximum size of the B-cell population.
   * - :math:`g_\text{1d}`
     - ``g_1d``
     - ``0.56``
     - Fraction of Memory Cells in the final responders population, when measurement is done 1 day after boost injection.
   * - :math:`g_\text{4d}`
     - ``g_4d``
     - ``0.``
     - Fraction of Memory Cells in the final responders population, when measurement is done 4 days after the last injection.
   * - :math:`\text{Ag-sel}`
     - ``B_sel``
     - ``False``
     - Whether selection for antigen binding is included in the model.
   * - :math:`\epsilon_\text{Ag}`
     - ``eps_B``
     - ``-13.59``
     - Threshold binding energy for antigen-binding selection.




Create and modify your parameters dictionary
____________________________________________

A dictionary containing the standard value for all parameters can be obtained with the function :func:`st_par`. Values of single parameters can be then modified using their labels::

    # dictionary with the standard value of model parameters
    par = am_sim.st_par()

    # modify one specific parameter
    par['mu_i'] = -13.


.. _exp_range:

Experimental sensitivity range
______________________________

For a good comparison between experimental data and results of the simulations some of the functions of the library take into account the experimental sensitivity range. The limit of this range are set through the parameters :obj:`low_en_exp_cutoff` and :obj:`high_en_exp_cutoff` to be respectively 0.1 and 500 in nano-Molar units, or equivalently -23.03 and -14.51 in :math:`k_B T` energy units. They can be modified by changing their value in :file:`model_parameters.py` `file <param_file_>`_. Cells with lower affinity than 500 nM are considered as non-detectable, while cells with higher affinity than 0.1 nM are detected, but their affinity is by default set to 0.1 nM.



.. rubric:: Footnotes:

.. [1] As inferred in [Molari2020]_ scenario C.

.. _param_file: https://github.com/mmolari/affinity_maturation/blob/master/am_sim/model_parameters.py
