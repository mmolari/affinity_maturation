Simulating Affinity Maturation: a short tutorial
================================================

Here you can find some short documented examples on how to use the `affinity maturation <https://en.wikipedia.org/wiki/Affinity_maturation>`_ simulation library, available in the following `repository <https://github.com/mmolari/affinity_maturation>`_. The model implemented by the code is described in [Molari2020]_. Please don't forget to cite it if you use the code.

The two main topics of this small tutorial are:

- how to simulate an **immunization scheme**, that is a vaccination strategy potentially consisting of multiple antigen injections separated by some interval of time. In particular in our model we try to capture the effect of varying the injected **antigen dose** and the **time delay** between subsequent injections. The aim of the simulation is to predict the quality of immunization, in terms of **antibody affinity** for the antigen.
- how to **infer model parameters** from experimental data, making use of our implementation of the `parallel tempering <https://en.wikipedia.org/wiki/Parallel_tempering>`_ technique.


.. module:: am_sim

Preliminary setup
_________________

All the code is written in ``python3`` (version ``3.7.4``). Running the code also requires the installation of these additional dependencies: ``numpy``, ``scipy`` and ``pandas``. ``matplotlib`` is also required for the plots in this tutorial.

All of the functions and classes needed to perform the immunization simulations are defined in the module :mod:`am_sim`. In order to run the examples shown here one needs to include the ``am_sim`` folder from the `repository`_ in the execution path, and import the module with::

   import am_sim


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   parameters
   pop_fun
   germ_cent
   immscheme
   dataset
   inference


.. rubric:: Bibliography

.. [Molari2020] *Quantitative modeling of the effect of antigen dosage on the distribution of B-cell affinities in maturating germinal centers*, by Marco Molari, Klaus Eyer, Jean Baudry, Simona Cocco, and Rémi Monasson


.. _repository: https://github.com/mmolari/affinity_maturation
