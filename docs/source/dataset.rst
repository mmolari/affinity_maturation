Import your dataset
===================

Experimental measurements collected under a particular condition are collected in a :class:`~am_sim.dataset` object. This object also includes information on the immunization scheme used.

Initialization of a dataset object is done by specifying a list of affinity measurements (in nano-Molar units) collected as a list of lists [#f1]_ ::

   mouse_1 = [21, 26, 12] # experimental measurements for mouse 1
   mouse_2 = [25, 34, 13, 45] # experimental measurements for mouse 2

   aff_meas = [mouse_1, mouse_2] # list of affinity measurements for each mouse

   D_inj = [1.,1.] # protocol consisting of two injections of 1 microgram of Ag each
   T_delay = [28, 28] # 4 weeks of delay between the two injections, cells measured 4 weeks later
   mp = '1d' # string describing the measurement protocol used. Here cells are measured 1 day after boost

   # creating the dataset objects
   ds = am_sim.dataset(aff_meas, D_inj, T_delay, mp)


Once created datasets can be quickly saved and loaded using the appropriate :meth:`~am_sim.dataset.save` and :meth:`~am_sim.dataset.load` class methods. They take as argument the file name, and results are stored in .csv format::

   # saving the dataset on a file
   ds.save('my_dataset.csv')

   # loading the dataset from a file
   ds = am_sim.dataset.load('my_dataset.csv')

The resulting file contains a header with details of the immunization scheme used. Measurement from every mouse are separated in different columns, and represented as dissociation constants in nano-Molar units.


Class :class:`dataset` description
__________________________________

Here is a non-complete description, including only the main methods of the class. For a more complete description one can refer directly to the `code <dataset_file_>`_.

.. currentmodule:: am_sim

.. class:: dataset

   Class that contains experimental affinity measurements and details on the immunization protocol applied.

   .. method:: __init__(nM_affinity_table, D_inj, T_delay, meas_prot)

      Class initializer. It returns the newly-created dataset and takes the following arguments:

      :param list nM_affinity_table: experimental affinity measurements in nanoMolar units. They are organized in sub-lists, one per mouse tested.
      :param list D_inj: list of injected antigen dosages in micrograms, one per injection performed.
      :param list T_delay: list of time delays in days. This represents the time delay between two injections, and the time delay before affinity measurement.
      :param str meas_prot: string representing the measurement protocol used. Either ``1d`` for measurements performed one day after boost, or ``4d`` for measurements 4 days after the last injection. More measurement protocols can be implemented if needed.


   .. method:: save(filename)

      Saves the dataset in the path specified by ``filename``. Notice that ``filename`` must have ``.csv`` format.

   .. staticmethod:: load(filename)

      Static method of :class:`dataset` class used to load the dataset from the save-file ``filename``.

   .. method:: all_en()

      Method that return the cumulative list of all the affinity measurements from all the mice. Measurements are reported as binding energies.

   .. method:: mean_en()

      Method that returns the average binding energy of the cumulative list of measurements from all the mice.


.. [#f1] All the affinity measurements must fall inside the specified :ref:`experiemental sensitivity range <exp_range>`.

.. _dataset_file: https://github.com/mmolari/affinity_maturation/blob/master/am_sim/dataset.py
