Run the inference on your data
==============================

In our library we include an implementation of the `Paralle Tempering <par_temp_wiki_>`_ algorithm to perform the inference of model parameters on experimental data. More details on the implementation can be found in [Molari2020]_.

Run the inference algorithm
___________________________

The inference is performed through the use of the class :class:`~am_sim.parallel_tempering`. An object of this class is created by passing to the constructor the search parameters (see the signature of the :meth:`~am_sim.parallel_tempering.__init__` method). The search is then performed by executing the method :meth:`~am_sim.parallel_tempering.search` or :meth:`~am_sim.parallel_tempering.search_parallel`. The only difference between the two methdos is that the latter is parallelized with the use of the ``multiprocessing`` library.

Here is an example of use::

    # ---------- setup the search arguments
    # load the experimental datasets
    dset_list = [ds_1, ds_2, ds_3]
    # initial values of the parameters in the search. Parameters are contained in the parameters dictionary
    par = am_sim.st_par()
    # number of layers in the parallel tempering algorithm
    n_layers = 10
    # maximum number of algorithm iterations
    T_max = 10000
    # keys of the parameters to infer, as specified in the parameters dictionary
    pars_to_mutate = ['mu_i', 'sigma_i', 'k_consumption_per_day', 'alpha_C']
    # folder in which to save the results
    save_folder = 'my_results_folder'

    # ---------- initialize inference object
    partemp_inference = am_sim.parallel_tempering(dset_list, par, n_layers, T_max, pars_to_mutate, save_folder)

    # ---------- perform the search
    partemp_inference.search()

As it runs on the terminal the iteration number is printed. At the beginning the method initializes a new folder in which results are saved, as specified by the ``save_folder`` argument. This folder must be empty or non-existent to avoid overwriting existing results. The algorithm produces many files. The two main ones are :file:`search_setup.txt`, that contains details on the search setup and the initial values of the parameters, and :file:`t_*_search_history.csv`, where at the place of ``*`` the iteration time at which the file was produced is printed. This file is progressively updated during the search every 100 iterations (see parameter ``save_every`` in :meth:`~am_sim.parallel_tempering.__init__`). This file consist of a ``.csv`` table in which each entry corresponds to a set of values of the model parameters for a particular layer and a particular time. Each entry also includes the value of the log-likelihood associate to the parameters set, the round at which the parameters set was saved and details of the parallel tempering layer that generated it.

.. _ifmain_note:

.. note::

    When using the :meth:`~am_sim.parallel_tempering.search_parallel` method, as an additional caution the main script should be preceded by the ``if __name__ == '__main__'`` clause, to make sure that the code is executed only by the main thread (see for example the :file:`parallel_tempering.py` `file <github_partemp_>`_ )::

        # do all the imports
        import am_sim

        # part to be executed by the main thread
        if __name__ == '__main__':

            # ... your code here ...

            # run parallel tempering search
            partemp_inference.search()



Recover the parameters maximum-likelihood estimate
__________________________________________________

The **maximum-likelihood value of the parameters** can be found by extracting the parameters set with maximum likelihood from the :file:`t_*_search_history.csv` file. This can for example easily be done with Pandas_::

    import pandas as pd

    # load the file as a pandas DataFrame
    df = pd.read_csv('my_results_folder/t_final_search_history.csv', index_col=0)

    # find the index of the maximum-likelihood value of the parameters
    idx = df['logl'].argmax()

    # select the table entry corresponding to this index, and convert it to a parameters dictionary
    max_likl_par = df.loc[idx].to_dict()


:class:`~am_sim.parallel_tempering` class description
_____________________________________________________

Here is a non-exhaustive description of the class, containing only the main methods. For a complete description please refer directly to the `code <partemp_file_>`_.

.. currentmodule:: am_sim

.. class:: parallel_tempering

    This class contains the implementation of the `Parallel Tempering <par_temp_wiki_>`_ algorithm, used to infer the value of model parameters from experimental data.

    .. method:: __init__(dset_list, par_i, n_layers, T_max, pars_to_mutate, save_folder, [beta_list=None, mut_strength_list=None, mut_single_list=None, save_every=100])

        Class initializer. It takes as arguments the setup parameters of the search.

        :param list dset_list: list of :class:`am_sim.dataset` objects containing the experimental measurements for various immunization schemes.
        :param dict par_i: initial value of the :ref:`parameters dictionary <par_dict>` in the search.
        :param int n_layers: number of layers of the parallel tempering algorithm.
        :param int T_max: Maximum number of iterations.
        :param list pars_to_mutate: dictionary keys of the parameters to mutate, as represented in the :ref:`parameters dictionary <par_dict>`. Notice that the search is not implemented for all parameters.
        :param save_folder: folder in which the results of the search are saved. To avoid overwriting the folder must be either empty or non-existent.
        :param beta_list: optional parameter. If specified it overwrites the list of inverse temperatures for each layer.
        :type beta_list: list or None
        :param mut_strength_list: optional argument. If specified it overwrites the list of mutation strength for each layer.
        :type mut_strength_list: list or None
        :param mut_single_list: optional argument. If specified it overwrites the list of which layer mutates one parameter at a time.
        :type mut_single_list: list or None
        :param int save_every: number of iterations between two successive updates of the results file :file:`t_*_search_history.csv`.

        .. note::

            In the current version of the library not all parameters can be inferred. The limiting factor is the implementation of the parameters variation in the search algorithm. This variation can be extended to include additional parameters by modifying the implementation of the :func:`generate_variated_par` function in the :file:`am_sim/inference_utils.py` `file <https://github.com/mmolari/affinity_maturation/blob/master/am_sim/inference_utils.py>`_.

    .. method:: search()

        This method performs the parallel tempering search. Results are saved in the specified ``save_folder`` argument passed to the class initializer. They are stored in a ``.csv`` file that is iteratively updated during the search, named :file:`t_*_search_history.csv`. This function also saves the search setup parameters in a :file:`search_setup.txt` file.

    .. method:: search_parallel()

        This method is analogous to the :meth:`~am_sim.parallel_tempering.search` method, the only difference being that it is parallelized with the use of the ``multiprocessing`` library. Because of this one should add some additional specifications when using it in a script, see `this note <ifmain_note_>`_.


.. _par_temp_wiki: https://en.wikipedia.org/wiki/Parallel_tempering

.. _Pandas: https://pandas.pydata.org/docs/getting_started/index.html

.. _partemp_file: https://github.com/mmolari/affinity_maturation/blob/master/am_sim/inference.py

.. _github_partemp: https://github.com/mmolari/affinity_maturation/blob/master/parallel_tempering.py
