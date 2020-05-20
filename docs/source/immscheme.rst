Simulating an immunization scheme
=================================

The function :func:`~am_sim.simulate_immscheme` is used to evaluate the predicted outcome of an immunization scheme. This function takes as arguments details of the immunization scheme, the values of the model parameters, and wether the simulation should be done according to the stochastic or deterministic model. [#f1]_ Here is an example::

     # function arguments, detailing the model parameters, the simulation kind and the immunizatio scheme used
     par = am_sim.st_par() # standard value of the model parameters
     sim_type = 'deterministic' # deterministic simulation
     D_inj = [1.,1.] # two injections of 1 microgram of antigen
     T_delay = [14,28] # injections are separated by a 2 weeks interval, and measurement is performed 4 weeks after the last injection
     meas_prot = '1d' # cell harvesting 1 day after boost injection.

     # simulate the immunization scheme and return the predicted responder population
     resp_pop = simulate_immscheme(sim_type, D_inj, T_delay, meas_prot, par)

The function has the following signature:

.. currentmodule:: am_sim

.. function:: simulate_immscheme(sim_type, D_inj, T_delay, meas_prot, par)

    This function simulates an immunization scheme and returns the predicted responders population. The simulation can be either stochastic or deterministic, depending on the value of the ``sim_type`` parameter.

    :param str sim_type: either ``stochastic`` or ``deterministic``.
    :param D_inj: list of injected antigen dosages (in micrograms), one per injection performed.
    :param T_delay: list of time delays between injections and before measurement (in days).
    :param str meas_prot: measurement protocol. Either '1d' or '4d' according to weather the measurement is performed 1 day after boost or 4 days after injection.
    :param dict par: dictionary containing the value of the parameters of the model.

    :returns: the population of responders. This is either a :class:`~am_sim.stoch_pop` or a :class:`~am_sim.det_pop` object, depending on the simulation type.

.. rubric:: Footnotes:

.. [#f1] For the difference between the stochastic and deterministic model see [Molari2020]_.
