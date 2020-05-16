def extract_mean_en_and_rhaff(rp_list):
    '''
    function that for a list of responder populations extract the mean binding
    energy and the high-affinity fraction, and returns them in two lists.

    Args:
    - rp_list (list): list of responders populations. They can be either
        stoch_pop or det_pop objects

    Returns:
    - mean_en_list (list of floats): list of mean binding energy for each
        population, evaluated keeping into account the experimental sensitivity
        range.
    - r_haff_list (list of floats): list of high affinity fraction for each
        population, evaluated keeping into account the experimental sensitivity
        range.
    '''
    mean_en_list = [rp.mean_en_exp() for rp in rp_list]
    r_haff_list = [rp.r_haff_exp() for rp in rp_list]
    return mean_en_list, r_haff_list
