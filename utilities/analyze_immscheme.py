def extract_mean_en_and_rhaff(rp_list):
    # function that for a list of responder populations extract the mean binding
    # energy and the high-affinity fraction, and returns them in two lists
    mean_en_list = [rp.mean_en_exp() for rp in rp_list]
    r_haff_list = [rp.r_haff_exp() for rp in rp_list]
    return mean_en_list, r_haff_list
