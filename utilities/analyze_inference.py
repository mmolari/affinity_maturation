def best_par_in_df(search_df):
    '''
    This function takes as input a pandas dataframe containing the inference
    results (saved in a file having 'search_history.csv' ending) and returns
    the highest-likelihood parameters set that was found during the search.

    Args:
    - search_df (pandas DataFrame object): dataframe containing the results of
        the inference procedure, produced by the 'parallel_tempering' class.

    Returns:
    - best_par: model parameters dictionary containing the highest-likelihood
        parameters set found during the inference.
    '''
    best_idx = search_df['logl'].argmax()
    best_par = search_df.loc[best_idx].to_dict()
    return best_par
