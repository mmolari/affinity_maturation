'''
This function implements immunization schemes. Be general or specific? Maybe
specific...

An immunization scheme is defined by 3 things:
- C_inj (list)
- T_delay (list)
- measure (how the measurement is performed, 1d vs 4d)

We are interested in knowing:
- sometimes the evolution
- during the inference only the outcome
'''


class immscheme:
    def __init__(self, D_inj, T_delay, measure, par):
        # conversion Dose -> Concentration
        pass


def one_injection_immscheme(par, C_inj):
    pass


def two_injections_immscheme(par, C_inj, T_delay):
    pass
