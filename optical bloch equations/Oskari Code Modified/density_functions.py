import numpy as np

def generate_density_matrix(QN, states_pop, pops):
    """
    Function for generating the density given the list of quantum numbers that define the basis,
    a list of the states that are populated and a list of the populations in the states.

    inputs:
    QN = list of state objects that defines the basis
    states_pop = states that are populated
    pops = populations in the states

    outputs:
    rho = density matrix
    """

    #Initialize the density matrix
    rho = np.zeros((len(QN), len(QN)), dtype = complex)

    #Loop over the populated states and set diagonal elements of rho
    for state, pop in zip(states_pop, pops):
        i = QN.index(state)
        rho[i,i] = pop

    return rho
