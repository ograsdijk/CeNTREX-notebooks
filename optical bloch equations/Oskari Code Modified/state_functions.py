import numpy as np
from matrix_element_functions import ED_ME_coupled

def find_exact_states(states_approx, H, QN):
    """
    Function for finding the closest eigenstates corresponding to states_approx

    inputs:
    states_approx = list of State objects
    H = Hamiltonian whose eigenstates are used (should be diagonal in basis QN)
    QN = List of State objects that define basis for H

    returns:
    states = eigenstates of H that are closest to states_aprox in a list
    """
    states = []
    for state_approx in states_approx:
        i = find_state_idx_from_state(H, state_approx, QN)
        states.append(QN[i])

    return states

def find_state_idx_from_state(H, reference_state, QN):
    """
    This function determines the index of the state vector most closely corresponding
    to an input state

    H = Hamiltonian whose eigenstates the input state is compared to
    refernce_state = state whose index needs to be determined
    idx = index of the eigenstate that most closely corresponds to the input
    """

    #Determine state vector of reference state
    reference_state_vec = reference_state.state_vector(QN)

    #Find eigenvectors of the given Hamiltonian
    E, V = np.linalg.eigh(H)



    overlaps = np.dot(np.conj(reference_state_vec),V)
    probabilities = overlaps*np.conj(overlaps)

    idx = np.argmax(probabilities)

    return idx

def find_closest_state(H, reference_state, QN):
    """
    Function that finds the eigenstate of the Hamiltonian H that most closely corresponds to reference state.

    inputs:
    H = hamiltonian whose eigenstates reference_state is compared to
    reference_state = state which we want find in eigenstates of H

    returns:
    state = eigenstate of H closest to reference state
    """

    #Make a state vector for reference state
    reference_state_vec = reference_state.state_vector(QN)

    #Find eigenvectors of the given Hamiltonian
    E, V = np.linalg.eigh(H)

    #Find out which of the eigenstates of H corresponds to reference state
    state_index = find_state_idx(reference_state_vec,V,n=1)

    #Find state vector of state corresponding to reference
    state_vec = V[:,state_index:state_index+1]

    #return the state
    state = matrix_to_states(state_vec,QN)[0]

    return state

def calculate_BR(excited_state, ground_states, tol = 1e-5):
    """
    Function that calculates branching ratios from the given excited state to the given ground states

    inputs:
    excited_state = state object representing the excited state that is spontaneously decaying
    ground_states = list of state objects that should span all the states to which the excited state can decay

    returns:
    BRs = list of branching ratios to each of the ground states
    """

    #Initialize container fo matrix elements between excited state and ground states
    MEs = np.zeros(len(ground_states), dtype = complex)

    #loop over ground states
    for i, ground_state in enumerate(ground_states):
        MEs[i] = ED_ME_mixed_state(ground_state.remove_small_components(tol = tol),excited_state.remove_small_components(tol = tol))

    #Calculate branching ratios
    BRs = np.abs(MEs)**2/(np.sum(np.abs(MEs)**2))

    return BRs

def ED_ME_mixed_state(bra, ket, pol_vec = np.array([1,1,1]), reduced = False):
    """
    Calculates electric dipole matrix elements between mixed states

    inputs:
    bra = state object
    ket = state object
    pol_vec = polarization vector for the light that is driving the transition (the default is useful when calculating branching ratios)

    outputs:
    ME = matrix element between the two states
    """
    ME = 0
    bra = bra.transform_to_omega_basis()
    ket = ket.transform_to_omega_basis()
    for amp_bra, basis_bra in bra.data:
        for amp_ket, basis_ket in ket.data:
            ME += amp_bra.conjugate()*amp_ket*ED_ME_coupled(basis_bra, basis_ket, pol_vec = pol_vec, rme_only = reduced)

    return ME
