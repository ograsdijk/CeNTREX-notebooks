import numpy as np
import scipy.constants as constants
from state_functions import ED_ME_mixed_state
from matrix_element_functions import ED_ME_coupled

def optical_coupling_matrix(QN, ground_states, excited_states, pol_vec = np.array([0,0,1]), reduced = False):
    """
    Function that generates the optical coupling matrix for given ground and excited states

    inputs:
    QN = list of states that defines the basis for the calculation
    ground_states = list of ground states that are coupled to the excited states (i.e. laser is resonant)
    excited_states = list of excited states that are coupled to the ground states

    outputs:
    H = coupling matrix
    """

    #Initialize the coupling matrix
    H = np.zeros((len(QN),len(QN)), dtype = complex)

    #Start looping over ground and excited states
    for ground_state in ground_states:
        i = QN.index(ground_state)
        for excited_state in excited_states:
            j = QN.index(excited_state)

            #Calculate matrix element and add it to the Hamiltonian
            H[i,j] = ED_ME_mixed_state(ground_state, excited_state, pol_vec = pol_vec, reduced = reduced)

    #Make H hermitian
    H = H + H.conj().T

    return H

def calculate_power_needed(Omega, ME, fwhm = 0.0254, D_TlF = 13373921.308037223):
    """
    Function to calculate the power required to get peak Rabi rate Omega
    for a transition with given matrix element with a gaussian spatial profile
    """

    #Calculate the electric field required (in V/m)
    E =  Omega/(ME*D_TlF) * 100

    #Convert E to peak intensity
    c = constants.c
    epsilon_0 = constants.epsilon_0
    I = 1/2 * c * epsilon_0 * E**2

    #Convert FWHM to standard deviation
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))

    #Convert power to amplitude of the Gaussian
    P = I * (2*np.pi *sigma**2)

    return P

def laser_field(x, z0 = 0, fwhm = 0.0254, power = 1):
    """
    Function that calculates the electric field at position x due to a
    laser with a Gaussian intensity profile defined by its width (fwhm)
    and total power.

    inputs:
    x = position where electric field is to be evaluated (meters)
    z0 = position of the center of the microwave beam
    fwhm = full-width-half-maximum of the intensity profile
    power = output power in watts

    returns:
    E = magnitude of electric field at x
    """

    #Convert FWHM to standard deviation
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))

    #Convert power to amplitude of the Gaussian
    I0 = power/(2*np.pi *sigma**2)

    #Get the value of z where the field needs to be evaluated
    z = x[2]

    #Calculate intensity at (0,0,z)
    I_z = I0 * np.exp(-1/2*((z-z0)/sigma)**2)

    #Calculate electric field from intensity (in V/m)
    c = constants.c
    epsilon_0 = constants.epsilon_0
    E = np.sqrt(2*I_z/(c*epsilon_0))

    #Return electric field in V/cm
    return E/100

def sidebands(beta,omega_sb = 2*np.pi*1.6e6, n_cut = 15):
    """
    This function provides the modified amplitude for a laser that is being phase modulated by an EOM

    inputs:
    beta = modulation depth
    omega_sb = frequency of modulation
    n_cut = order at which the sideband expansion is cut off

    outputs:
    sideband_amps(t) = lambda function of time which gives the amplitudes of the sidebands
    """

    sideband_amps = lambda t: (jv(0, beta) + np.sum([jv(k,beta)*(np.exp(1j*k*omega_sb*t) + (-1)**k*np.exp(-1j*k*omega_sb*t))
                                                    for k in range(1, n_cut)]))

    return sideband_amps
