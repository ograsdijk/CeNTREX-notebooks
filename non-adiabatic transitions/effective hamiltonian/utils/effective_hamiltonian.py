import itertools
import numpy as np

def index_range(J):
    return 4*(J)**2, 4*(J+1)**2

def generate_Q0a(H0, P0, J, Jmax):
    E = np.linalg.eigh(H0)[0]
    idx_min, idx_max = index_range(J)
    Q0 = np.eye(4*(Jmax+1)**2) - P0
    for j in range(Jmax+1):
        if J != j:
            idmin, idmax = index_range(j)
            Q0[idmin:idmax,idmin:idmax] /= np.mean(E[idx_min:idx_max]) - np.mean(E[idmin:idmax])
    return Q0

def compress_variables(d):
    """
    Compressing variable permutations into one matrix, e.g. xy, yx -> xy.
    """
    ax = list(d.keys())
    for axes in ax:
        for c in itertools.permutations(axes):
            c = ''.join(c)
            if c != axes:
                if not isinstance(d.get(c), type(None)):
                    d[axes] += d[c]
                    ax.remove(c)
                    del d[c]
    return d

def generate_effective_hamiltonian(H0, variables, Hpertubations, J = 2, order = 2, Jmax = 6):
    """"
    Function to generate the effective Hamiltonian in the specified J state
    
    H0            : unperturbed Hamiltonian  
    variables     : list of variables corresponding to the list of Hpertubations
    Hpertubations : list of perturbing Hamiltonians
    J             : J state to create effective Hamiltonian in
    order         : to which order to calculate the pertubation
    Jmax          : maximum J state in the unperturbed Hamiltonian
    
    See Brown-Carrington Chapter 7 a description of the derivation.
    """
    # determining the indices of the block corresponding to J
    idx_min, idx_max = index_range(J)
    
    # operator selecting only the J block matrix elements
    P0 = np.eye(4*(Jmax+1)**2)
    mask = np.zeros(P0.shape, dtype = bool)
    mask[idx_min:idx_max, idx_min:idx_max] = True
    P0[~mask] = 0
    
    # operator projecting off J block elements onto the J block
    Q0 = generate_Q0a(H0, P0, J, Jmax)
    
    # calculating the contributing U terms 
    U = dict([(0,dict([('', P0)]))])
    for o in range(1, order+1):
        U[o] = dict()
        for axis, Hpert in zip(variables, Hpertubations):
            for axes, mat in U[o-1].items():
                product = Q0@(Hpert@mat)
                if np.abs(product).sum() == 0:
                    continue
                if axis+axes in U[o]:
                    U[o][axis+axes] += product
                else:
                    U[o][axis+axes] = product
            
        for axis, Hpert in zip(variables, Hpertubations):
            for p in range(1,o):
                for axes_p, mat_p in U[p].items():
                    for axes, mat in U[o-p-1].items():
                        product = Q0@(mat_p@Hpert@mat)
                        if np.abs(product).sum() == 0:
                            continue
                        if axis+axes_p+axes in U[o]:
                            U[o][axis+axes_p+axes] -= product
                        else:
                            U[o][axis+axes_p+axes] = -product
        
        # going through the variables to compress all permutations into one  matrix, e.g. xy, yx -> xy 
        U[o] = compress_variables(U[o])
    
    # generating the effective Hamiltonian
    Heff = dict()
    for o in range(order+1):
        Heff[o] = dict()
        for axis, Hpert in zip(variables, Hpertubations):
            for axes, mat in U[o].items():
                product = P0@Hpert@mat
                if np.abs(product).sum() == 0:
                    continue
                if axis+axes in Heff[o]:
                    Heff[o][axis+axes] += (P0@Hpert@mat)[idx_min:idx_max, idx_min:idx_max]
                else:
                    Heff[o][axis+axes] = (P0@Hpert@mat)[idx_min:idx_max,idx_min:idx_max]
    
        # going through the variables to compress all permutations into one  matrix, e.g. xy, yx -> xy 
        Heff[o] = compress_variables(Heff[o])
                
    return U, Heff