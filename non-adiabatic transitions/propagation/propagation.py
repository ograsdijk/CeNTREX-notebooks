from scipy.linalg import expm
import numpy as np

def calculate_state_probabilities(phi, hamiltonian):
    probabilities = []
    for p in np.linalg.eigh(hamiltonian)[1].T:
        probabilities.append(p.conj()@phi)
    return np.abs(probabilities)**2

def calculate_state_probabilities_eigvecs(phi, eigvecs):
    probabilities = []
    for p, v in zip(phi, eigvecs):
        probabilities.append(p.conj()@v)
    probabilities = np.abs(probabilities)**2
    return probabilities


def propagate_exp_HamArray(HamArray, phi0, dt):
    phi = np.zeros([len(HamArray), len(phi0)], dtype = complex)
    phi[0] = phi0

    for idH, H in enumerate(HamArray[:-1]):
        phi[idH+1] = expm(-1j*H*dt)@phi[idH]
    return phi

def propagate_exp_adaptive(ham, phi0, traj, E, B, dtAdaptive = (1e-7, 1e-6)):
    phi = phi0
    ti = 0
    while ti <= traj[-1]:
        xi = traj[0](ti)
        yi = traj[1](ti)
        zi = traj[2](ti)
        if (zi > -0.36) & (zi < -0.260):
            dt = dtAdaptive[0]
        elif (zi < 0.36) & (zi > 0.260):
            dt = dtAdaptive[0]
        else:
            dt = dtAdaptive[1]
        phi = (expm(-1j*ham(*E(xi,yi,zi), *B(xi,yi,zi))*dt)@phi)
        ti += dt

    return phi

def propagate_exp_adaptive_save(ham, phi0, traj, E, B, dtAdaptive = (1e-7, 1e-6),
                                nstep = 10):
    phi = phi0
    tList = [0]
    phiList = [phi0]
    HList = []
    ti = 0
    idSave = 0
    while ti <= traj[-1]:
        xi = traj[0](ti)
        yi = traj[1](ti)
        zi = traj[2](ti)
        if (zi > -0.36) & (zi < -0.260):
            dt = dtAdaptive[0]
        elif (zi < 0.36) & (zi > 0.260):
            dt = dtAdaptive[0]
        else:
            dt = dtAdaptive[1]
        H = ham(*E(xi,yi,zi), *B(xi,yi,zi))
        phi = (expm(-1j*H*dt)@phi)
        ti += dt
        if idSave == nstep:
            tList.append(ti)
            phiList.append(phi)
            HList.append(H)
            idSave = 0
        idSave += 1

    tList.append(ti)
    phiList.append(phi)
    HList.append(H)
    return tList, HList, phiList

def propagate_exp_adaptive_RF(ham, phi0, traj, E, B, RFfreq, RFamp, dtAdaptive = (5e-9, 1e-6)):
    phi = phi0
    ti = 0
    f = RFfreq*2*np.pi
    while ti <= traj[-1]:
        xi = traj[0](ti)
        yi = traj[1](ti)
        zi = traj[2](ti)
        if (zi > -0.45) & (zi < -0.250):
            dt = dtAdaptive[0]
            Ev = E(xi,yi,zi)
            Ex, Ey, Ez = Ev[0], Ev[1], Ev[2]
            Ey += np.cos(f*ti)*RFamp
        elif (zi < 0.45) & (zi > 0.250):
            dt = dtAdaptive[0]
            Ev = E(xi,yi,zi)
            Ex, Ey, Ez = Ev[0], Ev[1], Ev[2]
            Ey += np.cos(f*ti)*RFamp
        else:
            Ev = E(xi,yi,zi)
            Ex, Ey, Ez = Ev[0], Ev[1], Ev[2]
            dt = dtAdaptive[1]
        
        phi = (expm(-1j*ham(Ex,Ey,Ez, *B(xi,yi,zi))*dt)@phi)
        ti += dt

    return phi

def propagate_adaptive_save(ham, traj, E, B, level, dtAdaptive = (1e-7, 1e-6),
                            nstep = 10):
    xi,yi,zi = traj[0](0), traj[1](0), traj[2](0)
    H0 = ham(*E(xi,yi,zi), *B(xi,yi,zi))
    phi0 = np.linalg.eigh(H0)[1][:,level]
    t, H, phi = propagate_exp_adaptive_save(ham, phi0, traj, E, B, dtAdaptive, nstep)
    return np.array(t), [H0]+H, phi

def propagate_adaptive(ham, traj, E, B, level, dtAdaptive = (1e-7, 1e-6)):
    xi,yi,zi = traj[0](0), traj[1](0), traj[2](0)
    H0 = ham(*E(xi,yi,zi), *B(xi,yi,zi))
    phi0 = np.linalg.eigh(H0)[1][:,level]
    phi = propagate_exp_adaptive(ham, phi0, traj, E, B, dtAdaptive)
    return H0, phi, phi0

def propagate_adaptive_RF(ham, traj, E, B, level, RFfreq, RFamp, dtAdaptive = (5e-9, 1e-6)):
    xi,yi,zi = traj[0](0), traj[1](0), traj[2](0)
    H0 = ham(*E(xi,yi,zi), *B(xi,yi,zi))
    phi0 = np.linalg.eigh(H0)[1][:,level]
    phi = propagate_exp_adaptive_RF(ham, phi0, traj, E, B, RFfreq, RFamp, dtAdaptive)
    return H0, phi, phi0

def propagate_fieldcst_HamArray(Hamiltonian, level, tmax, dt, E, B,
                                   return_full = False):
    tArray = np.arange(0,tmax+dt/2,dt)
    HamList = [None]*tArray.size
    for idT, ti in enumerate(tArray):
        HamList[idT] = Hamiltonian.create_hamiltonian(*E, *B)
    phi0 = np.linalg.eigh(HamList[0])[1][:,level]
    phi = propagate_exp_HamArray(HamList, phi0, dt)

    return tArray, HamList, phi, phi0

def propagate_trajectory_HamArray(Hamiltonian, traj, level, dt, fE, B):
    tArray = np.arange(0,traj[-1],dt)
    HamList = [None]*tArray.size
    for idT, ti in enumerate(tArray):
        xi = traj[0](ti)
        yi = traj[1](ti)
        zi = traj[2](ti)
        Exi = fE[0](xi,yi,zi)
        Eyi = fE[1](xi,yi,zi)
        Ezi = fE[2](xi,yi,zi)
        HamList[idT] = Hamiltonian.create_hamiltonian(Exi, Eyi, Ezi, *B)
    phi0 = np.linalg.eigh(HamList[0])[1][:,level]
    phi = propagate_exp_HamArray(HamList, phi0, dt)

    return tArray, HamList, phi, phi0
