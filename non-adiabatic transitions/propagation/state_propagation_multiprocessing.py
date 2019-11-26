import h5py
import numpy as np
from scipy.linalg import expm
import pickle
from numpy.polynomial import polynomial

def create_hamiltonian(Ex,Ey,Ez,Bx,By,Bz, Hff_m, HSx_m, HSy_m, HSz_m, HZx_m, HZy_m, HZz_m):
    HamE =  Hff_m + \
            Ex*HSx_m  + Ey*HSy_m + Ez*HSz_m + \
            Bx*HZx_m  + By*HZy_m + Bz*HZz_m
    return HamE

def propagate_multi(args):
    start_level, B, Hterms, traj_range, dt, fname_traj, fname_pot, q = args
    
    with open(fname_traj, 'rb') as rfp:
        interpolateTrajectories = pickle.load(rfp)
    with open(fname_pot, 'rb') as rfp:
        interpolate_list = pickle.load(rfp)
    c, intfun, cP, intfunP = interpolate_list
    
    x_cder = polynomial.polyder(c, axis = 0)
    x_cPder = polynomial.polyder(cP, axis = 0)
    y_cder = polynomial.polyder(c, axis = 1)
    y_cPder = polynomial.polyder(cP, axis = 1)

    intEx = lambda x,y,z: (intfunP(z)*polynomial.polyval2d(x,y,x_cPder)+intfun(z)*polynomial.polyval2d(x,y,x_cder))/100
    intEy = lambda x,y,z: (intfunP(z)*polynomial.polyval2d(x,y,y_cPder)+intfun(z)*polynomial.polyval2d(x,y,y_cder))/100
    intEz = lambda x,y,z: (intfunP.derivative()(z)*polynomial.polyval2d(x,y,cP)+intfun.derivative()(z)*polynomial.polyval2d(x,y,c))/100
    
    phi0_list = []
    phi_list = []
    H0_list = []
    H_list = []
    idT_list = []
    dt0 = dt[0]
    dt1 = dt[1]
    for idT in range(*traj_range):
        traj = interpolateTrajectories[idT]
        xi = traj[0](0)
        yi = traj[1](0)
        zi = traj[2](0)
        xs, ys, zs = xi, yi, zi
        Exi = intEx(xi,yi,zi)
        Eyi = intEy(xi,yi,zi)
        Ezi = intEz(xi,yi,zi)
        H0 = create_hamiltonian(intEx(xi,yi,zi), intEy(xi,yi,zi), intEz(xi,yi,zi), *B, *Hterms)
        phi0 = np.linalg.eigh(H0)[1][:,start_level]

        ti = 0
        H = 0
        phi = phi0
        while ti <= traj[-1]:
            xi = traj[0](ti)
            yi = traj[1](ti)
            zi = traj[2](ti)
            if (zi > -0.36) & (zi < -0.250):
                dt = dt0
            elif (zi < 0.36) & (zi > 0.250):
                dt = dt0
            else:
                dt = dt1
            Exi = intEx(xi,yi,zi)
            Eyi = intEy(xi,yi,zi)
            Ezi = intEz(xi,yi,zi)
            H = create_hamiltonian(Exi, Eyi, Ezi, *B, *Hterms)
            phi = expm(-1j*H*dt)@phi
            ti += dt
        
        idT_list.append(idT)
        phi0_list.append(phi0)
        phi_list.append(phi)
        H0_list.append(H0)
        H_list.append(H)
        q.put(1)
    q.put(None)
    return idT_list, phi0_list, phi_list, H0_list, H_list