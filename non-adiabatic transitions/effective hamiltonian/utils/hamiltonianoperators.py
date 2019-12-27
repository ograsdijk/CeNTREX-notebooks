import numpy as np
from numpy import sqrt
from .TlF import *

def Hrot(psi):
    return Brot * J2(psi)

def Hc1(psi):
    return c1 * ( com(I1z,Jz,psi) + .5*(com(I1p,Jm,psi)+com(I1m,Jp,psi)) )

def Hc2(psi):
    return c2 * ( com(I2z,Jz,psi) + .5*(com(I2p,Jm,psi)+com(I2m,Jp,psi)) )

def Hc4(psi):
    return c4 * ( com(I1z,I2z,psi) + .5*(com(I1p,I2m,psi)+com(I1m,I2p,psi)) )

def Hc3a(psi):
    return 15*c3/c1/c2 * com(Hc1,Hc2,psi) / ((2*psi.J+3)*(2*psi.J-1))

def Hc3b(psi):
    return 15*c3/c1/c2 * com(Hc2,Hc1,psi) / ((2*psi.J+3)*(2*psi.J-1))

def Hc3c(psi):
    return -10*c3/c4/Brot * com(Hc4,Hrot,psi) / ((2*psi.J+3)*(2*psi.J-1))

def Hff(psi):
    return Hrot(psi) + Hc1(psi) + Hc2(psi) + Hc3a(psi) + Hc3b(psi) \
            + Hc3c(psi) + Hc4(psi)

def HZx(psi):
    if psi.J != 0:
        return -mu_J/psi.J*Jx(psi) - mu_Tl/psi.I1*I1x(psi) - mu_F/psi.I2*I2x(psi)
    else:
        return -mu_Tl/psi.I1*I1x(psi) - mu_F/psi.I2*I2x(psi)

def HZy(psi):
    if psi.J != 0:
        return -mu_J/psi.J*Jy(psi) - mu_Tl/psi.I1*I1y(psi) - mu_F/psi.I2*I2y(psi)
    else:
        return -mu_Tl/psi.I1*I1y(psi) - mu_F/psi.I2*I2y(psi)

def HZz(psi):
    if psi.J != 0:
        return -mu_J/psi.J*Jz(psi) - mu_Tl/psi.I1*I1z(psi) - mu_F/psi.I2*I2z(psi)
    else:
        return -mu_Tl/psi.I1*I1z(psi) - mu_F/psi.I2*I2z(psi)

def R10(psi):
    amp1 = sqrt((psi.J-psi.mJ)*(psi.J+psi.mJ)/(8*psi.J**2-2))
    ket1 = BasisState(psi.J-1, psi.mJ, psi.I1, psi.m1, psi.I2, psi.m2)
    amp2 = sqrt((psi.J-psi.mJ+1)*(psi.J+psi.mJ+1)/(6+8*psi.J*(psi.J+2)))
    ket2 = BasisState(psi.J+1, psi.mJ, psi.I1, psi.m1, psi.I2, psi.m2)
    return State([(amp1,ket1),(amp2,ket2)])

def R1m(psi):
    amp1 = -.5*sqrt((psi.J+psi.mJ)*(psi.J+psi.mJ-1)/(4*psi.J**2-1))
    ket1 = BasisState(psi.J-1, psi.mJ-1, psi.I1, psi.m1, psi.I2, psi.m2)
    amp2 = .5*sqrt((psi.J-psi.mJ+1)*(psi.J-psi.mJ+2)/(3+4*psi.J*(psi.J+2)))
    ket2 = BasisState(psi.J+1, psi.mJ-1, psi.I1, psi.m1, psi.I2, psi.m2)
    return State([(amp1,ket1),(amp2,ket2)])

def R1p(psi):
    amp1 = -.5*sqrt((psi.J-psi.mJ)*(psi.J-psi.mJ-1)/(4*psi.J**2-1))
    ket1 = BasisState(psi.J-1, psi.mJ+1, psi.I1, psi.m1, psi.I2, psi.m2)
    amp2 = .5*sqrt((psi.J+psi.mJ+1)*(psi.J+psi.mJ+2)/(3+4*psi.J*(psi.J+2)))
    ket2 = BasisState(psi.J+1, psi.mJ+1, psi.I1, psi.m1, psi.I2, psi.m2)
    return State([(amp1,ket1),(amp2,ket2)])

def HSx(psi):
    return -D_TlF * ( R1m(psi) - R1p(psi) )

def HSy(psi):
    return -D_TlF * 1j * ( R1m(psi) + R1p(psi) )

def HSz(psi):
    return -D_TlF * sqrt(2)*R10(psi)

def HI1R(psi):
    return com(I1z,R10,psi) + .5*(com(I1p,R1m,psi)+com(I1m,R1p,psi))

def HI2R(psi):
    return com(I2z,R10,psi) + .5*(com(I2p,R1m,psi)+com(I2m,R1p,psi))

def Hc3_alt(psi):
    return 5*c3/c4*Hc4(psi) - 15*c3/2*(com(HI1R,HI2R,psi)+com(HI2R,HI1R,psi))

def Hff_alt(psi):
    return Hrot(psi) + Hc1(psi) + Hc2(psi) + Hc3_alt(psi) + Hc4(psi)

def HMatElems(H, QN):
    result = np.empty((len(QN),len(QN)), dtype=complex)
    for i,a in enumerate(QN):
        for j,b in enumerate(QN):
            result[i,j] = (1*a)@H(b)
    return result

def eigenstates(Ex,Ey,Ez,Bx,By,Bz,epsilon=1e-6):
    # diagonalize the Hamiltonian
    H = Hff_m \
        + Ex*HSx_m + Ey*HSy_m + Ez*HSz_m \
        + Bx*HZx_m + By*HZy_m + Bz*HZz_m
    eigvals,eigvecs = np.linalg.eigh(H)

    # find the quantum numbers of the largest-|amplitude| components
    states = []
    for eigvec in eigvecs.T:
        # normalize the largest |amplitude| to 1
        eigvec = eigvec / np.max(np.abs(eigvec))
        # find indices of the largest-|amplitude| components
        major = np.abs(eigvec) > epsilon

        # collect the major components into a State
        eigenstate = State()
        for amp,psi in zip(eigvec[major], QN[major]):
            eigenstate += amp * psi

        # sort the components by decreasing |amplitude|
        amps = np.array(eigenstate.data).T[0]
        cpts = np.array(eigenstate.data).T[1]
        cpts = cpts[np.argsort(np.abs(amps))]
        amps = amps[np.argsort(np.abs(amps))]
        sorted_state = State( data=np.array((amps,cpts)).T )
        states.append(sorted_state)

    return eigvals, np.array(states)

class Hamiltonian:
    def __init__(self, Jmax, I_Tl, I_F):
        self.Jmax = Jmax
        self.I_Tl = I_Tl
        self.I_F = I_F
        self.calculate_HMatElems()

    def calculate_HMatElems(self):
        self.QN = np.array([BasisState(J,mJ,self.I_Tl,m1,self.I_F,m2)
                      for J in range(self.Jmax+1)
                      for mJ in range(-J,J+1)
                      for m1 in np.arange(-self.I_Tl,self.I_Tl+1)
                      for m2 in np.arange(-self.I_F,self.I_F+1)])
        self.Hff_m = HMatElems(Hff, self.QN)
        self.HSx_m = HMatElems(HSx, self.QN)
        self.HSy_m = HMatElems(HSy, self.QN)
        self.HSz_m = HMatElems(HSz, self.QN)
        self.HZx_m = HMatElems(HZx, self.QN)
        self.HZy_m = HMatElems(HZy, self.QN)
        self.HZz_m = HMatElems(HZz, self.QN)

    def create_hamiltonian(self, Ex,Ey,Ez,Bx,By,Bz):
        multi = lambda a, b: np.einsum('i,jk->ijk', a, b)
        if type(Ex) == np.ndarray:
            HamE =  multi(np.ones(len(Ex)),self.Hff_m) + \
                    multi(Ex,self.HSx_m) + multi(Ey,self.HSy_m) + \
                    multi(Ez,self.HSz_m) + multi(Bx,self.HZx_m) + \
                    multi(By,self.HZy_m) + multi(Bz,self.HZz_m)
        else:
            HamE =  self.Hff_m + \
                    Ex*self.HSx_m  + Ey*self.HSy_m + Ez*self.HSz_m + \
                    Bx*self.HZx_m  + By*self.HZy_m + Bz*self.HZz_m
        return HamE

    def spectrum(self,  Ex, Ey, Ez, Bx, By, Bz):
        H = self.create_hamiltonian( Ex, Ey, Ez, Bx, By, Bz)
        return np.linalg.eigh(H)

    def eigenstates(self, Ex, Ey, Ez, Bx, By, Bz, epsilon=1e-6):
        H = self.create_hamiltonian( Ex, Ey, Ez, Bx, By, Bz)
        eigenvalues, eigenvectors = np.linalg.eigh(H)
        states = []
        for eigenvectors in eigenvectors.T:
            # normalize the largest |amplitude| to 1
            eigenvectors = eigenvectors / np.max(np.abs(eigenvectors))
            # find indices of the largest-|amplitude| components
            major = np.abs(eigenvectors) > epsilon

            # collect the major components into a State
            eigenstate = State()
            for amp,psi in zip(eigenvectors[major], self.QN[major]):
                eigenstate += amp * psi

            # sort the components by decreasing |amplitude|
            amps = np.array(eigenstate.data).T[0]
            cpts = np.array(eigenstate.data).T[1]
            cpts = cpts[np.argsort(np.abs(amps))]
            amps = amps[np.argsort(np.abs(amps))]
            sorted_state = State( data=np.array((amps,cpts)).T )
            states.append(sorted_state)
        return eigenvalues, np.array(states)

    def largest_eigenstate(self, Ex,Ey,Ez,Bx,By,Bz,J=1,epsilon=.95):
        energies, states = self.eigenstates(Ex,Ey,Ez,Bx,By,Bz,epsilon=epsilon)
        print("E [kHz]\t J, mJ, m1, m2\n-----------------------")
        for i in reversed(np.arange((2*self.I_Tl+1)*(2*self.I_F+1)*J**2,(2*self.I_Tl+1)*(2*self.I_F+1)*(J+1)*(J+1))):
            print("%+0.2f" % ((energies[int(i)])/1e3), end='\t ')
            states[int(i)][0][1].print_quantum_numbers()

    def major_eigenstates(self, Ex,Ey,Ez,Bx,By,Bz,J=1,epsilon=.95):
        energies, states = self.eigenstates(Ex,Ey,Ez,Bx,By,Bz,epsilon=epsilon)
        for i in reversed(np.arange((2*self.I_Tl+1)*(2*self.I_F+1)*J**2,(2*self.I_Tl+1)*(2*self.I_F+1)*(J+1)*(J+1))):
            print("E = %+0.16e" % ((energies[int(i)])/1e3),"kHz",end='\n')
            for amp,psi in states[int(i)]:
                print("\t%+0.3f"%np.real(amp),"%+0.3fi"%np.imag(amp), end=' ')
                psi.print_quantum_numbers()

    def level_eigenstates(self, Ex,Ey,Ez,Bx,By,Bz, level,epsilon=.95):
        energies, states = self.eigenstates(Ex,Ey,Ez,Bx,By,Bz,epsilon=epsilon)
        for l in level:
            print('-'*23)
            print('Level = {0}'.format(l))
            print("E = %+0.16e" % ((energies[int(l)])/1e3),"kHz",end='\n')
            for amp,psi in states[int(l)]:
                print("\t%+0.3f"%np.real(amp),"%+0.3fi"%np.imag(amp), end=' ')
                psi.print_quantum_numbers()
