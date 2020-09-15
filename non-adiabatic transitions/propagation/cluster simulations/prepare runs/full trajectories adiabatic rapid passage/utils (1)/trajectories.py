from scipy.io import loadmat
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def load_trajectories(fname):
    """
    Load trajectory coordinates from matlab file
    """
    return loadmat(fname)['x_alive'][0]

def load_velocities(fname):
    """
    Load trajectory velocities from matlab fiile
    """
    return loadmat(fname)['v_alive'][0]

def get_time_position_velocity(pos, vel):
    """
    Reconstruct time, position and velocity arrays from
    matlab data
    """
    nrPoints = len(pos[0][0])
    nrTrajectories = len(pos)

    xTrajectory = np.empty([nrPoints, nrTrajectories])
    yTrajectory = np.empty([nrPoints, nrTrajectories])
    zTrajectory = np.empty([nrPoints, nrTrajectories])
    for idP in range(nrTrajectories):
        xTrajectory[:,idP] = pos[idP][0]
        yTrajectory[:,idP] = pos[idP][1]
        zTrajectory[:,idP] = pos[idP][2]

    xVelocity = np.empty([nrPoints, nrTrajectories])
    yVelocity = np.empty([nrPoints, nrTrajectories])
    zVelocity = np.empty([nrPoints, nrTrajectories])
    for idP in range(nrTrajectories):
        xVelocity[:,idP] = vel[idP][0]
        yVelocity[:,idP] = vel[idP][1]
        zVelocity[:,idP] = vel[idP][2]

    vTot = np.sqrt(xVelocity**2+yVelocity**2+zVelocity**2)
    distance = np.sqrt(np.diff(xTrajectory, axis = 0)**2+
                       np.diff(yTrajectory, axis = 0)**2+
                       np.diff(zTrajectory, axis = 0)**2)

    t = np.zeros(vTot.shape)
    t[1:] = np.cumsum(distance/vTot[:-1], axis = 0)
    return t, (xTrajectory, yTrajectory, zTrajectory), (xVelocity, yVelocity, zVelocity)

def make_interpolate(t,x,y,z):
    """
    Interpolate trajectories
    """
    intTrajX = interp1d(t, x, kind = 'linear')
    intTrajY = interp1d(t, y, kind = 'quadratic')
    intTrajZ = interp1d(t, z, kind = 'linear')

    # rescale trajectories to start before the lens in homogeneous field
    # and end after the lens in a homogeneous field
    rootfun = lambda t: intTrajZ(t) - (0.8052-0.15)
    tstart = brentq(rootfun, 0.0, 0.01)
    rootfun = lambda t: intTrajZ(t) - (0.8052+0.6+0.15)
    tstop = brentq(rootfun, 0.0, 0.01)

    intTrajX_rescale = lambda t: intTrajX(t+tstart)
    intTrajY_rescale = lambda t: intTrajY(t+tstart)
    intTrajZ_rescale = lambda t: intTrajZ(t+tstart)-(0.8052-0.15)-0.45
    return intTrajX_rescale, intTrajY_rescale, intTrajZ_rescale, tstop-tstart
