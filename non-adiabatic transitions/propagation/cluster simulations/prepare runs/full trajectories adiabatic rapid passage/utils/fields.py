import numpy as np
from numpy.polynomial import polynomial
from scipy.interpolate import interp1d, UnivariateSpline

def load_potential(fname):
    """
    Load potential from ANSYS txt file
    """
    rfp = open(fname, 'r')
    tmp = rfp.read()
    rfp.close()
    tmp = tmp.split('\n')[2:-1]
    tmp = [[float(tis) for tis in ti.split()] for ti in tmp]

    x           = np.array([row[0] for row in tmp])
    y           = np.array([row[1] for row in tmp])
    z           = np.array([row[2] for row in tmp])
    potential   = np.array([row[3] for row in tmp])

    return x,y,z,potential

def create_grid_potential(x,y,z,data):
    """
    Create a grid from the unstructured ANSYS data
    """
    x_vals, x_idx = np.unique(x, return_inverse=True)
    y_vals, y_idx = np.unique(y, return_inverse=True)
    z_vals, z_idx = np.unique(z, return_inverse=True)

    datax_grid = np.empty([max(x_vals.shape), max(y_vals.shape), max(z_vals.shape)])
    datax_grid.fill(np.nan)
    datax_grid[x_idx,y_idx,z_idx] = data
    return (x_vals, y_vals, z_vals), datax_grid

def polyfit2d(x, y, f, deg):
    """
    Fit a 2D polynomial
    """
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f, rcond=None)[0]
    return c.reshape(deg+1)

def find_nearest(array, value):
    """
    Find index of nearest value in array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def fit_potential_z(x_grid, y_grid, z_grid, pot_grid, r, zmask):
    idx = find_nearest(x_grid, np.sqrt(r**2/2))
    idy = find_nearest(y_grid, np.sqrt(r**2/2))

    maskz = z_grid > zmask
    f = np.abs(pot_grid[idx, idy, :].flatten())[maskz]
    f -= f.min()
    f /= f.max()

    l = maskz.sum()
    z_spliced = np.zeros(l*2)
    ind_sort = np.argsort(-z_grid[maskz])
    z_spliced[:l] = -z_grid[maskz][ind_sort]
    z_spliced[-l:] = z_grid[maskz]
    f_spliced = np.zeros(l*2)
    f_spliced[:l] = f[ind_sort]
    f_spliced[-l:] = f

    intfun = UnivariateSpline(z_spliced, f_spliced, s = 0)
    return intfun

def fit_potential_xy(x, y, z, potential, r, zmask):
    mask = (np.sqrt(x**2+y**2) > r)
    potential[mask] = np.nan

    maskz = np.abs(z[np.invert(mask)] - zmask) < 1e-4
    c = polyfit2d(x[np.invert(mask)][maskz], y[np.invert(mask)][maskz],
                  potential[np.invert(mask)][maskz], [3,3])

    return c

def fit_potential_xyz(fname, rxy, rz, zmaskxy, zmaskz, scale = 1):
    x,y,z,potential = load_potential(fname)
    potential *= scale
    (x_grid, y_grid, z_grid), pot_grid = create_grid_potential(x,y,z,potential)
    intfun = fit_potential_z(x_grid,y_grid,z_grid,pot_grid,rz,zmaskz)
    c = fit_potential_xy(x,y,z,potential, rxy, zmaskxy)

    return c, intfun



def combine_quadrupole_plate(c,intfun,cP,intfunP):
    x_cder = polynomial.polyder(c, axis = 0)
    x_cPder = polynomial.polyder(cP, axis = 0)
    y_cder = polynomial.polyder(c, axis = 1)
    y_cPder = polynomial.polyder(cP, axis = 1)

    Ex = lambda x,y,z: -(intfunP(z)*polynomial.polyval2d(x,y,x_cPder)+intfun(z)*polynomial.polyval2d(x,y,x_cder))/100
    Ey = lambda x,y,z: -(intfunP(z)*polynomial.polyval2d(x,y,y_cPder)+intfun(z)*polynomial.polyval2d(x,y,y_cder))/100
    Ez = lambda x,y,z: -(intfunP.derivative()(z)*polynomial.polyval2d(x,y,cP)+intfun.derivative()(z)*polynomial.polyval2d(x,y,c))/100
    Emag = lambda x,y,z: np.sqrt(Ex(x,y,z)**2+Ey(x,y,z)**2+Ez(x,y,z)**2)
    return Ex, Ey, Ez, Emag

def combine_quadrupole_plate_factor(c, intfun, cP, intfunP, factor):
    """
    Change the homogeneous field strenght by a factor
    """
    x_cder = polynomial.polyder(c, axis = 0)
    x_cPder = polynomial.polyder(cP, axis = 0)
    y_cder = polynomial.polyder(c, axis = 1)
    y_cPder = polynomial.polyder(cP, axis = 1)

    Ex = lambda x,y,z: -(factor*(intfunP(z)*polynomial.polyval2d(x,y,x_cPder))+intfun(z)*polynomial.polyval2d(x,y,x_cder))/100
    Ey = lambda x,y,z: -(factor*(intfunP(z)*polynomial.polyval2d(x,y,y_cPder))+intfun(z)*polynomial.polyval2d(x,y,y_cder))/100
    Ez = lambda x,y,z: -(factor*(intfunP.derivative()(z)*polynomial.polyval2d(x,y,cP))+intfun.derivative()(z)*polynomial.polyval2d(x,y,c))/100
    Emag = lambda x,y,z: np.sqrt(Ex(x,y,z)**2+Ey(x,y,z)**2+Ez(x,y,z)**2)
    return Ex, Ey, Ez, Emag

def potential_plate_z_offset(fname, r, zmask, d):
    x,y,z,potential = load_potential(fname)
    (x_grid, y_grid, z_grid), pot_grid = create_grid_potential(x,y,z,potential)
    idx = find_nearest(x_grid, np.sqrt(r**2/2))
    idy = find_nearest(y_grid, np.sqrt(r**2/2))

    z_grid = z_grid - d
    maskz = z_grid > zmask
    f = np.abs(pot_grid[idx, idy, :].flatten())[maskz]
    f -= f.min()
    f /= f.max()

    l = maskz.sum()
    z_spliced = np.zeros(l*2+1)
    ind_sort = np.argsort(-z_grid[maskz])
    z_spliced[:l] = -z_grid[maskz][ind_sort]
    z_spliced[-l:] = z_grid[maskz]
    f_spliced = np.zeros(l*2+1)
    f_spliced[:l] = f[ind_sort]
    f_spliced[-l:] = f

    intfun = UnivariateSpline(z_spliced, f_spliced, s = 0)
    return intfun
