import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from .set_fontsize import set_fontsize
from scipy.interpolate import griddata


def create_grid(x,y,z):
    dx = np.min(np.abs(np.diff(x)[np.diff(x) != 0]))
    dy = np.min(np.abs(np.diff(y)[np.diff(y) != 0]))
    
    xi = np.linspace(x.min(), x.max(), len(np.unique(x)))
    yi = np.linspace(y.min(), y.max(), len(np.unique(y)))
    
    X, Y = np.meshgrid(xi, yi)
    Z = griddata((x,y), z, (X,Y), method = 'linear')
    return xi,yi,X,Y,Z

def plot_pmesh(x,y,z, levels = [], parameters = {'width':0, 'val':''},
              xlabel = '', ylabel = '', scale = 'linear', shading = 'flat',
              limits = None):
       
    if limits:
        (x1,y1), (x2,y2) = limits
        maskx = (x > x1) & (x < x2)
        masky = (y > y1) & (y < y2)
        mask = maskx & masky
        x = x[mask]
        y = y[mask]
        z = z[mask]
        
    xi,yi,X,Y,Z = create_grid(x,y,z)
        
    fig, ax = plt.subplots(figsize = (12,10))
    if scale == 'log':
        cax = ax.pcolormesh(X,Y,Z, shading = shading, 
                           norm=LogNorm(vmin=Z.min(), vmax=Z.max()))
    else:
        cax = ax.pcolormesh(X,Y,Z, shading = shading)
    vmin, vmax = cax.get_clim()
    cbar = fig.colorbar(cax)
    cbar.set_label('probability', fontsize = 15)
    if len(levels) > 0:
        cs = ax.contour(xi, yi, Z, colors = ['C3', 'k'], levels = levels, linewidths = 3)
        cbar.add_lines(cs)
    ax.set_xlabel(xlabel, fontsize = 15)
    ax.set_ylabel(ylabel, fontsize = 15)
    format_string = f'>{parameters["width"]}'
    ax.set_title(f'Exit probability; {parameters["val"]:{format_string}}')
    set_fontsize(ax, 15)
    set_fontsize(cbar.ax, 15)
    ax.set_aspect(xi.ptp()/yi.ptp())