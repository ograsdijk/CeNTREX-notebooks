from scipy.interpolate import griddata
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from .set_fontsize import set_fontsize
from matplotlib.animation import FuncAnimation

def create_grid_trajectories(interpolateTrajectories, probabilities):
    x0 = []
    y0 = []
    for idL in range(len(interpolateTrajectories)):
        traj = interpolateTrajectories[idL]
        x0.append(traj[0](0.000))
        y0.append(traj[1](0.000))
    x0 = np.array(x0)
    y0 = np.array(y0)
    
    dx = np.min(np.abs(np.diff(x0)[np.diff(x0) != 0]))
    dy = np.min(np.abs(np.diff(y0)[np.diff(y0) != 0]))
    
    xi = np.linspace(x0.min(), x0.max(), int(np.sqrt(500))*2)
    yi = np.linspace(y0.min(), y0.max(), int(np.sqrt(500))*2)
    X, Y = np.meshgrid(xi, yi)
    Z = griddata((x0,y0), probabilities, (X,Y), method = 'linear')
    return x0,y0,xi,yi,X,Y,Z
    
def plot_grid_trajectories(interpolateTrajectories, probabilities, levels = [0.01, 0.05],
                           parameters = {'width':0, 'val': ''}):
    x0,y0,xi,yi,X,Y,Z = create_grid_trajectories(interpolateTrajectories, probabilities)
    fig, ax = plt.subplots(figsize = (12,10))
    cax = ax.pcolormesh(X,Y,Z)
    vmin, vmax = cax.get_clim()
    ax.scatter(x0, y0, c = probabilities, edgecolors = 'k', s = 50, vmin = vmin, vmax = vmax, lw = 2)
    cbar = fig.colorbar(cax)
    cbar.set_label('probability', fontsize = 15)
    cs = ax.contour(xi, yi, Z, colors = ['C3', 'k'], levels = levels, linewidths = 3)
    cbar.add_lines(cs)
    ax.set_xlabel("y [m]", fontsize = 15)
    ax.set_ylabel("z [m]", fontsize = 15)
    format_string = f'>{parameters["width"]}'
    ax.set_title(f'Exit probability; {parameters["val"]:{format_string}}', fontsize = 15);
    ax.set_xticks([-0.008, -0.004, 0, 0.004, 0.008])
    set_fontsize(ax, 15)
    set_fontsize(cbar.ax, 15)
    ax.set_aspect(1/1.35)