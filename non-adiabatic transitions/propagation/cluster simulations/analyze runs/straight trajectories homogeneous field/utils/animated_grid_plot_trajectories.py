from .set_fontsize import *
import matplotlib.pyplot as plt
from .grid_plot_trajectories import *
import matplotlib.animation as animation

class AnimatedGridPlot:
    def __init__(self, trajectories, grouped_data, levels = [0.01, 0.05], parameters = ''):
        self.trajectories = trajectories
        self.grouped_data = grouped_data
        self.levels       = levels
        self.parameters    = parameters
        
        self.figsize = (12,10)
        
        self.prob_max = grouped_data[3].max().max()
    
    def create_animation(self, fname, interval = 2000):
        self.init_plot(min(self.grouped_data.groups.keys()))
        anim = animation.FuncAnimation(self.fig, self.animate, frames = self.grouped_data.groups.keys(), 
                                       interval = interval, blit = True)
        anim.save(fname)
                                       
    def init_plot(self, group):
        data = self.grouped_data.get_group(group).values
        probabilities = data[:,-1]
        x0,y0,xi,yi,X,Y,Z = create_grid_trajectories(self.trajectories[data[:,1].astype(int)],
                                                     probabilities)
        
        self.fig, self.ax = plt.subplots(figsize = self.figsize)
        
        cax = self.ax.pcolormesh(X,Y,Z,vmin=0,vmax=self.prob_max)
        vmin, vmax = cax.get_clim()
        scatter = self.ax.scatter(x0, y0, c = probabilities, edgecolors = 'k', s = 50, vmin = vmin, vmax = vmax, lw = 2)
        cbar = self.fig.colorbar(cax)
        cbar.set_label('probability', fontsize = 15)
        cs = self.ax.contour(xi, yi, Z, colors = ['C3', 'k'], levels = self.levels, linewidths = 3)
        cbar.add_lines(cs)
        self.ax.set_xlabel("y [m]", fontsize = 15)
        self.ax.set_ylabel("z [m]", fontsize = 15)
        format_string = f'>{self.parameters["width"]}'
        val = eval(self.parameters['val'])
        format_string = f'>{self.parameters["width"]}'
        self.ax.set_title(f'Exit probability; {val:{format_string}}', fontsize = 15);
        self.ax.set_xticks([-0.008, -0.004, 0, 0.004, 0.008])
        set_fontsize(self.ax, 15)
        set_fontsize(cbar.ax, 15)
        self.ax.set_aspect(1/1.35)
        
        self.pmesh, self.scatter, self.cs, self.cbar = cax, scatter, cs, cbar
    
    def animate(self, group):
        data = self.grouped_data.get_group(group).values
        probabilities = data[:,-1]
        x0,y0,xi,yi,X,Y,Z = create_grid_trajectories(self.trajectories[data[:,1].astype(int)],
                                                     probabilities)
        self.pmesh = self.ax.pcolormesh(X,Y,Z,vmin=0,vmax=self.prob_max)
        vmin, vmax = self.pmesh.get_clim()
        for c in self.cs.collections:
            c.remove()
        self.cs = self.ax.contour(xi, yi, Z, colors = ['C3', 'k'], levels = self.levels, linewidths = 3)
        self.scatter.remove()
        self.scatter = self.ax.scatter(x0, y0, c = probabilities, edgecolors = 'k', 
                                       s = 50, vmin = vmin, vmax = vmax, lw = 2)
        
        val = eval(self.parameters['val'])
        format_string = f'>{self.parameters["width"]}'
        self.ax.set_title(f'Exit probability; {val:{format_string}}', fontsize = 15)
        return self.pmesh, self.scatter,  