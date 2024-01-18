

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

## Constant variables
A = 6.022e23 # Avogadro's Number (float)
M = 66.5 # kg/mol BSA (float)
time = np.arange(0,20e-9,1e-12)+1e-12 # "100 ns" time in 1 ps steps, don't start at zero (numpy array)
temp = 298 # (int)
kB = 1.38e-23 # J/K = kg*m^2/s^2/K (float)
R = 8.314 # J/mol/K = kg*m^2/s^2/K/mol (float)

D = 7.1e-6 # 7.1 nm (float)


ff_k = 1e-6 # (float)
ff_n = -3 # (int)
sigma = np.sqrt((R*temp)/M)# m/s = um/us = nm/ns (float)


class Spectrophotometer:
    def __init__(self, N=20, D=7.1e-6, width=10, height=1, edge_ff_k=2e1, edge_ff_n=-6, particle_ff_k=1e1, particle_ff_n=-3, laser_steps=1000, eps=1e-7, figure=True):
        self.width, self.height = width, height
        self.lower_bound = np.array([[[0, -0.5*height]]])
        self.upper_bound = np.array([[[width, 0.5*height]]])

        self.n_laser_steps = laser_steps
        self.laser_steps = np.arange(0, width, self.n_laser_steps)
        self.laser_hits = np.zeros(self.n_laser_steps)
        
        assert int(N) > 0, 'Must include at least one particle'
        self.N = int(N) # number of particles
        self.D = float(D) # diameter of the particles, in mm. 7.1e-6 = 7.1 nm
        self.coordinates = np.random.uniform(self.lower_bound, self.upper_bound, size=(1,self.N,2))
        self.orientation = np.random.normal(size=(1,self.N,2)) # direction, will need to be normalized
        
        self.edge_ff_k = edge_ff_k
        self.edge_ff_n = edge_ff_n
        self.particle_ff_k = particle_ff_k
        self.particle_ff_n = particle_ff_n
        
        self.time = [0.0] # observed time steps
        self.history = [np.array(self.coordinates)] # Record of past coordinate positions 
        self.visible = [np.all(np.abs(self.coordinates[:,1])-self.D>0)]
        
        
        self.eps = eps # Small number used for numerical stability near zero
        self.absorbance = (~self.visible[0])*self.eps

        self.measure_absorbance()
        if figure:
            self.generate_figure()

    def normalize_orientation(self):
        self.orientation = self.orientation/(np.sqrt(np.sum(self.orientation**2, axis=-1, keepdims=True))+self.eps)
        
    def apply_edge_ff(self):
        lower_distances = self.edge_ff_k * (self.coordinates - self.lower_bound)
        upper_distances = self.edge_ff_k * (self.upper_bound - self.coordinates)
        self.orientation += np.power(lower_distances,self.edge_ff_n) - np.power(upper_distances,self.edge_ff_n)
        
    def apply_particle_ff(self, dt):
        x, y = np.split(self.coordinates, 2, axis=-1) # both arrays now (1,N,1)
        delta = np.concatenate((x-x.transpose((0,2,1)), y-y.transpose((0,2,1))), axis=0).transpose((1,2,0))
        d = np.sum(delta**2, axis=-1, keepdims=True)
        weight = np.power(self.particle_ff_k*d+self.eps, self.particle_ff_n)
        self.orientation += np.sum(weight*delta*dt, axis=1)        
        
    def measure_absorbance(self, total_time=10, time_step=0.001, use_particle_ff=False): # Times in ns, 0.001 = 1 ps

        # Evolve the system
        current_time = self.time[-1]
        time_steps = np.arange(0,total_time,time_step)+time_step
        
        velocities = stats.chi.rvs(df=2, scale=sigma, size=(len(time_steps),self.N,1))
        for t, v in zip(time_steps, velocities):
            # v is a (N,1) matrix ready to be multiplied by the normalized orientation
            dt = t - current_time
            dd = v*dt # v*dt is in m, so multiply by 1e3 to get to mm
            
            pre_sign = np.sign(self.coordinates)
            self.apply_edge_ff()
            if use_particle_ff:
                self.apply_particle_ff(dt)
            self.normalize_orientation()
            self.coordinates += self.orientation*dd
            post_sign = np.sign(self.coordinates)
            if np.any(pre_sign != post_sign) or np.any(np.abs(self.coordinates[:,1])-self.D<0):
                self.visible.append(False)
                (post_sign-pre_sign != 0)*self.coordinates
            else:
                self.visible.append(True)
            self.history.append(np.array(self.coordinates))
            self.time.append(t)
            current_time = t

        # Calculate absorbance
        transmittance = np.sum(self.visible)/len(self.visible)
        self.absorbance = -np.log10(transmittance)
            
        return self.absorbance
        

    def slice_history(self, frame, tail_n):
        history_slice = np.concatenate(self.history[max([0,frame-tail_n]):frame], axis=0)
        return np.transpose(history_slice, (1,2,0))
        

    def generate_figure(self):
        fig, ax = plt.subplots(figsize=(self.width,self.height), dpi=150)
        ax.axis('off')
        border = Rectangle(self.lower_bound[0][0], self.width, self.height, facecolor='#FFFFFF', #'#81D0F2', 
                           edgecolor='#444444', linewidth=2)
        ax.add_patch(border)
        
        pnt_s = 6
        pnt_color = '#994444'
        pnt_alpha = 0.8

        tail_n = 10000
        tail_style = 'dotted'
        tail_width = 1.0
        tail_color = '#444444'
        tail_alpha = 0.05

        ax.scatter([-0.35],[0], c='#cc4444', s=400)
        ax.plot([-0.2,self.width+0.2], [0]*2, c='#cc4444', linewidth=2)
        scatter = ax.scatter( *self.history[-1].T, s=pnt_s, c=pnt_color, alpha=pnt_alpha )
        lines = [ax.plot(*z, linestyle=tail_style, linewidth=tail_width, color=tail_color, alpha=tail_alpha) 
                         for z in self.slice_history(len(self.time), tail_n)]
        ax.plot([self.width+0.2]*2, [0.8*self.lower_bound[0,0,-1],0.8*self.upper_bound[0,0,-1]], c='#444444', linewidth=3)
        ax.text(self.width+0.25, 0, f'Abs = {self.absorbance:.3}', fontsize=10)
        ax.set_title('CyberSpec 6765', style='italic', weight='bold', size=12)
        
        
        return fig
        
        
        