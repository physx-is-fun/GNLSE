from libraries import *
from variables import *

# Defining a class for the simulation parameters
# Class for holding info about the laser params
class LASER_config:
    def __init__(self,lambda0, tau0, repetition_frequency, average_power, beam_waist):
        self.lambda0 = lambda0
        self.omega0 = 2 * np.pi * c / lambda0
        self.frequency0 = self.omega0 / 2 * np.pi
        self.k0 = self.omega0 / c
        self.tau0 = tau0
        self.repetition_frequency = repetition_frequency
        self.average_power = average_power
        self.pulse_energy = average_power/repetition_frequency
        self.peak_power = self.pulse_energy / tau0
        self.beam_waist = beam_waist
        self.area=(np.pi/4)*self.beam_waist**2
        self.peak_intensity = self.peak_power / self.area
        #self.amplitude = np.sqrt((8*self.pulse_energy)/(((np.pi*2*np.log(2))**3/2)*self.tau0*self.beam_waist**2))
        self.amplitude = np.sqrt(self.peak_power)  

# Class for holding info about the fiber
class FIBER_config:
    def __init__(self, alpha_dB_per_m, beta2, beta3, effective_mode_diameter, nonlinear_refractive_index, lambda0, B1, B2, B3, C1, C2, C3):
        self.alpha_dB_per_m = alpha_dB_per_m
        self.alpha_nepers_per_m = alpha_dB_per_m / 4.343
        self.beta2 = beta2
        self.beta3 = beta3
        self.effective_mode_diameter = effective_mode_diameter
        self.effective_mode_area=(np.pi/4)*effective_mode_diameter**2
        self.nonlinear_refractive_index = nonlinear_refractive_index
        self.omega0 = 2 * np.pi * c / lambda0
        self.gamma = (nonlinear_refractive_index * self.omega0) / (c * self.effective_mode_area)
        self.B1 = B1
        self.B2 = B2
        self.B3 = B3
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3

# Class for holding info about the compressor
class COMPRESSOR_config:
    def __init__(self, GDD_mirror, N_bouncing):
        self.GDD_mirror = GDD_mirror
        self.N_bouncing = N_bouncing
        self.GDD_compressor = GDD_mirror * N_bouncing

# Class for holding info about the simulation
class SIMULATION_config:
    def __init__(self, nt, t_max, nx, x_max, ny, y_max, nz, z_max, lambda0):
        # --- Temporal grid ---
        self.frequency0 = 2 * np.pi * c / lambda0
        self.nt = nt
        self.t_max = t_max
        self.t = np.linspace(-t_max, t_max, nt)
        self.t_fs = self.t * 1e15
        self.dt = self.t[1] - self.t[0]
        self.w = fftshift(fftfreq(nt, self.dt))* 2 * np.pi 
        self.f_PHz = self.w * 1e-15
        self.f = self.w / (2 * np.pi)
        self.f_rel = self.f + self.frequency0
        self.f_PHz_rel = self.f_PHz + self.frequency0 * 1e-15

        # --- Transverse x spatial grid ---
        self.nx = nx
        self.x_max = x_max 
        self.x = np.linspace(-x_max, x_max, nx)
        self.dx = self.x[1] - self.x[0]
        self.kx = fftfreq(nx, self.dx)* 2 * np.pi 

        # --- Transverse y spatial grid ---
        self.ny = ny
        self.y_max = y_max
        self.y = np.linspace(-y_max, y_max, ny)
        self.dy = self.y[1] - self.y[0]
        self.ky = fftfreq(ny, self.dy)* 2 * np.pi 

        # --- Longitudinal z spatial grid ---
        self.nz = nz
        self.z_max = z_max 
        self.z = np.linspace(0, z_max, nz)
        self.dz = z_max / nz
        #self.dz = 1e-10
        #self.z = np.linspace(0, self.dz*nz, nz)                                            