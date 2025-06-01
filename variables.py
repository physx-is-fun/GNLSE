from libraries import *

# --- Physical constants ---
c = 3e8                                                  # Speed of light in vacuum [m/s]
n0 = 1.0                                                 # Refractive index of air
epsilon_0 = 8.854e-12                                    # Dielectric constant [F/m]

# --- Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP ---
lambda0 = 800e-9                                         # Pulse central wavelengt [m]
tau0 = 20e-15                                            # Pulse duration [s]
repetition_frequency = 85*1e6                            # Repetition frequency [Hz]
average_power = 600*1e-3                                 # Average power [W]

# --- Initialize 780HP single mode fiber from Thorlabs ---
alpha_dB_per_m=0.2e-3                                    # Power attenuation coeff in decibel per m. Usual value at 1550 nm is 0.2 dB/km
effective_mode_diameter=5e-6                             # Effective mode diameter [m] for 780HP single mode fiber @ 850 nm from THORLABS
nonlinear_refractive_index=2.7e-20                       # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
beam_waist = 100e-6                                      # Beam waist [m]
beta2=36.16                                              # Dispersion in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2*=(1e-27)                                           # Convert fs^2 to s^2 so everything is in SI units of fused silica @ 800nm
beta3=27.47                                              # Dispersion in fs^3/mm (units typically used when referring to beta3) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta3*=(1e-42)                                           # Convert f3^2 to s^3 and mm to m so everything is in SI units of fused silica @ 800nm
# Sellmeier coefficients for fused silica (Malitson, 1965)
B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
C1, C2, C3 = 0.0684043**2, 0.1162414**2, 9.896161**2  # [µm²]

# --- Initialize UMC10-15FS - Ø1" Dispersion-Compensating Mirror from Thorlabs --- https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=3746&pn=UMC10-15FS
GDD_mirror = -54e-30                                     # Group delay dispersion [s^2], -1.5 mm of Fused Silica
N_bouncing = 1                                           # Bouncing number

# --- Temporal grid ---
nt = 2**8                                                # Number of points in the temporal grid
t_max = 200e-15                                          # Maximum value in the time axis [s]

# --- Longitudinal z spatial grid ---
nz = 2**8                                                # Number of points in the z grid
z_max = 1e-3                                             # Maximum value in the z axis [m] 

# --- Transverse x spatial grid ---
nx = 2**6                                                # Number of points in the x grid
x_max = 1e-3                                             # Maximum value in the x axis [m] 

# --- Transverse y spatial grid ---                       
ny = 2**6                                                # Number of points in the y axis
y_max = 1e-3                                             # Number of points in the y axis