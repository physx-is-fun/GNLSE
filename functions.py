from variables import *
from libraries import *
from classes import *

'''
def getIntensity(amplitude):
    return (1/2)*n0*epsilon_0*c*np.abs(amplitude)**2
'''

def getIntensity(amplitude):
    return np.abs(amplitude)**2

'''
def getEnergy(A, simulation:SIMULATION_config, fiber:FIBER_config):
    I = getIntensity(A)
    I_t = np.trapz(I, dx=simulation.dt, axis=2)
    I_ty = np.trapz(I_t, dx=simulation.dy, axis=1)
    energy = np.trapz(I_ty, dx=simulation.dx, axis=0)
    return energy
'''
    
def getPhotonNumber(A, simulation:SIMULATION_config, fiber:FIBER_config, laser:LASER_config):
    """
    Calculate photon number using:
        N = ∭ |A(x, y, t)|^2 / (ħω) dx dy dt
    """
    I = getIntensity(A)  # shape: (nx, ny, nt)
    omega = simulation.w + laser.omega0  # shape: (nt,)
    photon_energy = hbar * omega  # shape: (nt,)

    # Broadcast photon energy to match I shape
    photon_energy = photon_energy[None, None, :]  # shape: (1, 1, nt)

    I_over_hw = I / photon_energy  # shape: (nx, ny, nt)

    # Integrate over t
    int_t = np.trapz(I_over_hw, simulation.t, axis=2)  # shape: (nx, ny)

    # Integrate over y
    int_y = np.trapz(int_t, simulation.y, axis=1)  # shape: (nx,)

    # Integrate over x
    photon_number = np.trapz(int_y, simulation.x, axis=0)  # scalar

    return photon_number


def GaussianPulse(time,duration_FWHM, X, Y, beam_waist_FWHM, amplitude):
    temporal_profile = np.exp(-4*np.log(2)*(time/duration_FWHM)**2)
    spatial_profile = np.exp(-4*np.log(2)*(X/beam_waist_FWHM)**2) * np.exp(-4*np.log(2)*(Y/beam_waist_FWHM)**2)
    return amplitude * temporal_profile * spatial_profile

# Getting the spectrum based on a given pulse
def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def getPulseFromSpectrum(time,spectrum_aplitude):
    dt=time[1]-time[0]
    pulse=ifft(ifftshift(spectrum_aplitude))/dt
    return pulse

# Getting FWHM based on a given pulse
# Find the FWHM of the frequency/time domain of the signal
def FWHM(X, Y):
    deltax = X[1] - X[0]
    half_max = max(Y) / 2.
    l = np.where(Y > half_max, 1, 0)
    return np.sum(l) * deltax

def estimate_stable_dz(gamma, P0, lambda0, n0, beta2=None, omega_max=None, C_nl=0.1, C_phi=0.1):
    """
    Estimate physically stable step size Δz for SSFM in NLSE.

    Parameters:
        gamma     : Nonlinear coefficient [1/(W·m)]
        P0        : Peak power [W]
        lambda0   : Central wavelength [m]
        n0        : Linear refractive index
        beta2     : (Optional) GVD coefficient [s^2/m]
        omega_max : (Optional) max angular frequency content [rad/s]
        C_nl      : Safety factor for nonlinear phase shift
        C_phi     : Safety factor for spatial phase shift

    Returns:
        dz_suggested : Estimated stable step size [m]
    """
    dz_nl = C_nl / (gamma * P0)
    dz_phi = C_phi * lambda0 / (2 * np.pi * n0)

    if beta2 is not None and omega_max is not None:
        dz_disp = 1 / (np.abs(beta2) * omega_max**2)
        return min(dz_nl, dz_phi, dz_disp)

    return min(dz_nl, dz_phi)

def estimate_grid_resolution(tau0, w0, omega_max=None):
    """
    Estimate temporal and spatial grid resolutions.

    Parameters:
        tau0      : Pulse duration (e.g., FWHM) [s]
        w0        : Beam waist (radius) [m]
        omega_max : Maximum angular frequency content (optional) [rad/s]

    Returns:
        dt_est    : Temporal grid spacing [s]
        dx_est    : Spatial grid spacing [m]
        dy_est    : Spatial grid spacing [m]
    """

    # Temporal resolution based on pulse duration
    dt_est = tau0 / 10

    # If max spectral content is known, enforce Nyquist
    if omega_max is not None:
        dt_nyquist = np.pi / omega_max
        dt_est = min(dt_est, dt_nyquist)

    # Spatial resolution based on beam waist
    dx_est = w0 / 10
    dy_est = dx_est

    return dt_est, dx_est, dy_est

def refractive_index(omega,fiber:FIBER_config):
    omega = np.where(np.abs(omega) < 1e-12, 1e-12, omega)  # Avoid divide-by-zero
    lambda_m = 2 * np.pi * c / omega  # [m]
    lambda_um = lambda_m * 1e6  # [µm]
    lambda2 = lambda_um**2
    # Initialize output array
    n2 = np.full_like(lambda_um, np.nan)

    # Valid domain for fused silica Sellmeier model
    #valid = (lambda_um > 0.2) & (lambda_um < 3.5)
    valid = (lambda_um > 0.1) & (lambda_um < 1.5)
    lambda2 = lambda_um[valid] ** 2

    # Compute n^2 using Sellmeier formula
    n2[valid] = 1 + (fiber.B1 * lambda2) / (lambda2 - fiber.C1) + \
                     (fiber.B2 * lambda2) / (lambda2 - fiber.C2) + \
                     (fiber.B3 * lambda2) / (lambda2 - fiber.C3)

    # Take square root, avoid NaNs
    n = np.sqrt(n2)
    return np.where(np.isnan(n), 0, n)

# --- Raman response setup ---
def raman_response(simulation:SIMULATION_config):
    '''
    K. J. Blow, D. Wood, Theoretical description of transient
    stimulated Raman scattering in optical fibers.  IEEE J. Quantum Electron.,
    25 (1989) 1159, https://doi.org/10.1109/3.40655.
    '''
    tau1 = 12.2e-15  # s
    tau2 = 32e-15  # s
    f_R = 0.18       # Raman fractional contribution
    
    t = simulation.t
    hR = np.zeros_like(t)

    # Only positive times contribute (causal)
    t_pos_mask = t >= 0
    t_pos = t[t_pos_mask]

    # Compute hR only for t >= 0
    hR_pos = ((tau1**2 + tau2**2) / (tau1 * tau2**2)) * np.exp(-t_pos / tau2) * np.sin(t_pos / tau1)

    # Assign
    hR[t_pos_mask] = hR_pos

    # Normalize over positive times only
    norm = np.trapz(hR_pos, t_pos)
    if norm != 0:
        hR /= norm
    else:
        raise ValueError("Normalization integral is zero, check time vector resolution!")
    
    return f_R, hR

# Defining the Simulation function
def Simulation(fiber:FIBER_config,simulation:SIMULATION_config,laser: LASER_config):
    
    # Initial pulse A(x, y, t)
    X, Y, T = np.meshgrid(simulation.x, simulation.y, simulation.t, indexing='ij')
    A0 = GaussianPulse(T,laser.tau0, X, Y, laser.beam_waist, laser.amplitude)
    A0 = A0.astype(np.complex128)

    # --- Operators ---
    W = simulation.w
    omega = W + laser.omega0
    k_omega = (omega / c) * refractive_index(omega,fiber)
    k0 = laser.omega0 / c * refractive_index(laser.omega0, fiber)
    domega = W[1] - W[0]
    beta1 = np.gradient(k_omega, domega)[np.argmin(np.abs(W))]
    Dispersion = np.exp(1j * (k_omega - k0 - beta1 * (W))*simulation.dz)
    #Dispersion = np.exp(1j * ((1/2) * fiber.beta2 * W**2 - (1/6) * fiber.beta3 * W**3) * simulation.dz)
    Dispersion = Dispersion.reshape(1, simulation.nt)  # Shape: (1, nt) for correct broadcast over time dimension

    KX, KY = np.meshgrid(simulation.kx, simulation.ky, indexing='ij')
    KX2KY2 = KX**2 + KY**2
    Diffraction = np.exp(-1j * KX2KY2[:, :, None] / (2 * laser.k0) * simulation.dz)  # Shape (nx, ny, 1)

    Loss = np.exp(-(fiber.alpha_nepers_per_m / 2) * simulation.dz)

    f_R, hR = raman_response(simulation)

    # --- Storage ---
    A_snapshots = []
    A_snapshots.append(A0.copy())
    PhotonNumber_values = []
    PhotonNumber_values.append(getPhotonNumber(A0.copy(),simulation,fiber,laser))

    # --- Main loop ---
    for step_z in range(0,simulation.nz - 1):

        A = A_snapshots[step_z]
        I = getIntensity(A)

        # Nonlinearity half-step for Kerr
        Nonlinearity1 = np.exp(1j * fiber.gamma * I * simulation.dz / 2)
        A *= Nonlinearity1

        # Self-steepening half-step: d/dt of (I * A_out)
        NL = I * A
        
        # Step 1: FFT in time only (axis=2), keep (x, y) in real space
        NL_fft_time = fftshift(fft(ifftshift(NL, axes=2), axis=2), axes=2)

        # Step 2: Multiply by iω (temporal derivative)
        NL_fft_time *= 1j * W[None, None, :]  # W is fftshifted already

        # Step 3: IFFT in time (back to time domain, still spatially resolved)
        dNL_dt = fftshift(ifft(ifftshift(NL_fft_time, axes=2), axis=2), axes=2)

        # Self-steepening correction (half-step)
        Nonlinearity2 = 1j * (fiber.gamma / laser.omega0) * dNL_dt * simulation.dz / 2
        A += Nonlinearity2

        # --- Raman Term (half-step)---
        # Apply 1D convolution along time axis for all (x, y)
        raman_conv = convolve1d(I, hR, axis=2, mode='constant')
        raman_factor = (1 - f_R) * I + f_R * raman_conv * simulation.dt
        Nonlinearity3 = np.exp(1j * fiber.gamma * raman_factor * simulation.dz / 2)
        A *= Nonlinearity3

        # Step 1: Apply spatial FFT (x, y) → (kx, ky)
        A_fft_spatial = fftshift(fft2(ifftshift(A, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))

        # Step 2: Apply temporal FFT (t) → (w)
        A_fft_full = fftshift(fft(ifftshift(A_fft_spatial, axes=(2)), axis=2), axes=(2))

        # Step 3: Apply linear operator (dispersion)
        A_fft_full *= Dispersion * Loss * Diffraction  # Must be aligned with W (fftshifted)

        # Step 4: Inverse temporal FFT
        A_ifft_time = fftshift(ifft(ifftshift(A_fft_full, axes=(2)), axis=2), axes=(2))

        # Step 5: Inverse spatial FFT
        A_out = fftshift(ifft2(ifftshift(A_ifft_time, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))

        # Compute nonlinear terms (Kerr + self-steepening)
        
        I = getIntensity(A_out)

        # Nonlinearity half-step for Kerr
        Nonlinearity1 = np.exp(1j * fiber.gamma * I * simulation.dz / 2)
        A_out *= Nonlinearity1

        # Self-steepening half-step: d/dt of (I * A_out)
        NL = I * A_out
        
        # Step 1: FFT in time only (axis=2), keep (x, y) in real space
        NL_fft_time = fftshift(fft(ifftshift(NL, axes=2), axis=2), axes=2)

        # Step 2: Multiply by iω (temporal derivative)
        NL_fft_time *= 1j * W[None, None, :]  # W is fftshifted already

        # Step 3: IFFT in time (back to time domain, still spatially resolved)
        dNL_dt = fftshift(ifft(ifftshift(NL_fft_time, axes=2), axis=2), axes=2)

        # Self-steepening correction (half-step)
        Nonlinearity2 = 1j * (fiber.gamma / laser.omega0) * dNL_dt * simulation.dz / 2
        A_out += Nonlinearity2

        # --- Raman Term (half-step)---
        # Apply 1D convolution along time axis for all (x, y)
        raman_conv = convolve1d(I, hR, axis=2, mode='constant')
        raman_factor = (1 - f_R) * I + f_R * raman_conv * simulation.dt
        Nonlinearity3 = np.exp(1j * fiber.gamma * raman_factor * simulation.dz / 2)
        A *= Nonlinearity3

        A_snapshots.append(A_out.copy())

        PhotonNumber_values.append(getPhotonNumber(A0.copy(),simulation,fiber,laser))

        delta = int(round(step_z*100/nz)) - int(round((step_z-1)*100/nz))
        if delta == 1:
            print(str(int(round(step_z*100/nz))) + " % ready")
    # return results
    return A_snapshots, PhotonNumber_values

def savePlot(fileName):
    if not os.path.isdir('results/'):
        os.makedirs('results/')
    plt.savefig('results/%s.png'%(fileName))

def plotFirstAndLastPulse(Pulse,simulation:SIMULATION_config):
    # Initial vs Final Pulse along z-t, both at x = nx//2, y = ny//2
    # Initial and final pulse intensities at the center (x = nx//2, y = ny//2)
    initial_pulse = Pulse[0]
    initial_pulse = getIntensity(initial_pulse[simulation.nx // 2, simulation.ny // 2 , :])
    initial_pulse_maximum = np.max(initial_pulse)
    initial_pulse /= initial_pulse_maximum # Normalize
    final_pulse = Pulse[-1]
    final_pulse = getIntensity(final_pulse[simulation.nx // 2, simulation.ny // 2, :])
    final_pulse /= initial_pulse_maximum # Normalize
    # Create a figure for z-t evolution of initial and final pulses
    plt.figure(figsize=(10, 4))
    # Plot the intensity as a function of time for both the initial and final pulses
    plt.plot(simulation.t_fs, initial_pulse, label='Initial pulse')
    plt.plot(simulation.t_fs, final_pulse, label='Final pulse')
    plt.xlabel('Time [fs]')
    plt.ylabel('Intensity [a.u.]')
    plt.title('Initial vs final pulse at beam center')
    plt.legend()
    plt.tight_layout()
    savePlot('Initial vs final pulse at beam center')
    plt.show()

def plotFirstAndLastSpectrum(Pulse,simulation:SIMULATION_config):
    # Initial vs Final spectrum along z-t, both at x = nx//2, y = ny//2
    # Initial pulse and spectrum
    initial_pulse = Pulse[0][simulation.nx // 2, simulation.ny // 2, :]
    initial_spectrum = getIntensity(fftshift(fft(initial_pulse)))
    initial_spectrum_maximum = np.max(initial_spectrum)
    initial_spectrum /= initial_spectrum_maximum  # Normalize
    # Final pulse and spectrum
    final_pulse = Pulse[-1][simulation.nx // 2, simulation.ny // 2, :]
    final_spectrum = getIntensity(fftshift(fft(final_pulse)))
    final_spectrum /= initial_spectrum_maximum  # Normalize
    # Plot the comparison of initial and final spectrum
    plt.figure(figsize=(10, 4))
    plt.plot(simulation.f_PHz_rel*2*np.pi, initial_spectrum, label='Initial spectrum')
    plt.plot(simulation.f_PHz_rel*2*np.pi, final_spectrum, label='Final spectrum')
    plt.xlabel('Angular frequency [PHz]')
    plt.ylabel('Intensity [a.u.]')
    plt.title('Initial vs final spectrum at beam center')
    plt.legend()
    plt.tight_layout()
    savePlot('Initial vs final spectrum at beam center')
    plt.show()

def plotPeakIntensity(Pulse,simulation:SIMULATION_config):
    # --- Plot the peak intensity evolution ---
    Peak_values = [np.max(getIntensity(A[simulation.nx // 2, simulation.ny // 2, :])) for A in Pulse]
    plt.figure(figsize=(10, 4))
    plt.plot(simulation.z, Peak_values)
    plt.xlabel('Propagation distance [m]')
    plt.ylabel('Peak Intensity [W/m^2]')
    plt.title('Peak intensity evolution at center')
    savePlot('Peak intensity evolution at center')
    plt.tight_layout()
    plt.show()

def plotXYBeamprofile(Pulse,simulation:SIMULATION_config):
    # --- Beam profile in x-y plane at final z and final time t ---
    final_field = Pulse[-1]  # shape: (nx, ny, nt)
    # Extract intensity slice at final time
    intensity_xy = getIntensity(final_field[:, :, -1])
    intensity_xy /= np.max(intensity_xy)  # normalize
    X, Y = np.meshgrid(simulation.x,simulation.y)
    plt.figure(figsize=(10, 4))
    plt.contourf(X, Y, intensity_xy, levels=100, cmap='inferno')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Beam profile in x–y plane (at final z, final t)')
    plt.colorbar(label='Intensity [a.u.]')
    savePlot('XY Beamprofile')
    plt.tight_layout()
    plt.show()

def plotPulseEvolution(Pulse,simulation:SIMULATION_config):
    # --- Pulse Evolution at Beam Center (z–t View) ---
    # Build 2D array: rows = z steps, columns = time
    zt_matrix = np.array([getIntensity(A[simulation.nx // 2, simulation.ny // 2, :]) for A in Pulse])  # shape: (nz, nt)
    zt_matrix /= np.max(zt_matrix)  # normalize
    # Plot (z, t) image
    extent_zt = [simulation.t_fs[0], simulation.t_fs[-1], 0, simulation.z_max]
    plt.figure(figsize=(10, 4))
    plt.imshow(zt_matrix, extent=extent_zt, aspect='auto', origin='lower', cmap='inferno')
    plt.xlabel('Time [fs]')
    plt.ylabel('Propagation distance z [m]')
    plt.title('Pulse evolution at beam center (z–t view)')
    plt.colorbar(label='Intensity [a.u.]')
    savePlot('Pulse evolution')
    plt.tight_layout()
    plt.show()

def plotSpectrumEvolution(Pulse,simulation:SIMULATION_config):
    # --- Spectrum Evolution at Beam Center (z–t View) ---
    # FFT along time at each z, at the center (x, y)
    spectrum_zw = [
        getIntensity(fftshift(fft(A[simulation.nx // 2, simulation.ny // 2, :])))
        for A in Pulse
    ]
    spectrum_zw = np.array(spectrum_zw)
    # Normalize (optional)
    spectrum_zw /= np.max(spectrum_zw)
    # Plot
    W, Z = np.meshgrid(simulation.f_PHz_rel*2*np.pi, simulation.z)
    plt.figure(figsize=(10, 4))
    plt.contourf(W, Z, spectrum_zw, levels=40, cmap='inferno')
    plt.xlabel('Angular frequency [PHz]')
    plt.ylabel('Propagation distance z [m]')
    plt.title('Spectrum evolution at beam center (z–frequency view)')
    plt.colorbar(label='Intensity [a.u.]')
    savePlot('Spectrum evolution')
    plt.tight_layout()
    plt.show()

def plotWavelength(Pulse, simulation: SIMULATION_config, laser: LASER_config):
    # Initial vs Final spectrum as function of wavelength along z-t, both at x = nx//2, y = ny//2
    # Compute temporal FFT of the initial pulse at (x=nx//2, y=ny//2)
    initial_pulse = Pulse[0][simulation.nx // 2, simulation.ny // 2, :]
    final_pulse = Pulse[-1][simulation.nx // 2, simulation.ny // 2, :]
    initial_spectrum = getIntensity(fftshift(fft(initial_pulse)))  # shape: (nt,)
    final_spectrum = getIntensity(fftshift(fft(final_pulse)))  # shape: (nt,)
    # Define angular frequency and corresponding wavelength
    omega = simulation.w + laser.omega0  # [rad/s], shape: (nt,)
    omega = np.where(np.abs(omega) < 1e-12, 1e-12, omega)
    wavelength = 2 * np.pi * c / omega  # [m], shape: (nt,)
    # Convert spectral intensity I(ω) to I(λ) using: I(λ) ∝ I(ω) * dω/dλ ∝ I(ω) * (2πc / λ²)
    initial_spectrum_lambda = initial_spectrum * 2 * np.pi * c / wavelength**2
    final_spectrum_lambda = final_spectrum * 2 * np.pi * c / wavelength**2
    # Apply wavelength window
    wavelength0 = simulation.lambda0
    valid = (wavelength > 0.5 * wavelength0) & (wavelength < 1.5 * wavelength0)
    # Select valid portion
    wavelength_valid = wavelength[valid]
    initial_spectrum_valid = initial_spectrum_lambda[valid]
    final_spectrum_valid = final_spectrum_lambda[valid]
    # Plot
    plt.figure(figsize=(10, 4))
    plt.plot(wavelength_valid * 1e9, initial_spectrum_valid / np.max(initial_spectrum_valid), label="Initial Spectrum")  # Optional: convert to nm
    plt.plot(wavelength_valid * 1e9, final_spectrum_valid / np.max(initial_spectrum_valid), label="Final Spectrum")  # Optional: convert to nm
    plt.title("Intensity vs. Wavelength")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Normalized Intensity")
    plt.legend()
    savePlot("Intensity as function of wavelength")
    plt.show()

def plotPhotonNumberValues(PhotonNumber_values,simulation: SIMULATION_config):
    # --- Plot the peak intensity evolution ---
    plt.figure(figsize=(10, 4))
    plt.plot(simulation.z, PhotonNumber_values)
    plt.xlabel('Propagation distance [m]')
    plt.ylabel('Photon number [count]')
    plt.title('Photon number conservation')
    savePlot('Photon number conservation')
    plt.tight_layout()
    plt.show()