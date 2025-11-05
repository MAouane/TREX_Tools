import numpy as np
import matplotlib.pyplot as plt

''' Constants and functions to convert between ToF, Energy in meV and Wavelengths in Å'''
TOF_CONSTANT = 252.78  # µs·Å-1·m-1 for conversion between lambda and ToF
ENERGY_CONSTANT = 81.8042  # meV·Å**2
H_OVER_MN = 3956  # m/s

# Conversion functions
def energy_to_wavelength(E):
    return np.sqrt(ENERGY_CONSTANT / E)

def wavelength_to_energy(lam):
    return ENERGY_CONSTANT / lam**2


def calculate_detector_coverage(E_i, low_angle, high_angle, num_points=1000):
    """
    Calculate and plot the detector coverage in Q, Delta E space for a given Ei, low angle, and high angle.

    Parameters:
    E_i: float
        Incident energy in meV.
    low_angle: float
        Low detector angle in degrees.
    high_angle: float
        High detector angle in degrees.
    num_points: int
        Number of points used for smooth curves in the plot.
    """
    # Constant for the expression k = sqrt(E / C)
    C = 2.072  # meV·Å^2

    # Calculate incident wavevector, ki (Å^-1)
    k_i = np.sqrt(E_i / C)

    # Generate energy transfer range from -Ei to +Ei
    if E_i > 12:
        delta_E = np.linspace(-0.30*E_i, E_i, num_points)
    else:
        delta_E = np.linspace(-0.50*E_i, E_i, num_points)

    # Initialize arrays for Q values at low and high angles
    Q_low = np.zeros_like(delta_E)
    Q_high = np.zeros_like(delta_E)

    # Calculate Q vectors for both low and high angles
    for i, dE in enumerate(delta_E):
        E_f = E_i - dE  # Final energy in meV
        if E_f >= 0:  # Only compute for valid final energies
            k_f = np.sqrt(E_f / C)  # Final wavevector (Å^-1)

            theta_low = np.radians(low_angle)  # Convert low angle to radians
            theta_high = np.radians(high_angle)  # Convert high angle to radians

            # Q calculation for low and high detector angles
            Q_low[i] = np.sqrt(k_i**2 + k_f**2 - 2 * k_i * k_f * np.cos(theta_low))
            Q_high[i] = np.sqrt(k_i**2 + k_f**2 - 2 * k_i * k_f * np.cos(theta_high))
        else:
            Q_low[i] = np.nan  # Assign NaN for invalid values
            Q_high[i] = np.nan

    # Plotting the detector coverage
    if high_angle == 144:
        plt.plot(Q_low, delta_E, label=f'Full Coverage - λi = {np.round(energy_to_wavelength(E_i),2)}Å', color='red', linestyle='-')
        plt.plot(Q_high, delta_E, color='red', linestyle='-') #label=f'High Angle: {high_angle}°',
    else:
        plt.plot(Q_low, delta_E, label=f'Day 1 Coverage - λi = {np.round(energy_to_wavelength(E_i),2)}Å', color='blue', linestyle='--')
        plt.plot(Q_high, delta_E, color='blue', linestyle='--') #label=f'High Angle: {high_angle}°',

    # Add labels, title, and grid to the plot
    if E_i >= 12:
        plt.title('Thermal Regime - Detector Coverage in Q, $\Delta E$ Space', fontsize=15)
    else:
        plt.title('Cold Regime - Detector Coverage in Q, $\Delta E$ Space', fontsize=15)
        
    plt.xlabel(r'Momentum Transfer Q ($\AA^{-1}$)', fontsize=15)
    plt.ylabel(r'$\Delta E$ (meV)', fontsize=15)
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)  # Line at Delta E = 0
    plt.legend(fontsize = 15)




plt.figure(figsize=(16, 10))

Ei = [120, 81.8, 30, 12]#, 3.86, 1.28]
for i in range(len(Ei)):
    calculate_detector_coverage(Ei[i], 1, 144, 2000)
    calculate_detector_coverage(Ei[i], 1, 72, 2000)
plt.grid()
plt.minorticks_on()
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.tight_layout()
plt.show()

plt.figure(figsize=(16, 10))

lambda_i = [3.5, 4.6, 6.5, 8]
for i in range(len(lambda_i)):
    Ei = wavelength_to_energy(lambda_i[i])
    calculate_detector_coverage(Ei, 1, 144, 2000)
    calculate_detector_coverage(Ei, 1, 72, 2000)


plt.grid()
plt.minorticks_on()
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.tight_layout()
plt.show()