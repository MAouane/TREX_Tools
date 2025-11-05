import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap

''' Constants and functions to convert between ToF, Energy in meV and Wavelengths in Å'''
TOF_CONSTANT = 252.78  # µs·Å-1·m-1 for conversion between lambda and ToF
ENERGY_CONSTANT = 81.8042  # meV·Å**2
H_OVER_MN = 3956  # m/s

''' T-REX Specific Distances '''
L_sample = 163.8   # Source to sample in m
L_detector = 3.0   # Sample to detector in m
L_total = L_sample + L_detector
LM = 162.0  # Distance to the M-chopper

# Conversion functions
def energy_to_wavelength(E):
    return np.sqrt(ENERGY_CONSTANT / E)

def wavelength_to_energy(lam):
    return ENERGY_CONSTANT / lam**2

def Lechner(Ld, lambda_f, lambda_i, tau_M, tau_P, Ld_mm_conversion=1e-3):
    '''
    Calculate energy resolution using Lechner Formula.
    Inputs:
        - Ld: Any flight path uncertainties, in mm
        - lambda_f: Scattered wavelength in Å
        - lambda_i: Incident wavelength in Å
        - tau_M: M-chopper opening time (seconds)
        - tau_P: P-chopper opening time (seconds)
    Outputs:
        - Returns delta E / FWHM in meV
    '''
    # Convert Ld from mm to m
    Ld = Ld * Ld_mm_conversion

    A = (0.2041*tau_M) * (L_PM + L_MS + L_SD * (lambda_f**3 / lambda_i**3))  # [A] = s * m
    B = (0.2887*tau_P) * (L_MS + L_SD * (lambda_f**3 / lambda_i**3))         # [B] = s * m
    C = mn * L_PM * lambda_f * Ld / h                                                                 # [C] = s * m

    # Final energy resolution calculation
    return ((h**3.)/(mn**2)) * np.sqrt((A*A) + (B*B) + (C*C)) * 1./(L_SD * L_PM * (lambda_f**3.)) * 6.21 * 1e21  # Convert Joules -> meV

import scipy.constants as cte

# Constants
h = cte.h                          # Planck constant
mn = cte.m_n                       # Neutron mass
L_SD = 3.0                         # Sample to detector distance in m
L_PM = LM - 108.0                  # P- to M-chopper distance in m
#print(L_PM)
L_MS = 1.8                         # M-chopper to sample distance in m

# Figure parameters
Ld = 20.0  # Flight path uncertainty in mm
steps = np.arange(1, 0., -0.05)  # Range of energy loss
Lambda_Reps = [1.65]#[1, 1.2, 1.4, 1.65]
#Lambda_Reps = [2.6, 3.5, 5, 6.5]
n_lines = len(Lambda_Reps)

fig, ax = plt.subplots(figsize=(16, 8))
#cmap = plt.get_cmap('coolwarm_r')

# Define setups for the two conditions
setups = [
    {"fM": 252, "M_opening":4.4, "P_opening":35, "P_num":2, "tau_factor_m": 2, "tau_factor_p": 2, "label_desc": "2 M-Chopper disks - High Flux", "cmap": "inferno", 'linestyle': '--', 'marker': 'x'},
    {"fM": 252, "M_opening":2.5, "P_opening":20, "P_num":2, "tau_factor_m": 2, "tau_factor_p": 2, "label_desc": "2 M-Chopper disks - High Resolution", "cmap": "jet", 'linestyle': '--', 'marker': 'o'}
]

print('-------------------------------------------------------------------------')
print('Inelastic resolution calculation for a given wavelength')
print('-------------------------------------------------------------------------')
print('Script uses the Lechner formula for the resolution calculation')
print('-------------------------------------------------------------------------')
print('Script will warn you if your pulse/width ratio is not optimised for flux/resolution')
print('-------------------------------------------------------------------------')
print('High resolution openings of the PS- and M-choppers used in this script')
print('-------------------------------------------------------------------------')
print('Triangular transmission function at the M-chopper: τM is multiplied by 0.2041')
print('-------------------------------------------------------------------------')
print('Top-Hat transmission function at the PS-chopper: τP is multiplied by 0.2887')
print('-------------------------------------------------------------------------')


for setup in setups:
    fM = setup["fM"]
    P_opening = setup["P_opening"]
    M_opening = setup["M_opening"]
    tau_factor_m = setup["tau_factor_m"]
    tau_factor_p = setup["tau_factor_p"]
    linestyle = setup["linestyle"]
    label_desc = setup["label_desc"]
    marker = setup["marker"]
    cmap_ = setup["cmap"]
    fP = 0.75 * fM  # P-Chopper frequency
    cmap = plt.get_cmap(cmap_)

    # Calculate chopper opening times
    tau_M = M_opening / tau_factor_m / 360. / fM  # M-chopper opening time in s
    if tau_factor_m == 2:
        tau_P = tau_M * (L_PM + L_MS + L_SD)/(L_MS + L_SD) #P_opening / tau_factor_p / 360. / fP         # P-chopper opening time in s
    else:
        tau_P = 0.5 * tau_M * (L_PM + L_MS + L_SD)/(L_MS + L_SD) #P_opening / tau_factor_p / 360. / fP         # P-chopper opening time in s
        
    for i in range(len(Lambda_Reps)):
        dE_vals_Matched = []
        Res_vals_Matched = []
        color = cmap(i / n_lines)

        # Generate label only once for each line (based on wavelength)
        #label = f'λ = {round(Lambda_Reps[i], 2)}Å - {label_desc}, Optimised'
        label = f'T-REX - fM = {fM}Hz, {label_desc}'

        for j in range(len(steps)):
            lambda_f = np.sqrt(1./steps[j]) * Lambda_Reps[i]
            lambda_i = Lambda_Reps[i]
            Resolution_FWHM_Matched = np.sqrt(8*np.log(2))*Lechner(Ld, lambda_f * 1e-10, lambda_i * 1e-10, tau_M, tau_P)

            # Append dE_vals and Res_vals
            dE_vals_Matched.append(wavelength_to_energy(lambda_i) - wavelength_to_energy(lambda_f))
            Res_vals_Matched.append(Resolution_FWHM_Matched)

        # Plot line with label
        ax.plot(dE_vals_Matched, Res_vals_Matched, linestyle, color=color, label=label, marker = marker, markersize = 5)


print('Ratio between opening times: ', tau_P/tau_M)
print('Distance Ratio on T-REX: ', (L_SD+L_MS+L_PM)/(L_SD+L_MS))
if tau_P/tau_M == (L_SD+L_MS+L_PM)/(L_SD+L_MS):
    print('Optimised Performance for Flux/Resolution Trade-Off')
else:
    print('Unoptimised Instrument. Check the PS- and M-choppers Opening Times')
print('-------------------------------------------------------------------------')

# Final adjustments
ax.set_xlabel("Energy Transfer (meV)", fontsize=14)
ax.set_ylabel("FWHM (meV)", fontsize=14)
ax.legend(fontsize = 15)  # Automatically use labels from ax.plot()
ax.set_title(f"T-REX Performance Comparison - λi = {Lambda_Reps[0]}Å - Optimised" , fontsize = 15)
#ax.set_title(f"T-REX Performance - Thermal Regime - fM = {fM}Hz", fontsize = 15)
#ax.set_title(f"T-REX Performance - Cold Regime - fM = {fM}Hz", fontsize = 15)
plt.grid(True)
plt.minorticks_on()
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.tight_layout()
plt.show()