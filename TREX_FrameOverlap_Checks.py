# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

''' Constants and functions to convert between ToF, Energy in meV and Wavelengths in Å'''
TOF_CONSTANT = 252.78  # µs·Å-1·m-1 for conversion between lambda and ToF
ENERGY_CONSTANT = 81.8042  # meV·Å**2
H_OVER_MN = 3956  # m/s
BANDWIDTH = 0.845      # Half-bandwidth for wavelength range

''' T-REX Specific Distances to M-chopper, Sample and Detector for accurate calculations '''
L_sample = 163.8   #source to sample in m
L_detector = 3.0   #sample to detector in m
L_total = L_sample + L_detector
L_fan = 151.0 #source to FAN chopper in m
LM = 162.0  #Distance to the M-chopper
T_offset = 1700

#Conversion functions
def energy_to_wavelength(E):
    return np.sqrt(ENERGY_CONSTANT / E)

def wavelength_to_energy(lam):
    return ENERGY_CONSTANT / lam**2

def tof_to_lambda(tof):
    return tof / (TOF_CONSTANT * L_sample)

def lambda_to_tof_fan(lam):
    return T_offset + TOF_CONSTANT * lam * L_sample

def lambda_to_tof_sample(lam):
    return T_offset + TOF_CONSTANT * lam * L_sample
    
def lambda_to_tof_detector(lam):
    return T_offset + TOF_CONSTANT * lam * L_total

def LambdaLoss(energy_loss, lam):
    scale_factor = np.sqrt(1 / energy_loss)
    return np.array(lam) * scale_factor

def generate_wavelengths(central_wavelength, fM):    
    """    Generate evenly spaced wavelengths based on fM, ensuring λ_c is included.    
    Parameters:        
    central_wavelength (float): Central wavelength (Å).        
    fM (float): Chopper frequency (Hz; between 14 and 336).    
    Returns:        
    list: List of wavelengths in Å.    """    
    lambda_min = central_wavelength - BANDWIDTH    
    lambda_max = central_wavelength + BANDWIDTH    
    delta_lambda = 3956. / (fM * LM)    
    wavelengths_below = []    
    current_lambda = central_wavelength    
    while current_lambda > lambda_min:        
        current_lambda -= delta_lambda        
        if current_lambda >= lambda_min:            
            wavelengths_below.append(current_lambda)    
    wavelengths_above = []    
    current_lambda = central_wavelength    
    while current_lambda < lambda_max:        
        current_lambda += delta_lambda        
        if current_lambda <= lambda_max:            
            wavelengths_above.append(current_lambda)    
    wavelengths_below.reverse()    
    return wavelengths_below + [central_wavelength] + wavelengths_above


print('-------------------------------------------------------------------------')
print('Script to determine the possible Eis with RRM and possible frame overlap for T-REX \n- Change the M-chopper frequency at line 64\n- Change your central wavelength in line 51')

print('-------------------------------------------------------------------------')
print('The ToFs at sample position and at the detectors have been compared with McStas outputs for a few settings!')



central_wavelength = 3.5 #Central wavelength for your setting
if central_wavelength > 7.0 or central_wavelength < 0.5:
    print('-------------------------------------------------------------------------')
    print('Well... Maybe not?') #Users can run something above 7Å if they want, it's just CSPEC and MIRACLES are optimised for longer wavelengths. No 'hard' limits
else:
    print('-------------------------------------------------------------------------')
    print(f'Your central wavelength is: {central_wavelength} Å')
    
# Rounding all of the floating points to .2f format to avoid any and all floating errors... Not setting the .2f messes up the subpulses for some reason
delta_lambda = round((3956)*(1/(14*166.8)), 2)    # Total wavelength bandwidth of T-REX, should be in the 1.7Å range
print(f'Natural bandwidth of the instrument is: {delta_lambda}Å')
print(f'Wavelength range is delimited by: +/-{delta_lambda/2.}Å')

fM = 14*12 # M-Chopper frequency in Hz (must be multiple of 14 between 14–336). The multiple defines how many reps you get on the sample
fP = fM*0.75 # M-Chopper frequency in Hz. Added as part of the contingency


frame_period_us = 71000  # µs (71 ms) - One ESS pulse duration. To check if longer wavelengths cross into the the fast neutron of the following frame

if fM % 14 != 0 or not (14 <= fM <= 336):
    raise ValueError("fM (M-chopper frequency) must be a multiple of 14 between 14 and 336Hz")
'''
if fP > 140:
    raise ValueError(f"fP (PS-chopper frequency) must be below 140Hz for the Al disks. Your frequency of {fP}Hz exceeds it")
''' 

print(f'For your chosen frequency of {fM}Hz, you get {int(fM/14)} Ei Reps')
print('-------------------------------------------------------------------------')
print(f'For your chosen frequency of {fM}Hz, Tthe PS chopper is running at: {fP}Hz')
print('-------------------------------------------------------------------------')

d_lambda = (H_OVER_MN / (fM * LM)) #Step size for subpulses
print(f'Step size for this setting is: {d_lambda}Å')
lambda_min = central_wavelength - delta_lambda / 2.
lambda_max = central_wavelength + delta_lambda / 2.

print(f'Your limits are {lambda_min}Å and {lambda_max}Å')
print('-------------------------------------------------------------------------')

reps = int(fM/14)
Lambda_Reps = generate_wavelengths(central_wavelength,fM)
Lambda_Reps = np.array(Lambda_Reps)

print('Your unfiltered incident wavelengths (in Å) from RRM are: \n', Lambda_Reps)
print('-------------------------------------------------------------------------')

print('Your unfiltered incident energies (in meV) from RRM are: \n', wavelength_to_energy(Lambda_Reps))
print('-------------------------------------------------------------------------')

for i in range(len(Lambda_Reps)):
    if Lambda_Reps[i] < 0.67:
        print('-------------------------------------------------------------------------')
        print(f'{Lambda_Reps[i]}Å is outside the optimal range of T-REX and ESS')
        print('-------------------------------------------------------------------------')

ToF_FanPos = lambda_to_tof_fan(Lambda_Reps)
ToF_SamplePos = lambda_to_tof_sample(Lambda_Reps)
ToF_DetectorPos_Lambda_i = lambda_to_tof_detector(Lambda_Reps)

#print('RRM ToFs at sample position in seconds are: ', ToF_SamplePos*1e-6)
#print('RRM ToFs at detector position in seconds are: ', ToF_DetectorPos_Lambda_i*1e-6)


print('The limit for energy loss is set to 80% of each thermal Ei\nFor cold energies, it is set to 20%\nGain is set to 80% as a test')

Lambda_RRM_Loss = np.zeros_like(Lambda_Reps)
E_RRM_Loss = np.zeros_like(Lambda_Reps)
ToF_RRM_Loss = np.zeros_like(Lambda_Reps)
ToF_RRM_Gain = np.zeros_like(Lambda_Reps)
ToF_RRM_Loss_Detector = np.zeros_like(Lambda_Reps)
ToF_RRM_Gain_Detector = np.zeros_like(Lambda_Reps)
Gain_fraction = 0.5  # 50% gain in neutron energy for both thermal and cold wavelengths

''' Originally, I set all cold wavelength to only have 20% of Loss since you'd only look at low energy transfers/QENS but changed it to be 85% like the thermal wavelengths '''

for i in range(len(Lambda_Reps)):    
    ToF_RRM_Gain_Detector[i] = TOF_CONSTANT * Lambda_Reps[i] * (np.sqrt(1. / (1. + Gain_fraction))) * L_detector 
    ToF_RRM_Gain[i] = ToF_RRM_Gain_Detector[i] + ToF_SamplePos[i]
    
    if Lambda_Reps[i] < 2.6:
        ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.15) * Lambda_Reps[i] * L_detector 
    else:
        ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.15) * Lambda_Reps[i] * L_detector 
    
    ToF_RRM_Loss[i] = ToF_RRM_Loss_Detector[i] + ToF_SamplePos[i]

Lambda_RRM_Loss = tof_to_lambda(ToF_RRM_Loss)
E_RRM_Loss = wavelength_to_energy(Lambda_RRM_Loss)

print('To compare the ToFs with McStas, the same offset was added to all ToFs')
print('-------------------------------------------------------------------------')
print('ToFs at the FAN chopper position for this setting in s are: \n', ToF_FanPos*1e-6)
print('-------------------------------------------------------------------------')
print('ToFs at the sample position for this setting in s are: \n', ToF_SamplePos*1e-6)
print('-------------------------------------------------------------------------')
print('Elastic ToF at the detectors in s: \n', ToF_DetectorPos_Lambda_i*1e-6)
print('-------------------------------------------------------------------------')
print('Inelastic (Loss) ToF at the detectors in s: \n', ToF_RRM_Loss*1e-6)
print('-------------------------------------------------------------------------')
print('Inelastic (Gain) ToF at the detectors in s: \n', ToF_RRM_Gain*1e-6)
print('-------------------------------------------------------------------------')


y_vals = [0, L_detector]
fig, ax = plt.subplots(figsize=(12, 6))

cmap = plt.get_cmap('jet_r')
n_lines = len(Lambda_Reps)

for i in range(n_lines):
    color = cmap(i / n_lines)
    
    #Elastic line
    ax.plot([ToF_SamplePos[i]*1e-3, ToF_DetectorPos_Lambda_i[i]*1e-3], y_vals, color=color, linewidth=2, label = f'Ei: {wavelength_to_energy(Lambda_Reps[i]):.2f}meV')
    ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_DetectorPos_Lambda_i[i]*1e-3+71], y_vals, color=color, linewidth=2)
    
    #Inelastic (loss) line
    ax.plot([ToF_SamplePos[i]*1e-3, ToF_RRM_Loss[i]*1e-3], y_vals, color=color, linestyle='dashed', linewidth=1.2)
    ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_RRM_Loss[i]*1e-3+71], y_vals, color=color, linestyle='dashed', linewidth=1.2)
    
    #Inelastic (loss) line
    ax.plot([ToF_SamplePos[i]*1e-3, ToF_RRM_Gain[i]*1e-3], y_vals, color=color, linestyle='dashed', linewidth=1.2)
    ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_RRM_Gain[i]*1e-3+71], y_vals, color=color, linestyle='dashed', linewidth=1.2)

ax.vlines(71, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*2, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*3, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*4, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*5, 0, 3, 'k', linestyle = 'solid')

ax.set_xlabel("Time-of-Flight (ms)")
ax.set_ylabel("Distance (m)")
ax.legend()
ax.set_title(f"Unfiltered Time–Distance Diagram (Elastic and Inelastic Neutron Paths) - {central_wavelength}Å")
plt.xlim(ToF_RRM_Gain[0]*1e-3-5, ToF_RRM_Loss[-1]*1e-3+75)
plt.grid(True)
plt.tight_layout()
plt.show()



''' Up to here, the script showed you all of the possible reps for a given M-chopper frequency.
    The following part, looks at your plots and tells you which ones need supressing with the FAN chopper on T-REX '''

print('-------------------------------------------------------------------------')
print('Checking for frame overlap between subpulses and subsequent ESS frames')
print('-------------------------------------------------------------------------')

subpulse_overlapping_indices = set() 

for i in range(len(Lambda_Reps) - 1): 
    #Elastic ToF of next subpulse
    elastic_next = ToF_DetectorPos_Lambda_i[i+1]
    
    #Current subpulse ToFs
    elastic_current = ToF_DetectorPos_Lambda_i[i]
    loss_current = ToF_RRM_Loss[i]
    
    #ToF (gain) of next subpulse
    gain_next = ToF_RRM_Gain[i+1]
    
    # First Check: Elastic of λ[i+1] within Loss of λ[i]
    if elastic_current <= elastic_next <= loss_current:
        subpulse_overlapping_indices.add(i+1)  # Suppress λ[i+1]
        print(f"Overlap detected: λi = {Lambda_Reps[i+1]:.2f} Å suppressed due overlap with λi = {Lambda_Reps[i]:.2f} Å")
        print(f"Overlap detected: Ei= {wavelength_to_energy(Lambda_Reps[i+1]):.2f} Å suppressed due to overlap with Ei = {wavelength_to_energy(Lambda_Reps[i]):.2f} Å")
    
    # Second Check: Loss of λ[i] overlaps with Gain of λ[i+1]
    if gain_next <= loss_current:
        subpulse_overlapping_indices.add(i+1)  # Suppress λ[i+1]
        print(f"Overlap detected: λi= {Lambda_Reps[i+1]:.2f} Å suppressed due to gain-loss overlap with λ = {Lambda_Reps[i]:.2f} Å")
        print(f"Overlap detected: Ei= {wavelength_to_energy(Lambda_Reps[i+1]):.2f} Å suppressed due to overlap with Ei = {wavelength_to_energy(Lambda_Reps[i]):.2f} Å")

# Filtering the unwanted lambdas out!

valid_subpulse_indices = [i for i in range(len(Lambda_Reps)) if i not in subpulse_overlapping_indices]

Lambda_Reps_filtered_subpulses = Lambda_Reps[valid_subpulse_indices]
ToF_SamplePos_filtered_subpulses = ToF_SamplePos[valid_subpulse_indices]
ToF_DetectorPos_Lambda_i_filtered_subpulses = ToF_DetectorPos_Lambda_i[valid_subpulse_indices]
ToF_RRM_Loss_filtered_subpulses = ToF_RRM_Loss[valid_subpulse_indices]
ToF_RRM_Gain_filtered_subpulses = ToF_RRM_Gain[valid_subpulse_indices]



overlapping_indices = set()
frame_period_us = 71000 

#Within two consecutive frames of the chosen setting, check which two wavelengths overlap.
for i in range(len(Lambda_Reps_filtered_subpulses)):
    E_i = wavelength_to_energy(Lambda_Reps_filtered_subpulses[i])
    elastic_i = ToF_DetectorPos_Lambda_i_filtered_subpulses[i]
    inelastic_i = ToF_RRM_Loss_filtered_subpulses[i]
    for j in range(len(Lambda_Reps_filtered_subpulses)):
        if i == j:
            continue
        E_j = wavelength_to_energy(Lambda_Reps_filtered_subpulses[j])
        elastic_j = ToF_DetectorPos_Lambda_i_filtered_subpulses[j] + frame_period_us
        inelastic_j = ToF_RRM_Loss_filtered_subpulses[j] + frame_period_us
        if (elastic_i <= inelastic_j and inelastic_i >= elastic_j):
            if E_j > 140:
                print(f"Cross-frame overlap: λ={Lambda_Reps_filtered_subpulses[i]:.2f}Å (frame 0) " f"and λ={Lambda_Reps_filtered_subpulses[j]:.2f}Å (frame +1), but energy of λ_j = {E_j:.1f} meV > 140: Keeping both.")
            else:
                idx_to_suppress = i if Lambda_Reps_filtered_subpulses[i] > Lambda_Reps_filtered_subpulses[j] else j
                overlapping_indices.add(idx_to_suppress)
                print(f"Frame overlap detected: λi={Lambda_Reps_filtered_subpulses[i]:.2f}Å and λi={Lambda_Reps_filtered_subpulses[j]:.2f}Å: Suppressing the longer wavelength.")
                print(f"Frame overlap detected: Ei={wavelength_to_energy(Lambda_Reps_filtered_subpulses[i]):.2f}meV and Ei={wavelength_to_energy(Lambda_Reps_filtered_subpulses[j]):.2f}meV: Suppressing the lower energy.")



#New list of lambdas that only keeps the indices that have no overlap after checking
valid_indices = [i for i in range(len(Lambda_Reps_filtered_subpulses)) if i not in overlapping_indices]

# Apply filter to all the ToF arrays before plotting
Lambda_Reps_filtered = Lambda_Reps_filtered_subpulses[valid_indices]
ToF_SamplePos_filtered = ToF_SamplePos_filtered_subpulses[valid_indices]
ToF_DetectorPos_Lambda_i_filtered = ToF_DetectorPos_Lambda_i_filtered_subpulses[valid_indices]
ToF_RRM_Loss_filtered = ToF_RRM_Loss_filtered_subpulses[valid_indices]
ToF_RRM_Gain_filtered = ToF_RRM_Gain_filtered_subpulses[valid_indices]


y_vals = [0, L_detector]
fig, ax = plt.subplots(figsize=(12, 6))

cmap = plt.get_cmap('jet_r')
n_lines = len(Lambda_Reps_filtered)

for i in range(n_lines):
    color = cmap(i / n_lines)
    
    #Elastic line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_DetectorPos_Lambda_i_filtered[i]*1e-3], y_vals, color=color, linewidth=2, label = f'Ei: {wavelength_to_energy(Lambda_Reps_filtered[i]):.2f}meV')
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_DetectorPos_Lambda_i_filtered[i]*1e-3+71], y_vals, color=color, linewidth=2)
    
    #Inelastic (loss) line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_RRM_Loss_filtered[i]*1e-3], y_vals, color=color, linestyle='dashed', linewidth=1.2)
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_RRM_Loss_filtered[i]*1e-3+71], y_vals, color=color, linestyle='dashed', linewidth=1.2)

    #Inelastic (gain) line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_RRM_Gain_filtered[i]*1e-3], y_vals, color=color, linestyle='dashed', linewidth=1.2)
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_RRM_Gain_filtered[i]*1e-3+71], y_vals, color=color, linestyle='dashed', linewidth=1.2)

ax.vlines(71, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*2, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*3, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*4, 0, 3, 'k', linestyle = 'solid')
ax.vlines(71*5, 0, 3, 'k', linestyle = 'solid')

ax.set_xlabel("Time-of-Flight (ms)")
ax.set_ylabel("Distance (m)")
ax.legend()
ax.set_title(f"Filtered Time–Distance Diagram (Elastic and Inelastic Neutron Paths) - {central_wavelength}Å")
plt.xlim(ToF_RRM_Gain_filtered[0]*1e-3-5, ToF_RRM_Loss_filtered[-1]*1e-3+75)
plt.grid(True)
plt.tight_layout()
plt.show()

# Updated section to suppress either odd or even subpulses
print('-------------------------------------------------------------------------')
print('Suppression of either odd or even subpulses based on user preference')
print('-------------------------------------------------------------------------')

# User selection for suppression
suppress_odd = True  # Set this to True to suppress odd subpulses, False to suppress even ones

if suppress_odd:
    print("Suppressing odd subpulses...")
    valid_indices = [i for i in range(len(Lambda_Reps)) if i % 2 == 0]  # Keep indices divisible by 2 (even indices)
else:
    print("Suppressing even subpulses...")
    valid_indices = [i for i in range(len(Lambda_Reps)) if i % 2 != 0]  # Keep odd indices

# Filter subpulses according to the selected indices
Lambda_Reps_filtered = Lambda_Reps[valid_indices]
ToF_SamplePos_filtered = ToF_SamplePos[valid_indices]
ToF_DetectorPos_Lambda_i_filtered = ToF_DetectorPos_Lambda_i[valid_indices]
ToF_RRM_Loss_filtered = ToF_RRM_Loss[valid_indices]
ToF_RRM_Gain_filtered = ToF_RRM_Gain[valid_indices]

# Print the results for confirmation
print('Filtered incident wavelengths (in Å): \n', Lambda_Reps_filtered)
print('Filtered incident energies (in meV): \n', wavelength_to_energy(Lambda_Reps_filtered))

# Plotting the Time–Distance Diagram after filtering
y_vals = [0, L_detector]
fig, ax = plt.subplots(figsize=(12, 6))
cmap = plt.get_cmap('jet_r')
n_lines = len(Lambda_Reps_filtered)

for i in range(n_lines):
    color = cmap(i / n_lines)
    # Elastic line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_DetectorPos_Lambda_i_filtered[i]*1e-3],
            y_vals,
            color=color,
            linewidth=2,
            label=f'Ei: {wavelength_to_energy(Lambda_Reps_filtered[i]):.2f} meV')
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_DetectorPos_Lambda_i_filtered[i]*1e-3+71],
            y_vals,
            color=color,
            linewidth=2)
    # Inelastic (loss) line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_RRM_Loss_filtered[i]*1e-3],
            y_vals,
            color=color,
            linestyle='dashed',
            linewidth=1.2)
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_RRM_Loss_filtered[i]*1e-3+71],
            y_vals,
            color=color,
            linestyle='dashed',
            linewidth=1.2)
    # Inelastic (gain) line
    ax.plot([ToF_SamplePos_filtered[i]*1e-3, ToF_RRM_Gain_filtered[i]*1e-3],
            y_vals,
            color=color,
            linestyle='dashed',
            linewidth=1.2)
    ax.plot([ToF_SamplePos_filtered[i]*1e-3+71, ToF_RRM_Gain_filtered[i]*1e-3+71],
            y_vals,
            color=color,
            linestyle='dashed',
            linewidth=1.2)

# Draw vertical lines for frame periods
ax.vlines(71, 0, 3, 'k', linestyle='solid')
ax.vlines(71*2, 0, 3, 'k', linestyle='solid')
ax.vlines(71*3, 0, 3, 'k', linestyle='solid')
ax.vlines(71*4, 0, 3, 'k', linestyle='solid')
ax.vlines(71*5, 0, 3, 'k', linestyle='solid')

# Plot settings
ax.set_xlabel("Time-of-Flight (ms)")
ax.set_ylabel("Distance (m)")
ax.legend()
ax.set_title(f"Time–Distance Diagram After Subpulse Suppression - {central_wavelength}Å")
plt.xlim(ToF_RRM_Gain_filtered[0]*1e-3-5, ToF_RRM_Loss_filtered[-1]*1e-3+75)
plt.grid(True)
plt.tight_layout()
plt.show()