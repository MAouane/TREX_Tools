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
L_PS = 108.0 #Source to PS choppers in m
LM = 162.0  #Distance to the M-chopper
T_offset = 1700

#Conversion functions
def energy_to_wavelength(E):
    return np.sqrt(ENERGY_CONSTANT / E)

def wavelength_to_energy(lam):
    return ENERGY_CONSTANT / lam**2

def tof_to_lambda(tof):
    return tof / (TOF_CONSTANT * L_sample)

def lambda_to_tof_PS(lam):
    return T_offset + TOF_CONSTANT * lam * L_PS

def lambda_to_tof_fan(lam):
    return T_offset + TOF_CONSTANT * lam * L_fan

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



central_wavelength = [2.4] #Central wavelength for your setting

#if central_wavelength > 5.0 or central_wavelength < 0.5:
#    print('-------------------------------------------------------------------------')
#    print('Well... Maybe not?') 
#else:
#    print('-------------------------------------------------------------------------')
#    print(f'Your central wavelength is: {central_wavelength} Å')
    
# Rounding all of the floating points to .2f format to avoid any and all floating errors... Not setting the .2f messes up the subpulses for some reason
delta_lambda = round((3956)*(1/(14*166.8)), 2)    # Total wavelength bandwidth of T-REX, should be in the 1.7Å range
print(f'Natural bandwidth of the instrument is: {delta_lambda}Å')
print(f'Wavelength range is delimited by: +/-{delta_lambda/2.}Å')

fM = 252-14*2 # M-Chopper frequency in Hz (must be multiple of 14 between 14–336). The multiple defines how many reps you get on the sample
fP = fM*0.75 # M-Chopper frequency in Hz. Added as part of the contingency


frame_period_us = 71000  # µs (71 ms) - One ESS pulse duration. To check if longer wavelengths cross into the the fast neutron of the following frame

if fM % 14 != 0 or not (14 <= fM <= 336):
    raise ValueError("fM (M-chopper frequency) must be a multiple of 14 between 14 and 336Hz")
'''
if fP > 140:
    raise ValueError(f"fP (PS-chopper frequency) must be below 140Hz for the Al disks. Your frequency of {fP}Hz exceeds it")
''' 
for j in range(len(central_wavelength)):
#    fig, ax = plt.subplots(figsize=(12, 6))
    fig, ax = plt.subplots(figsize=(15, 8))

    print(f'For your chosen frequency of {fM}Hz, you get {int(fM/14)} Ei Reps')
    print('-------------------------------------------------------------------------')
   # print(f'For your chosen frequency of {fM}Hz, Tthe PS chopper is running at: {fP}Hz')
   # print('-------------------------------------------------------------------------')

    d_lambda = (H_OVER_MN / (fM * LM)) #Step size for subpulses
    print(f'Step size for this setting is: {d_lambda}Å')
    lambda_min = central_wavelength[j] - delta_lambda / 2.
    lambda_max = central_wavelength[j] + delta_lambda / 2.

    print(f'Your limits are {lambda_min}Å and {lambda_max}Å')
    print('-------------------------------------------------------------------------')

    reps = int(fM/14)
    Lambda_Reps = generate_wavelengths(central_wavelength[j],fM)
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

    ToF_PSPos = lambda_to_tof_PS(Lambda_Reps)
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
    Gain_fraction = 0.50  # 50% gain in neutron energy for both thermal and cold wavelengths

    ''' Originally, I set all cold wavelength to only have 20% of Loss since you'd only look at low energy transfers/QENS but changed it to be 85% like the thermal wavelengths '''

    for i in range(len(Lambda_Reps)):    
        
        if Lambda_Reps[i] < 10.0:
            ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.20) * Lambda_Reps[i] * L_detector
        else:
            ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.80) * Lambda_Reps[i] * L_detector
        
        ToF_RRM_Gain_Detector[i] = ToF_SamplePos[i].copy() #TOF_CONSTANT * Lambda_Reps[i] * (np.sqrt(1. / (1. + Gain_fraction))) * L_detector 
        ToF_RRM_Gain[i] = ToF_SamplePos[i].copy() #ToF_RRM_Gain_Detector[i] + ToF_SamplePos[i]
        
      #' if Lambda_Reps[i] < 2.5:
      #      ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.20) * Lambda_Reps[i] * L_detector
      #  else:
      #      ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.50) * Lambda_Reps[i] * L_detector
        
       # ToF_RRM_Gain_Detector[i] = TOF_CONSTANT * Lambda_Reps[i] * (np.sqrt(1. / (1. + Gain_fraction))) * L_detector 
       # ToF_RRM_Gain[i] = ToF_RRM_Gain_Detector[i] + ToF_SamplePos[i]

       # else:
       #     ToF_RRM_Gain_Detector[i] = TOF_CONSTANT * Lambda_Reps[i] * (np.sqrt(1. / (1. + Gain_fraction))) * L_detector 
       #     ToF_RRM_Gain[i] = ToF_RRM_Gain_Detector[i] + ToF_SamplePos[i]
            
        
        #ToF_RRM_Loss_Detector[i] =  TOF_CONSTANT * np.sqrt(1./0.80) * Lambda_Reps[i] * L_detector 
        
        ToF_RRM_Loss[i] = ToF_RRM_Loss_Detector[i] + ToF_SamplePos[i]

        Lambda_RRM_Loss = tof_to_lambda(ToF_RRM_Loss)
        E_RRM_Loss = wavelength_to_energy(Lambda_RRM_Loss)

    print('To compare the ToFs with McStas, the same offset was added to all ToFs')
    print('-------------------------------------------------------------------------')
    print('ToFs at the PS chopper position for this setting in s are: \n', ToF_PSPos*1e-6)
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

    cmap = plt.get_cmap('jet_r')
    n_lines = len(Lambda_Reps)

    for i in range(n_lines):
        if Lambda_Reps[i] < 0.7:
            color = 'grey'
        else:
            color = cmap(i / n_lines)
        
        #Elastic line
        ax.plot([ToF_SamplePos[i]*1e-3, ToF_DetectorPos_Lambda_i[i]*1e-3], y_vals, color=color, linewidth=2, label = f'λi: {Lambda_Reps[i]:.2f}Å')
#        ax.plot([ToF_SamplePos[i]*1e-3, ToF_DetectorPos_Lambda_i[i]*1e-3], y_vals, color=color, linewidth=2, label = f'Ei: {wavelength_to_energy(Lambda_Reps[i]):.2f}meV')
        ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_DetectorPos_Lambda_i[i]*1e-3+71], y_vals, color=color, linewidth=2)
        
        #Inelastic (loss) line
        ax.plot([ToF_SamplePos[i]*1e-3, ToF_RRM_Loss[i]*1e-3], y_vals, color=color, linestyle='dashed', linewidth=1.2)
        ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_RRM_Loss[i]*1e-3+71], y_vals, color=color, linestyle='dashed', linewidth=1.2)
        
        #Inelastic (loss) line
        ax.plot([ToF_SamplePos[i]*1e-3, ToF_RRM_Gain[i]*1e-3], y_vals, color=color, linestyle='-.', linewidth=1.2)
        ax.plot([ToF_SamplePos[i]*1e-3+71, ToF_RRM_Gain[i]*1e-3+71], y_vals, color=color, linestyle='-.', linewidth=1.2)

    #ax.vlines(71, 0, 3, 'k', linestyle = 'solid')
    #ax.vlines(71*2, 0, 3, 'k', linestyle = 'solid')
    #ax.vlines(71*3, 0, 3, 'k', linestyle = 'solid')
    #ax.vlines(71*4, 0, 3, 'k', linestyle = 'solid')
    #ax.vlines(71*5, 0, 3, 'k', linestyle = 'solid')

    ax.set_xlabel("Time-of-Flight (ms)", fontsize = 15)
    ax.set_ylabel("Distance (m)", fontsize = 15)
    ax.legend(bbox_to_anchor=(1.0, 1.0), fontsize = 15)
    ax.set_title(f"Time-Distance Diagram from Sample to Detector - fM = {fM}Hz - λc = {central_wavelength[j]}Å", fontsize = 15)
#    ax.set_title(f"Time-Distance Diagram from Sample to Detector - Maximum Chopper Speed - λc = {central_wavelength[j]}Å - fM = {fM}Hz", fontsize = 15)
    plt.xlim(ToF_RRM_Gain[0]*1e-3-5, ToF_RRM_Loss[-1]*1e-3+75)
    plt.grid(True)
    plt.minorticks_on()
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.tight_layout()
    plt.show()



print('Looking at the frame overlap between subsequent wavelengths in RRM.')
print('-------------------------------------------------------')


for i in range(len(ToF_DetectorPos_Lambda_i)):
    i_next = (i + 1) % len(ToF_DetectorPos_Lambda_i)  #Avoid running error of 'out of range'

    if i == len(ToF_DetectorPos_Lambda_i) - 1:
        x_next = ToF_DetectorPos_Lambda_i[0] + 71*1e3  #Shift to the second ESS pulse for comparison of frame overlap
        gain_next = ToF_RRM_Gain[0] + 71*1e3
        lambda_next = Lambda_Reps[0]
    else:
        x_next = ToF_DetectorPos_Lambda_i[i_next]
        gain_next = ToF_RRM_Gain[i_next]
        lambda_next = Lambda_Reps[i_next]

    print(f'Wavelength assessed: {Lambda_Reps[i]:.2f}Å or Ei = {wavelength_to_energy(Lambda_Reps[i]):.2f}meV')

    x_el = ToF_DetectorPos_Lambda_i[i]
    x_loss = ToF_RRM_Loss[i]
    
    if Lambda_Reps[i] < 0.7:
        print(f'{Lambda_Reps[i]:.2f}Å is outside the ESS range... ')
        
    if gain_next > x_loss and x_next > x_loss:
        print('Clean Rep! Use it!')
        print('-------------------------------------------------------')
    elif x_el <= x_next <= x_loss:
        print('Overlap detected!')
        print(f'Elastic line of {lambda_next:.2f}Å is inside the loss window of {Lambda_Reps[i]:.2f}Å. Impossible to use')
        print('-------------------------------------------------------')
    else:
        print('Overlap detected!')
        loss_window = np.abs(x_el - x_loss)
        fo_range = np.abs(x_loss - gain_next)
        clean_window = np.abs(loss_window - fo_range)
        ratio_loss = (clean_window/loss_window)
        print(f'Lost dynamical range for Ei = {wavelength_to_energy(Lambda_Reps[i]):.2f} is {(100.*(1-ratio_loss)):.2f}%')
        #print(f'Maximum amount of energy loss observable cleanly: {(80-(1-ratio_loss)*80):.2f}')
        print('-------------------------------------------------------')

# ===== Usable LOSS window for each wavelength =====
usable_loss_fraction = []
for i in range(len(ToF_DetectorPos_Lambda_i)):
    i_next = (i + 1) % len(ToF_DetectorPos_Lambda_i)
    if i == len(ToF_DetectorPos_Lambda_i) - 1:
        x_next = ToF_DetectorPos_Lambda_i[0] + 71*1e3
        gain_next = ToF_RRM_Gain[0] + 71*1e3
    else:
        x_next = ToF_DetectorPos_Lambda_i[i_next]
        gain_next = ToF_RRM_Gain[i_next]
    x_el = ToF_DetectorPos_Lambda_i[i]
    x_loss = ToF_RRM_Loss[i]
    if Lambda_Reps[i] < 0.7:
        usable_loss_fraction.append(0.0)
        continue
    if gain_next > x_loss and x_next > x_loss:
        usable_loss_fraction.append(1.0)
        continue
    elif x_el <= x_next <= x_loss:
        usable_loss_fraction.append(0.0)
        continue
    else:
        loss_window = np.abs(x_el - x_loss)
        fo_range = np.abs(x_loss - gain_next)
        clean_window = np.abs(loss_window - fo_range)
        frac = max(0.0, min(1.0, clean_window / loss_window))
        usable_loss_fraction.append(frac)

# ===== Usable GAIN window for each wavelength =====
usable_gain_fraction = []
for i in range(len(ToF_DetectorPos_Lambda_i)):
    i_prev = (i - 1) % len(ToF_DetectorPos_Lambda_i)
    if i == 0:
        x_prev_el = ToF_DetectorPos_Lambda_i[-1] - 71*1e3
        loss_prev = ToF_RRM_Loss[-1] - 71*1e3
    else:
        x_prev_el = ToF_DetectorPos_Lambda_i[i_prev]
        loss_prev = ToF_RRM_Loss[i_prev]
    gain_this = ToF_RRM_Gain[i]
    x_el = ToF_DetectorPos_Lambda_i[i]
    if Lambda_Reps[i] < 0.7:
        usable_gain_fraction.append(0.0)
        continue
    # Gain window: gain_this to x_el (gain_this < x_el)
    # Loss_prev overlapping with this window reduces fraction
    if loss_prev < gain_this:
        usable_gain_fraction.append(1.0)
    elif loss_prev > x_el:
        usable_gain_fraction.append(1.0)
    elif gain_this <= loss_prev <= x_el:
        gain_window = np.abs(x_el - gain_this)
        overlap = np.abs(loss_prev - gain_this)
        clean_window = gain_window - overlap
        frac = max(0.0, min(1.0, clean_window / gain_window))
        usable_gain_fraction.append(frac)
    else:
        usable_gain_fraction.append(0.0)

# ===== PLOT BOTH TOGETHER =====
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(Lambda_Reps, usable_loss_fraction, 'o-', color='dodgerblue', linewidth=2, label='Loss window')
ax.plot(Lambda_Reps, usable_gain_fraction, 's-', color='darkorange', linewidth=2, label='Gain window')
ax.set_xlabel('Incident Wavelength (Å)', fontsize=15)
ax.set_ylabel('Fraction of Usable Window', fontsize=15)
ax.set_title(f'Usable Energy-Loss and Energy-Gain Windows\n Central Wavelength = {central_wavelength[0]}Å, Chopper Speed: {fM}Hz', fontsize=15)
ax.set_ylim(-0.05, 1.05)
ax.grid(True)
ax.legend(fontsize=13)
ax.tick_params(labelsize=13)
plt.tight_layout()
plt.show()

