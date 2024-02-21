#! /usr/bin/env python

# compute_diabatic_sds.py
#  
# Routine takes as an input diabatic energies, oscillator strengths and diabatic coupling between two excited 
# states and computes four spectral densities. S1, S2, S1-S2 cross correlation and S1-S2 coupling.
#
# 1 required command-line argument:
#   $1: input file name containing the fluctuations and couplings (a 5 column file is expected,
#       column 1: S1 energies in eV, column 2: S1 oscillator Str., Column 3: S2 energies in eV, 
#       column 4: S2 oscillator Str., Column 5: S1-S2 coupling in eV
#   3 parameters are specified in the main function that might need to be adjusted by the user: 
#       time_step_in_fs=float (default: 2.0) the timestep between individual data points in the energy gap fluctuations provided as input
#       T=float (default 300.0) the temperature at which the original classical MD that produced the data was run
#       decay_constant=float (default=500.0): A decay constant (in fs) for a decaying exponential that is applied to all autocorrelation
#           functions in the time domain. This guarantees well-behaved Fourier transforms by forcing the correlation function to go to
#           zero in the long timescale limit. It is equivalent to a Lorentzian broadening applied to the entire range of the spectral
#           density. 
#
# Last updated 20240220 by Tim J. Zuehlsdorff
# Based on routines from the MolSpeckPy package by Tim J. Zuehlsorff, available on GitHub:
# https://github.com/tjz21/Spectroscopy_python_code


import os.path
import numpy as np
import math
import sys


# global constants 
kb_in_eV=8.617330350*10.0**(-5.0)
Ha_to_eV=27.211399
fs_to_Ha=2.418884326505*10.0**(-2.0)
hbar_in_eVfs=0.6582119514

# function computes a spectral density from by fourier transforming a correlation function from a time
# series and applying the harmonic  quantum correction factor.
def compute_spectral_dens(corr_func,kbT, sample_rate,time_step):
        # fourier transform correlation func and apply the harmonic quantum correction factor
        corr_freq=time_step*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(corr_func)))
        freqs=np.fft.fftshift(np.fft.fftfreq(corr_func.size,d=1.0/sample_rate))/Ha_to_eV

        spectral_dens=np.zeros((int((corr_freq.shape[-1]+1)/2),2))
        freqs=np.fft.fftshift(np.fft.fftfreq(corr_func.size,d=1.0/sample_rate))/Ha_to_eV

        counter=0
        shift_index=corr_freq.shape[-1]-spectral_dens.shape[0]
        while counter<spectral_dens.shape[0]:
                spectral_dens[counter,0]=freqs[counter+shift_index]
                spectral_dens[counter,1]=freqs[counter+shift_index]/(2.0*kbT)*corr_freq[counter+shift_index].real
                counter=counter+1

        return spectral_dens


def main(filename, time_step_in_fs=2.0, T=300.0, decay_constant=500.0):
    ''' Extract input parameters of S1, S2 and coupling energy gap fluctuations. Find the mean, build correlation
    functions and then finally construct spectral densities for S1, S2, cross correlations and the coupling. The 
    spectral densities are stored as text files with the first column denoting the frequency in Ha, the second column
    the intensity in Ha. 
    '''
    # read user input 
    print(
        f'Generating spectral densities:\n'
        f'Input file name: {filename}\n'
        f'Fluctuation timestep in fs: {time_step_in_fs}\n'
        f'Temperature of the underlying MD simulation: {T}\n'
        f'decay constant (in fs) applied to all correlation funcs: {decay_constant}'
        )

    # some derived constants
    tau=decay_constant/fs_to_Ha
    time_step=time_step_in_fs/fs_to_Ha   # time step of the fluctuations in a.u
    kbT=kb_in_eV*T/Ha_to_eV 
    sampling_rate_in_fs=1.0/time_step_in_fs # one vertical excitation energy calculated every 2 fs. 
    sample_rate=sampling_rate_in_fs*hbar_in_eVfs*math.pi*2.0

    # process inputs. Input file is 5 columns: 
    # S1 energy, S1 oscillator strength, S2 energy, S2 oscillator strength, S1/S2 coupling
    # column 1, 3, and 5 is expected in units of eV
    input_vals=np.genfromtxt(filename)
    energy_S1=input_vals[:,0]/Ha_to_eV
    energy_S2=input_vals[:,2]/Ha_to_eV
    coupling=input_vals[:,4]/Ha_to_eV
    oscillator_S1=input_vals[:,1]
    oscillator_S2=input_vals[:,3]

    # compute mean of fluctuations 
    mean_coupling=np.mean(coupling)
    mean_S1=np.mean(energy_S1)
    mean_S2=np.mean(energy_S2)
    print('Mean S1, Mean S2, mean coupling: \n',
        mean_S1,mean_S2,mean_coupling)

    energy_S1=energy_S1-mean_S1
    energy_S2=energy_S2-mean_S2
    coupling=coupling-mean_coupling

    # compute correlation function using numpy.correlate()
    c_S1=np.correlate(energy_S1,energy_S1,mode='full')
    c_S2=np.correlate(energy_S2,energy_S2,mode='full')
    c_S1_S2=np.correlate(energy_S1,energy_S2,mode='full')
    c_coupling=np.correlate(coupling,coupling,mode='full')

    # multiply correlation func by decaying exponential to guarantee well defined fourier transforms:
    eff_decay_length=tau/time_step
    current_index=-(c_S1.shape[-1]-1)/2*1.0
    counter=0
    while counter<c_S1.shape[-1]:
        c_S1[counter]=c_S1[counter]*math.exp(-abs(current_index)/eff_decay_length)/(c_S1.shape[-1]/2.0)
        c_S2[counter]=c_S2[counter]*math.exp(-abs(current_index)/eff_decay_length)/(c_S1.shape[-1]/2.0)
        c_S1_S2[counter]=c_S1_S2[counter]*math.exp(-abs(current_index)/eff_decay_length)/(c_S1_S2.shape[-1]/2.0)
        c_coupling[counter]=c_coupling[counter]*math.exp(-abs(current_index)/eff_decay_length)/(c_coupling.shape[-1]/2.0)
        current_index=current_index+1.0
        counter=counter+1

    # build SDs
    sd_s1=compute_spectral_dens(c_S1,kbT, sample_rate,time_step)
    sd_s2=compute_spectral_dens(c_S2,kbT, sample_rate,time_step)
    sd_coupling=compute_spectral_dens(c_coupling,kbT, sample_rate,time_step)
    sd_s1_s2=compute_spectral_dens(c_S1_S2,kbT, sample_rate,time_step)

    # store them as text files
    np.savetxt('spectral_dens_diabatic_S1.txt', sd_s1)
    np.savetxt('spectral_dens_diabatic_S2.txt', sd_s2)
    np.savetxt('spectral_dens_diabatic_cross_corr.txt', sd_s1_s2)
    np.savetxt('spectral_dens_diabatic_coupling.txt', sd_coupling)

    print('Done!')
    return [
            sd_s1, sd_s2, sd_s1_s2, sd_coupling
            ]

# Run `main` immediately if called as a script, but don't when loaded as a
# module
if __name__ == "__main__":
    args = sys.argv[1]
    print(args)
    result = main(args)
