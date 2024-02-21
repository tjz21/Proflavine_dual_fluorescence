import h5py
import numpy as np
import cmath
import math
from scipy import integrate
import sys

filename = str(sys.argv[1])
total_E=0.111 # Ha


def full_spectrum_integrant(response_func,E_val):
    integrant=np.zeros(response_func.shape[0])
    counter=0
    while counter<integrant.shape[0]:
        integrant[counter]=(response_func[counter,1]*cmath.exp(1j*response_func[counter,0]*E_val)).real
        counter=counter+1
    return integrant

def full_spectrum(response_func,steps_spectrum,start_val,end_val):
    spectrum=np.zeros((steps_spectrum,2))
    counter=0

    step_length=((end_val-start_val)/steps_spectrum)
    while counter<spectrum.shape[0]:
        E_val=start_val+counter*step_length
        prefac=1.0
        integrant=full_spectrum_integrant(response_func,E_val)
        spectrum[counter,0]=E_val
        spectrum[counter,1]=prefac*(integrate.simps(integrant,dx=response_func[1,0].real-response_func[0,0].real))
        counter=counter+1

    np.savetxt('full_emission_spectrum_proflavine_CAMB3LYP.dat',spectrum)


with h5py.File(filename, "r") as f:
    print("Keys: %s" % f.keys())
    conv_key = list(f.keys())[2]
    print(type(f[conv_key]))

    group=f[conv_key]
    for key in group.keys():
        print(key)
    dipole_real = group["dcf-re"][()]
    dipole_im = group["dcf-im"][()]
    time=group["times"][()]

    # construct response function.

    print(time)
    print(time.shape[0])
    response_func=np.zeros((time.shape[0],5))

    for i in range(response_func.shape[0]):
        response_func[i,0]=time[i]
        response_func[i,1]=dipole_real[i]
        response_func[i,2]=dipole_im[i]
        eff_cmplx=cmath.polar(dipole_real[i]+1j*dipole_im[i])
        response_func[i,3]=eff_cmplx[0]
        response_func[i,4]=eff_cmplx[1]


    # also save S1 pop:
    s1_pop=group["s1"][()]
    s2_pop=group["s2"][()]
    print(s1_pop)

    s1s2_func=np.zeros((time.shape[0],3))
    for i in range(s1s2_func.shape[0]):
        s1s2_func[i,0]=time[i]
        s1s2_func[i,1]=s1_pop[i]
        s1s2_func[i,2]=s2_pop[i]

    np.savetxt('populations.txt',s1s2_func)

    # now make sure the phase is a continuous function
    counter = 0
    phase_fac = 0.0
    while counter < response_func.shape[0] - 1:
        response_func[counter, 4] = response_func[counter, 4] + phase_fac
        if (
            abs(response_func[counter, 4] - phase_fac - response_func[counter + 1, 4]) > 0.7 * math.pi
        ):  # check for discontinuous jump.
            diff = response_func[counter + 1, 4] - (response_func[counter, 4] - phase_fac)
            frac = diff / math.pi
            n = int(round(frac))
            phase_fac = phase_fac - math.pi * n
        #response_func[counter - 1, 4] = response_func[counter - 1, 4] + phase_fac

        counter = counter + 1

    # now apply energy shift to phase:
    for i in range(response_func.shape[0]):
        response_func[i,4]=response_func[i,4]-total_E*response_func[i,0]   # mean energy shift.

    np.savetxt('dipole_response.txt',response_func)

    eff_response=np.zeros((response_func.shape[0],2),dtype=np.complex_)
    for i in range(eff_response.shape[0]):
        eff_response[i,0]=response_func[i,0]
        eff_response[i,1]=response_func[i,3]*cmath.exp(1j*response_func[i,4])

    spectrum_start=0.07
    spectrum_end=0.15

    num_points=2000

    full_spectrum(eff_response,num_points,spectrum_start,spectrum_end)
