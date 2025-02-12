# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:11:32 2024

@author: Micha
"""

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d

# Global dictionary to store loaded spectral data
spectral_data = {}
data = loadmat('spectralLIB.mat')
spectral_data['muadeoxy'] = data['muadeoxy'].flatten()
spectral_data['muafat'] = data['muafat'].flatten()
spectral_data['muamel'] = data['muamel'].flatten()
spectral_data['muaoxy'] = data['muaoxy'].flatten()
spectral_data['muawater'] = data['muawater'].flatten()
spectral_data['musp'] = data['musp'].flatten()
spectral_data['nmLIB'] = data['nmLIB'].flatten()

def func_mua_stratum(wavelength):
    B = 0  # Blood content
    S = 0  # Blood oxygen saturation
    W = 0.2  # Water content
    M = 0.01  # Melanin content
    F = 0  # Fat content
    return calc_mua(wavelength, S, B, W, F, M)
def func_mus_stratum(wavelength):
    aPrime = 40  # musPrime at 500 nm
    fRay = 0.4  # Fraction of scattering due to Rayleigh scattering
    bMie = 1  # Scattering power for Mie scattering
    g = 0.9  # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)

"Epidermis = good"
def func_mua_epidermis(wavelength):
    B = 0  # Blood content
    S = 0.39  # Blood oxygen saturation (estimated same as dermis)
    W = 0.70  # Water content
    M = 0.025  # Melanin content
    F = 0  # Fat content
    return calc_mua(wavelength, S, B, W, F, M)
def func_mus_epidermis(wavelength):
    aPrime = 66.7  # musPrime at 500 nm
    fRay = 0.29  # Fraction of scattering due to Rayleigh scattering
    bMie = 0.689  # Scattering power for Mie scattering
    g = 0.9  # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)



    aPrime = 42.4 # musPrime at 500 nm
    fRay = 0.62 # Fraction of scattering due to Rayleigh scattering
    bMie = 1 # Scattering power for Mie scattering
    g = 0.9 # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)

    aPrime = 40  # musPrime at 500 nm
    fRay = 0  # Fraction of scattering due to Rayleigh scattering
    bMie = 1  # Scattering power for Mie scattering
    g = 0.9  # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)
"dermis = good"
def func_mua_dermis(wavelength):
    B = 0.002 # Blood content
    S = 0.74 # Blood oxygen saturation
    W = 0.65 # Water content
    M = 0 # Melanin content
    F = 0 # Fat content
    return calc_mua(wavelength, S, B, W, F, M)
def func_mus_dermis(wavelength):
    aPrime = 43.6 # musPrime at 500 nm
    fRay = 0.41 # Fraction of scattering due to Rayleigh scattering
    bMie = 0.563 # Scattering power for Mie scattering
    g = 0.9 # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)

"need mua for hypodermis"
def func_mua_hypodermis(wavelength):
    B = 0.002 # Blood content
    S = 0.86 # Blood oxygen saturation
    W = 0.216 # Water content
    M = 0 # Melanin content
    F = 0.80 # Fat content
    return calc_mua(wavelength, S, B, W, F, M)
def func_mus_hypodermis(wavelength):
    aPrime = 21.35 # musPrime at 500 nm
    fRay = 0.1475 # Fraction of scattering due to Rayleigh scattering
    bMie = 0.537 # Scattering power for Mie scattering
    g = 0.9 # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)

def func_mua_blood(wavelength):
    B = 1 # Blood content
    S = 0.75 # Blood oxygen saturation
    W = 0.95 # Water content
    M = 0 # Melanin content
    F = 0 # Fat content
    return calc_mua(wavelength, S, B, W, F, M)
def func_mus_blood(wavelength):
    aPrime = 10 # musPrime at 500 nm
    fRay = 0 # Fraction of scattering due to Rayleigh scattering
    bMie = 1 # Scattering power for Mie scattering
    g = 0.9 # Scattering anisotropy
    return calc_mus(wavelength, aPrime, fRay, bMie, g)


def calc_mus(wavelength, aPrime, fRay, bMie, g):
    fMie = 1 - fRay
    musPrime = aPrime * (fRay * (wavelength / 500) ** (-4) + fMie * (wavelength / 500) ** (-bMie))
    return musPrime / (1 - g)
def calc_mua(wavelength, S, B, W, F, M):
    global spectral_data
    mua_deoxy = interp1d(spectral_data['nmLIB'], spectral_data['muadeoxy'])(wavelength)
    mua_fat = interp1d(spectral_data['nmLIB'], spectral_data['muafat'])(wavelength)
    mua_mel = interp1d(spectral_data['nmLIB'], spectral_data['muamel'])(wavelength)
    mua_oxy = interp1d(spectral_data['nmLIB'], spectral_data['muaoxy'])(wavelength)
    mua_water = interp1d(spectral_data['nmLIB'], spectral_data['muawater'])(wavelength)

    mua = B * S * mua_oxy + B * (1 - S) * mua_deoxy + W * mua_water + F * mua_fat + M * mua_mel #no bili and beta carotene
    return mua


#Assuming the required function definitions and imports are already in place
def calc_anisotropy(wavelength):
    g = 0.62 + 0.29 *wavelength* 10**(-3)
    return g
# Set the wavelength
wavelength = 635
g = calc_anisotropy(wavelength)  # Scattering anisotropy factor


# Calculate optical properties
properties = {
    'epidermis': {
        'mu_a': func_mua_epidermis(wavelength) / 10,  # mm^-1
        'mu_s': func_mus_epidermis(wavelength) / 10,  # mm^-1
        'mu_s_prime': func_mus_epidermis(wavelength) / 10 * (1 - g)  # mm^-1
    },
    'dermis': {
        'mu_a': func_mua_dermis(wavelength) / 10,  # mm^-1
        'mu_s': func_mus_dermis(wavelength) / 10,  # mm^-1
        'mu_s_prime': func_mus_dermis(wavelength) / 10 * (1 - g)  # mm^-1
    },
    'hypodermis': {
        'mu_a': func_mua_hypodermis(wavelength) / 10,  # mm^-1
        'mu_s': func_mus_hypodermis(wavelength) / 10,  # mm^-1
        'mu_s_prime': func_mus_hypodermis(wavelength) / 10 * (1 - g)  # mm^-1
    }
}

# Display properties in a readable format using pandas (optional for better visualization)
import pandas as pd

# Create a DataFrame for visualization
df_properties = pd.DataFrame(properties).T  # Transpose for better layout
df_properties.columns = ['μa (mm^-1)', 'μs (mm^-1)', 'μs\' (mm^-1)']

print(df_properties)


def parse_data(file_path):
    data = {}
    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.split()  # Split the line into components based on whitespace
            if len(parts) < 5:
                continue  # Skip lines that do not have enough data
            
            wavelength = int(parts[0])
            data[wavelength] = {
                'mua_98': float(parts[1]),
                'mua_0': float(parts[2]),
                'mus_98': float(parts[3]),
                'g': float(parts[4])
            }
    
    return data

# Example usage
file_path = 'blood_opticalproperties.txt'  # Specify the path to your data file
blood_data_all_wavelengths = parse_data(file_path)
wavelength_data =False

while wavelength_data == False:
    if blood_data_all_wavelengths[wavelength]: 
        
        blood_data = blood_data_all_wavelengths[wavelength]
        wavelength_data = True
        print("blood values", blood_data)
        print(wavelength), 
    else:
            wavelength += 1


