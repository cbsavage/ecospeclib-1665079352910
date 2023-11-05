#! /usr/bin/env python

# ***********************************************
# Author: Kenneth Gordon

# Date  :  5 Mar 2020
#
# Research: Create dap.in files for all pixels in
#           a large data grid
# ***********************************************
import os                                    # load os module
import numpy as np                           # load array module
import matplotlib.pyplot as plt              # load plotting module
import time                                  # load time module


# Load the unique pixels, which provides input data on surfaces, cloud
# optical thickness, and cloud top pressure
uniq = np.loadtxt('unique_pixels.dat', delimiter=',')

# Initiate loop through the unique pixels
for i in range(len(uniq)):
    pixel = uniq[i]
    
    land  = pixel[0]
    cot   = pixel[1]
    layer = pixel[2]
    
    
    # Initialize the dap.in file
    f = open("dap.in", "w+")
    
    # Start with the first four lines, which will be the same for every file
    f.write('INPUT FOR RADIATIVE TRANSFER PROGRAM DAP')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nGENERAL PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    
    # GENERAL PARAMETERS
    # Output file based on the location of the pixel in the grid
    f.write('\nOUT_Bc_surf%.0f_cot%.2f_layer%.0f                     :name of output file' % (land, cot, layer))
    
    # Stokes parameters (keep set to 3 for just P_s calculations)
    f.write('\n3                             :size of Stokes parameters (1, 3, or 4)')
    
    # Wavelength file
    f.write('\nwavs.dat                     :file with wavelength data')
    
    # Spectral resolution
    f.write('\n1                            :resolution of spectra (in lines of the wavelength data file)')
    
    # Absorption file
    f.write('\nO2_earth.dat                  :file with absorption data\n')
    
    
    # SURFACE PARAMETERS
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nSURFACE PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    # Wavelength dependence (keep set to 1)
    f.write('\n1                :wavelength dependent surface albedo (0=no,1=yes)')
    
    # Albedo
    f.write('\n0.00             :surface albedo')
    
    # Surface file
    # 0 = ocean
    # 1 = forest
    # 2 = grass/small vegetation
    # 3 = sand/desert
    # 4 = snow/ice
    if land == 0:
        f.write('\nocean.dat         :file with wavelength dependent albedo')
    elif land == 1:
        f.write('\nforest.dat         :file with wavelength dependent albedo')
    elif land == 2:
        f.write('\ngrass.dat         :file with wavelength dependent albedo')
    elif land == 3:
        f.write('\nsand.dat         :file with wavelength dependent albedo')
    elif land == 4:
        f.write('\nice.dat         :file with wavelength dependent albedo')
    
    # Fresnel surface
    f.write('\n0             :Fresnel interface (0=no,1=yes)\n')
    
    
    # GEOMETRIC PARAMETERS
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nGEOMETRIC PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    # Gauss points
    f.write('\n60                :number of Gauss points')
    
    # Gauss point interpolation
    f.write('\n1                 :interpolation in Gauss points (0=no,1=yes)')
    
    # Phase angles
    f.write('\n181                :number of phase angles (max. 180)\n')
    
    
    # AEROSOL PARAMETERS
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nAEROSOL PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    # Number of aerosol types
    f.write('\n1                :number of aerosol types (max. 6)')
    
    # Aerosol scattering file (which clouds to use)
    f.write('\naerBcorr.dat         :file with aerosol scattering coefficients files\n')
    
    
    # LAYER PARAMETERS
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nLAYER PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    # Wavelength at which the cloud optical thickness is specified
    if land == 0:
        f.write('\n0.858           :wavelength at which baer is specified (in microns)')
    elif land == 1 or land == 2 or land == 3:
        f.write('\n0.645           :wavelength at which baer is specified (in microns)')
    elif land == 4:
        f.write('\n1.240           :wavelength at which baer is specified (in microns)')
    
    # Number of levels
    f.write('\n16              :number of levels (from bottom to top)')
    
    # Reflecting surface
    f.write('\n1               :level of reflecting surface\n')
    
    
    # LAYERS TABLE
    
    # Column titles
    f.write('\nn   baer  aaer g1   g2   f1   nr  pres(b)  temp(K) O3 (ppm)    H2O (ppm)')
    f.write('\n------------------------------------------------------------------------')
    
    # Define the standard layers for a clear atmosphere, with no clouds included
    layer1  = '1  0.00  0.0  0.0  0.0  0.0  1  1.0130000 294.0  0.3041E-01  0.1890E+05'
    layer2  = '2  0.00  0.0  0.0  0.0  0.0  1  0.9020000 290.0  0.3353E-01  0.1384E+05'
    layer3  = '3  0.00  0.0  0.0  0.0  0.0  1  0.7100000 279.0  0.4219E-01  0.5990E+04'
    layer4  = '4  0.00  0.0  0.0  0.0  0.0  1  0.6280000 273.0  0.4830E-01  0.3820E+04'
    layer5  = '5  0.00  0.0  0.0  0.0  0.0  1  0.5910000 270.0  0.5930E-01  0.2954E+04'
    layer6  = '6  0.00  0.0  0.0  0.0  0.0  1  0.5540000 267.0  0.7346E-01  0.2223E+04'
    layer7  = '7  0.00  0.0  0.0  0.0  0.0  2  0.4110000 253.0  0.7996E-01  0.1052E+04'
    layer8  = '8  0.00  0.0  0.0  0.0  0.0  2  0.3720000 248.0  0.9122E-01  0.6466E+04'
    layer9  = '9  0.00  0.0  0.0  0.0  0.0  2  0.3320000 243.0  1.0903E-01  0.4057E+04'
    layer10 = '10  0.00  0.0  0.0  0.0  0.0  2  0.2990000 238.0  1.2133E-01  0.3455E+04'
    layer11 = '11  0.00  0.0  0.0  0.0  0.0  2  0.2680000 233.0  1.5059E-01  0.0883E+04'
    layer12 = '12  0.00  0.0  0.0  0.0  0.0  2  0.2380000 228.0  1.8252E-01  0.0974E+04'
    layer13 = '13  0.00  0.0  0.0  0.0  0.0  2  0.2130000 223.0  2.1761E-01  0.00029E+04'
    layer14 = '14  0.00  0.0  0.0  0.0  0.0  2  0.1890000 218.0  0.6793E-01  0.00024E+04'
    layer15 = '15  0.00  0.0  0.0  0.0  0.0  2  0.0250000 213.0  0.0002E-01  0.3935E-03'
    layer16 = '16  0.00  0.0  0.0  0.0  0.0  1  0.0000003 210.0  0.1212E-04  0.3233E-02'
    
    # Print the layers to the file
    if layer == 1:
        f.write('\n1  %.2f  0.0  0.0  0.0  0.0  1  1.0130000 294.0  0.3041E-01  0.1890E+05' % (cot))
    else:
        f.write('\n%s' % (layer1))
        
    if layer == 2:
        f.write('\n2  %.2f  0.0  0.0  0.0  0.0  1  0.9020000 290.0  0.3353E-01  0.1384E+05' % (cot))
    else:
        f.write('\n%s' % (layer2))
        
    if layer == 3:
        f.write('\n3  %.2f  0.0  0.0  0.0  0.0  1  0.7100000 279.0  0.4219E-01  0.5990E+04' % (cot))
    else:
        f.write('\n%s' % (layer3))
        
    if layer == 4:
        f.write('\n4  %.2f  0.0  0.0  0.0  0.0  1  0.6280000 273.0  0.4830E-01  0.3820E+04' % (cot))
    else:
        f.write('\n%s' % (layer4))
    f.write('\n%s' % (layer5))
    f.write('\n%s' % (layer6))
    f.write('\n%s' % (layer7))
    f.write('\n%s' % (layer8))
    f.write('\n%s' % (layer9))
    f.write('\n%s' % (layer10))
    f.write('\n%s' % (layer11))
    f.write('\n%s' % (layer12))
    f.write('\n%s' % (layer13))
    f.write('\n%s' % (layer14))
    f.write('\n%s' % (layer15))
    f.write('\n%s' % (layer16))
    
    # Close the file
    f.close()
    
    # Execute the DAP Code and submit a slurm batch job for each input file
    os.system("sbatch runscript.slurm")
    
    # Wait approximately 1 minute between each batch job submission, to allow initialization
    time.sleep(60)

