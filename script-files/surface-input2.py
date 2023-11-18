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
# import time                                  # load time module


# Load the unique pixels, which provides input data on surfaces

# change with listing of folder i'm working with
# list folder name
# once folder is done 

# Open and read test_list.txt:
# ls -1
# pipe into folder
# with open() won't work, listdir
with open( 'test_list.txt', 'r' ) as uniq:
   srfc_list = uniq.read().splitlines()

# Initiate loop through the unique surfaces:   
for srfc in srfc_list:
    print( srfc )
    
    # Initialize the dap.in file
    f = open("dap.in", "w+")

    # Start with the first four lines, which will be the same for every file
    f.write('INPUT FOR RADIATIVE TRANSFER PROGRAM DAP')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nGENERAL PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    
    # GENERAL PARAMETERS
    # Output file based on the location of the pixel in the grid
    f.write('\nOUT_%s_srfc                     :name of output file' % ( srfc ))
    
    # Stokes parameters (keep set to 3 for just P_s calculations)
    f.write('\n3                             :size of Stokes parameters (1, 3, or 4)')
    
    # Wavelength file
    f.write('\nwav_hres.dat                     :file with wavelength data')
    
    # Spectral resolution
    f.write('\n20                            :resolution of spectra (in lines of the wavelength data file)')
    
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
    # add folder i'm using before /{srfc}
    f.write(f'\n{srfc}         :file with wavelength dependent albedo') 
    
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
    f.write('\naerA.dat         :file with aerosol scattering coefficients files\n')
    
    
    # LAYER PARAMETERS
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    f.write('\nLAYER PARAMETERS')
    f.write('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    # Wavelength at which the cloud optical thickness is specified
    f.write('\n0.858           :wavelength at which baer is specified (in microns)')
    
    # Number of levels
    f.write('\n2               :number of levels (from bottom to top)')
    
    # Reflecting surface
    f.write('\n1               :level of reflecting surface\n')
    
    
    
    # Column titles
    f.write('\nn   baer  aaer g1   g2   f1   nr  pres(b)  temp(K) O3 (ppm)    H2O (ppm)')
    f.write('\n------------------------------------------------------------------------')
    
    # Define the standard layers for a clear atmosphere, with no clouds included
    # layer15 = '15  0.00  0.0  0.0  0.0  0.0  1  0.0250000 213.0  0.0002E-01  0.3935E-03'
    # layer16 = '16  0.00  0.0  0.0  0.0  0.0  1  0.0000003 210.0  0.1212E-04  0.3233E-02'
    
    # Define and print the layers to the file
    f.write('\n15  0.0  0.0  0.0  0.0  0.0  1  0.0250000 213.0  0.0002E-01  0.3935E-03')
    f.write('\n16  0.0  0.0  0.0  0.0  0.0  1  0.0000003 210.0  0.1212E-04  0.3233E-02')

    
    # Close the file:
    f.close()
    
    # Execute the DAP Code and submit a slurm batch job for each input file
    os.system("sbatch runscript.slurm")
    
    # Wait approximately 1 minute between each batch job submission, to allow initialization
    # time.sleep(60)


