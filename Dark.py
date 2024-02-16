#The following file provides a masterdark FITS file.
import numpy as np
import matplotlib
from astropy.io import fits
import os
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import reza
directory = "enter location to the dark files here"

number_of_darks = len(os.listdir(directory))
dark = np.zeros((1024 , 1024) , dtype = float)
dark_stds = np.zeros((number_of_darks), float)
dark_meds = np.zeros((number_of_darks), float)
dark_exps = np.zeros((number_of_darks), float)
i = 0

for file in os.listdir(directory):
   filepath = os.path.join(directory, file)
        # Open the fits file
   file = fits.open(filepath)
        # Get the data from the first extension
   data = file[0].data
   header = file[0].header
   dark += file[0].data
   dark_stds[i] = np.std(file[0].data)
   dark_meds[i] = np.median(file[0].data)
   dark_exps[i] = float(header['EXPTIME']) * 1e-5 # exposure in seconds
   file.close()

master_dark = dark/(number_of_darks) 
masterdark = sigma_clip(master_dark, sigma=4, cenfunc='median')

#statistics of the master dark (mean, std, max, min):
print(reza.statist(masterdark))
#plot the histogram:
plt.hist(masterdark)
plt.show()
#plot the masterdark array:
plt.imshow(masterdark)
plt.show()

#export the master dark as a fits file:
# Define the output FITS file name
output_filename = 'masterdark+exptime.fits'

# Create a PrimaryHDU (header/data unit) from your array
primary_hdu = fits.PrimaryHDU(masterdark)

# Create an HDUList and append the PrimaryHDU
hdul = fits.HDUList([primary_hdu])

# Write the HDUList to the FITS file
hdul.writeto(output_filename, overwrite=True)

