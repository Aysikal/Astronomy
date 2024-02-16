#The following file provides a masterflat and gain table FITS file.
import numpy as np
from astropy.io import fits
import os
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import reza
#flats directory
directory = r"C:/Users\AYSAN\Desktop/test"

# get the relevant masterdark:


number_of_flats = len(os.listdir(directory))
flat = np.zeros((1024 , 1024) , dtype = float)
flat_stds = np.zeros((number_of_flats), float)
flat_meds = np.zeros((number_of_flats), float)
flat_exps = np.zeros((number_of_flats), float)
i = 0

for file in os.listdir(directory):
   filepath = os.path.join(directory, file)
        # Open the fits file
   file = fits.open(filepath)
        # Get the data from the first extension
   data = file[0].data
   header = file[0].header
   flat += file[0].data
   flat_stds[i] = np.std(file[0].data)
   flat_meds[i] = np.median(file[0].data)
   flat_exps[i] = float(header['EXPTIME']) * 1e-5 # exposure in seconds
   file.close()

master_flat = flat/(number_of_flats) 
masterflat = reza.sigma_clip(master_flat,3 , 3)

#get gain table by dividing masterflat by its median value
gain_table = masterflat/np.median(masterflat)

# get masterflat and gain table as fits files:

flat_array = np.array(masterflat)
gain_array = np.array(gain_table)

plt.imshow(gain_array)
plt.title("gaintable for red filter" )
plt.colorbar()
plt.show()

plt.imshow(flat_array)
plt.title("masterflat for green filter" )
plt.colorbar()
plt.show()

#export masterflat and gaintable as fits files:
# Define the output FITS file name
output_filename1 = 'masterflat+color.fits'
output_filename2 = 'gaintable+color.fits

# Create a PrimaryHDU (header/data unit) from your array
primary_hdu1 = fits.PrimaryHDU(data=flat_array)
primary_hdu2 = fits.PrimaryHDU(data=gain_array)

# Create an HDUList and append the PrimaryHDU
hdul1 = fits.HDUList([primary_hdu1])
hdul2 = fits.HDUList([primary_hdu2])

# Write the HDUList to the FITS file
hdul1.writeto(output_filename1, overwrite=True)
hdul2.writeto(output_filename2, overwrite=True)
