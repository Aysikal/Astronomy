#This code was written by Aysan Hemmati
import numpy as np
from astropy.io import fits
import os
from astropy.stats import sigma_clip

#flats dictionary (4s)
directory = r"C:\Users\aysan\Desktop\university\Asto Lab\Aznaveh data\flat\flat fits"

# get relevant masterdark (4s)
masterdark = fits.open(r"C:\Users\aysan\Desktop\university\Asto Lab\Aznaveh data\flat\4s masterdark.fits")
dark_data = masterdark[0].data

#import flat images from dictionary and save as grey scale data minus masterdark
grey_data_minus_dark = []

for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # Open the fits file
        file = fits.open(filepath)
        # Get the data from the first extension
        data = file[0].data
        # Summing over 3 channels to get grey scale data, subtract masterdark 
        grey_minus_dark = data[0,:,:]+data[1,:,:]+data[2,:,:] - dark_data
        grey_data_minus_dark.append(grey_minus_dark)

#sigma clip stacked images
images = np.array(grey_data_minus_dark)
stacked_images = np.stack(images)
data_clipped = sigma_clip(stacked_images, sigma=3, cenfunc='median')

#get masterflat by using median of the sigma clipped data
masterflat = np.array(np.median(data_clipped,axis=0)) 

#get gain table by dividing masterflat by its median value
gain_table = masterflat/np.median(masterflat)

# get masterflat and gain table as fits files

hdu1 = fits.PrimaryHDU(gain_table)
hdulist = fits.HDUList([hdu1])
hdu1.writeto('gain_table.fits')

hdu2 = fits.PrimaryHDU(masterflat)
hdulist = fits.HDUList([hdu1])
hdu2.writeto('masterflat.fits')

