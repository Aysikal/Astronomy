#The following file provides a masterflat and gain table FITS file.
import numpy as np
from astropy.io import fits
import os
from astropy.stats import sigma_clip

#flats dictionary
directory = r"flat FITS address"

# get relevant masterdark
masterdark = fits.open(r"co-responding masterdark address")
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
        #if you have a specific filter on your images, simply adding the channels will suffice.
        grey_minus_dark = data[0,:,:]+data[1,:,:]+data[2,:,:] - dark_data
        #if not, use this ratio instead : [0.2989, 0.5870, 0.1140]
        #grey_minus_dark = np.dot(data.T, [0.2989, 0.5870, 0.1140]) - dark_data
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

