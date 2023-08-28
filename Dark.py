#The following file provides a masterdark FITS file.
import numpy as np
from astropy.io import fits
import os
from astropy.stats import sigma_clip

directory = r"Adress to the raw dark fits files"

grey_data = []

#get the dark data from dictionary and save as grey scale data

for filename in os.listdir(directory):
        
        filepath = os.path.join(directory, filename)
        # Open the fits file
        file = fits.open(filepath)
        # Get the data from the first extension
        data = file[0].data
        
        # Summing over 3 channels, if you have a specific filter on your images, simply adding the channels will suffice.
        grey_scale = data[0,:,:]+data[1,:,:]+data[2,:,:]
        #if not, use this ratio instead : [0.2989, 0.5870, 0.1140]
        # grayscale = np.dot(data.T, [0.2989, 0.5870, 0.1140])
        
        grey_data.append(grey_scale)

#get the grey data as a stacked image and sigmaclip the data
images = np.array(grey_data)
stacked_images = np.stack(images)
data_clipped = sigma_clip(stacked_images, sigma=3, cenfunc='median')
median =np.array(np.median(data_clipped,axis=0)) 

# get masterdark as a fits file
hdu = fits.PrimaryHDU(median)
hdulist = fits.HDUList([hdu])
hdu.writeto('masterdark name.fits')
