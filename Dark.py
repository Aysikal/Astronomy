#this code was written by Aysan Hemmati
import numpy as np
from astropy.io import fits
import os
from astropy.stats import sigma_clip

directory = r"C:\Users\aysan\Desktop/university\Asto Lab\Aznaveh data\dark\30s fits"

grey_data = []

#get the dark data from dictionary and save as grey scale data

for filename in os.listdir(directory):
        
        filepath = os.path.join(directory, filename)
        # Open the fits file
        file = fits.open(filepath)
        # Get the data from the first extension
        data = file[0].data
        # Summing over 3 channels
        grey_scale = data[0,:,:]+data[1,:,:]+data[2,:,:]
        grey_data.append(grey_scale)

#get the grey data as a stacked image and sigmaclip the data
images = np.array(grey_data)
stacked_images = np.stack(images)
data_clipped = sigma_clip(stacked_images, sigma=3, cenfunc='median')
median =np.array(np.median(data_clipped,axis=0)) 

# get masterdark as a fits file
hdu = fits.PrimaryHDU(median)
hdulist = fits.HDUList([hdu])
hdu.writeto('30s masterdark.fits')