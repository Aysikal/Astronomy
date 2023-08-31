# The following code calculates absorption data by fitting a weighted line to the data.
import numpy as np
from astropy.io import fits
import os
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry, CircularAnnulus
from scipy import ndimage
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
  
#Choose a star that appears in all of the images, using Starry night check wether the chosen star is a variable star or not. do not choose a variable star.

#importing dark and flat corrected fits images
directory = r"Address to your images"
images = []
for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # Open the fits file
        file = fits.open(filepath)
        # Get the data from the first extension
        data = file[0].data
        grey_data = data[0,:,:]+data[1,:,:]+data[2,:,:]
        images.append(grey_data)

#location of your chosen star; for example:
star_x = [ 3064, 3144, 3174, 3254, 3284, 3723 ,3429, 3463, 3493, 3553, 3566, 3616, 3633, 3700]
star_y = [ 1437, 1477, 1527, 1573, 1639, 1410 , 801, 825, 871, 928, 392, 442, 495, 562]


def apparent_magnitude(star_flux):
    magnitude =  - 2.5 * np.log10(star_flux)
    return magnitude

def find_optimal_flux(image_data, x, y, max_aperture_radius=16, step_size=0.5):
    y, x = int(round(y)), int(round(x))
    # Initialize variables for the best aperture and its flux
    best_aperture = None
    best_flux = 0.0
    best_snr = 0.0

    aperture_radii = []
    aperture_snr = []
    # Try different aperture sizes and find the one with the highest SNR
    for aperture_radius in np.arange(0.5, max_aperture_radius + step_size, step_size):
        aperture = CircularAperture((x, y), r=aperture_radius)
        phot_table = aperture_photometry(image_data, aperture)
        flux = phot_table['aperture_sum'][0]
       
        # Calculate the noise and SNR for the current aperture
        gain = 1.0
        readnoise = 1.0
        n_pixels = np.pi * aperture_radius ** 2
        noise = np.sqrt(flux * gain + n_pixels * readnoise ** 2)
        snr = flux * gain / noise
        
        aperture_radii.append(aperture_radius)
        # Calculate the inner and outer radii and the sky background flux
        inner_radius = aperture_radius * 12/5
        outer_radius = aperture_radius * 15/5
        sky_aperture = CircularAnnulus((x, y), r_in=inner_radius, r_out=outer_radius)
        sky_table = aperture_photometry(image_data, sky_aperture)
        sky_flux = sky_table['aperture_sum'][0]
        sky_area = np.pi * (outer_radius ** 2 - inner_radius ** 2)
        sky_background_flux = sky_flux / sky_area

        # Subtract the sky background flux from the star flux to get the total flux
        total_flux = flux - (abs(sky_background_flux) * n_pixels)
        # Update the best aperture if the current SNR is higher
        aperture_snr.append(snr)
        aperture_radii.append(aperture_radius)

        if snr > best_snr:
            best_snr = snr
            best_flux = total_flux
            best_aperture = aperture
    return best_flux, best_aperture , aperture_radii, aperture_snr

# a flux function with known radius
def flux(image_data, x, y, radius=12.5):
    y, x = int(round(y)), int(round(x))
    # Initialize variables for the best aperture and its flux
    best_flux = 0.0
    best_snr = 0.0
    aperture = CircularAperture((x, y), r=12.5)
    phot_table = aperture_photometry(image_data, aperture)
    flux = phot_table['aperture_sum'][0]

    # Calculate the noise and SNR for the current aperture
    gain = 1.0
    readnoise = 1.0
    n_pixels = np.pi * radius ** 2
    noise = np.sqrt(flux * gain + n_pixels * readnoise ** 2)
    snr = best_flux * gain / noise 


    # Calculate the inner and outer radii and the sky background flux
    inner_radius = radius * 12/5
    outer_radius = radius * 15/5
    sky_aperture = CircularAnnulus((x, y), r_in=inner_radius, r_out=outer_radius)
    sky_table = aperture_photometry(image_data, sky_aperture)
    sky_flux = sky_table['aperture_sum'][0]
    sky_area = np.pi * (outer_radius ** 2 - inner_radius ** 2)
    sky_background_flux = sky_flux / sky_area

    # Subtract the sky background flux from the star flux to get the total flux
    best_flux = flux - (abs(sky_background_flux) * n_pixels)
    # Update the best aperture if the current SNR is higher
    return best_flux

def calculate_error(star_flux, gain, readnoise, n_pixels):
    noise = np.sqrt(star_flux * gain + n_pixels * readnoise ** 2)
    snr = star_flux * gain / noise
    error = 1.08 / snr
    return error, snr

# finding chosen star in each image and saving it as a cutout 
centers = []
stars = []
for i in range(0,len(images):
        image = images[i]
        star = image[star_y[i]-30 : star_y[i]+30 , star_x[i]-30 : star_x[i]+30]
        stars.append(star)
        center = ndimage.measurements.center_of_mass(star)
        centers.append(center)
        #plt.figure()
        #plt.title(f'{i+1}')
        #plt.imshow(stars[i])
        #plt.colorbar
        #plt.show()

apertures = []
get the best radii from each image:

for i in range(0,len(images)):
    center = centers[i]
    best_flux , best_aperture = find_optimal_flux(stars[i], center[1], center[0] )
    best_flux.append(best_flux)
    apertures.append(best_aperture)

plot snr with respect to r:
center = centers[13]
r = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16]
flux1 , best_aperture , aperture_radii, snr = find_optimal_flux(stars[13], center[1], center[0])
  
plt.plot(r, snr, linestyle='-', color='blue')
plt.xlabel('Radius')
plt.ylabel('Signal-to-Noise Ratio (SNR)')
plt.title('SNR with respect to radius')
#plt.show()
#for this specific star, max snr = 639.62 which has a radius of 12 pixels

# radii from each best aperture
radii = [14, 10.5, 12.5, 14.5, 12.5 , 15.5, 12.5, 11.5 , 11, 12.5 , 12.5, 14, 13, 12, 12, 13.5 , 13.5, 14, 12.5, 14]
radius = np.median(radii)
radius = 12.5

#using said radius, calculate flux and magnitude for each image 
star_flux = []
magnitudes = []
errors = []
SNRs = []
for i in range(0,len(images)):
    center = centers[i]
    best_flux = flux(stars[i], center[1], center[0] )
    star_flux.append(best_flux)
    star_mag = apparent_magnitude(best_flux)
    magnitudes.append(float(star_mag))
    error , snr = calculate_error(best_flux,1,1,np.pi * 12.5*12.5)
    SNRs.append(snr)
    errors.append(errors)
    
# get star altitudes for each image from starrynight; for example:
z = [58.47883333, 57.82233333, 57.16266667, 56.50016667, 55.66666667, 54.865, 52.43833333, 51.76, 51.07966667, 50.3965, 49.71416667,49.03333333, 48.33, 47.6433333]
sec_z = []
for i in range(0,len(images)):
     rad = z[i]*(np.pi / 180)
     secz = 1/np.cos(rad)
     sec_z.append(secz)

# differential magnitudes
differential_mag = magnitudes.copy()
for i in range (0,14):
     if i == 0:
          differential_mag[i] = magnitudes[0] - minimum_mag
     else:
          differential_mag[i] = magnitudes[i] -minimum_mag

# error bar for each point is error = 1.0875/snr.
# errors are used as weights to get a weighted line fit.
y_error = []
for i in range(0,len(images)):
     err = 1.0875/SNRs[i]
     y_error.append(err)

def func(x, a, b):
    return (a * x) + b

x = np.asarray(sec_z)

# plotting differential magnitudes in respect to secz

# gives a and b values and their respective errors:
popt, pcov = curve_fit(func, x, differential_mag, sigma=y_error, absolute_sigma=True)
yfit = func(x, *popt)
print('Weighted fit parameters:', popt)
cov_error = np.sqrt(np.diag(pcov))
print('Covariant errors:',cov_error)
print('Covariance matrix:'); print(pcov)
plt.xlabel("sec z")
plt.ylabel("differential magnitude")
plt.errorbar(x, differential_mag, yerr = y_error ,fmt='none',ecolor = 'blue',color='yellow') 
# plotting a weighted fitted line:
plt.plot(x, yfit, label='Weighted fit (WLS)')
# plotting the data:
plt.plot(x, differential_mag, '.')
plt.legend(loc='lower center')
plt.show()

#data table:
d = {'z': np.array(z), 'sec_z': np.array(sec_z), 'instrumental magnitude': np.array(magnitudes), 'differential mag' : np.array(differential_mag), 'magnitude error': y_error,'SNR' : np.array(SNRs)}
df = pd.DataFrame(data=d)
#print(df)

# calculate regression:
x = np.array(sec_z).reshape((-1, 1))
y = np.array(differential_mag)

model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print(f"coefficient of determination (Regression): {r_sq}")
print(f"intercept: {model.intercept_}")
print(f"slope: {model.coef_}")

