# code was written by Aysan Hemmati

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy import units as u
from scipy import ndimage
from photutils.profiles import RadialProfile
from scipy.optimize import curve_fit

stars_fits = fits.open(r"stars .fits")
data = stars_fits[0].data
#info = stars_fits.info()
# show the fits file:
#plt.figure()
#plt.imshow(data, cmap='gray')
#plt.colorbar()
#plt.show() #Figure 1'
coordiantes = ((1593,1225),(1680,1263),(1238,1295),(1223,1495),(420,1370),(1350,1608),(1703,1783),(1770,1823),(1795,1903),(1923,1420),(62,1440),(732,205),(432,487),(1655,667),(1966,455),(465,602),(290,710),(1610,115),(1563,345),(1610,112))

stars = []
centers = []
for i in range (0,20):
        s = data[coordiantes[i][1]-30 : coordiantes[i][1]+30 , coordiantes[i][0]-30 : coordiantes[i][0]+30]
        max = np.max(s)
        norm = np.array(s/max)
        stars.append(norm)
        center = ndimage.measurements.center_of_mass(s)
        centers.append(center)
        
centers_tuple = tuple(centers) 

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 
plots = []
for i in range (0,20):
        s = data[coordiantes[i][1]-30 : coordiantes[i][1]+30 , coordiantes[i][0]-30 : coordiantes[i][0]+30]
        max = np.max(s)
        norm = np.array(s/max)
        stars.append(norm)
        center = ndimage.measurements.center_of_mass(s)
        x = radial_profile(s, center)
        plots.append(x)
        
y = plots[8]
x = range(0,len(y))
mean = sum(x * y) / sum(y)
        
sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

def Gauss(x, a,x0, sigma , b):
        return b +(a * np.exp(-(x-x0)**2 / (2 * sigma**2)))
                  
popt,pcov = curve_fit(Gauss, x, y)
print('fit parameters:', popt) 
plt.plot(x, y , 'b+:', label='data')
plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
plt.legend()
plt.title('Radial profile of star 20')
plt.xlabel('Radius (pixels)')
plt.ylabel('Brightness')
plt.show()