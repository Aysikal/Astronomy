import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import reza
from scipy.interpolate import griddata
import time
import scipy.optimize as optimize
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable


t1=time.time()
light_file = fits.open(r"C:/Users\AYSAN\Desktop/project\data\BW/R/1/light-r-2023_10_10-exp00.02.00.000-1x1_High_3.fit")
light = light_file[0].data

dark_pixel_locations = []

higher_limit = 150
lower_limit = 0
# Iterate over each row and column
for i in range(len(light)):
    for j in range(len(light[i])):
        if lower_limit < light[i][j] < higher_limit:
            dark_pixel_locations.append((i, j , light[i][j]))
    
print("Dark pixel counts (values less than",higher_limit, "and larger than" , lower_limit, "):",len(dark_pixel_locations))
array = np.array(dark_pixel_locations)

random_numbers = np.random.randint(0, len(dark_pixel_locations), size=100000)

new_list = []
for i in random_numbers:
    new = array[i]
    new_list.append(new)

xs = []
ys = []
zs = []
for i in range(0,len(new_list)):
    x = new_list[i][0]
    xs.append(x)
    y = new_list[i][1]
    ys.append(y)
    z = new_list[i][2]
    zs.append(z)

#heat map of data:
X, Y = np.meshgrid(np.linspace(np.min(xs), np.max(xs), 100),
                   np.linspace(np.min(ys), np.max(ys), 100))
Z = griddata((xs, ys), zs, (X, Y), method='linear', fill_value=0)

fig, ax = plt.subplots()
art = ax.pcolor(X, Y, Z, shading='auto')
plt.colorbar(art, ax=ax, label='z')
ax.set_xlabel('x')
ax.set_ylabel('y')
#plt.show()

#fit
xy = np.vstack((xs , ys)).T


def plane_fit(X , a , b , c , d ):
    return a*X[0] + b*X[1] + c*X[0]*X[1] + d

def plane_fit2(X , a , b , c , d , e , f):
    return a*X[0] + b*X[1] + c*X[0]*X[1] + d*(x**2) + e*(y**2) +f

guess = (1,1,1,1,1,1)
popt, pcov = optimize.curve_fit(plane_fit, np.transpose(xy), zs)

zfit = plane_fit(xy, *popt)
print(popt)

a = float(popt[0])
b = float(popt[1])
c = float(popt[2])
d = float(popt[3])
#e = float(popt[4])
#f = float(popt[5])

#plot the plane: 

X, Y = np.meshgrid(np.linspace(np.min(xs), np.max(xs), 1000),
                   np.linspace(np.min(ys), np.max(ys), 1000))
def plane(X , Y , a , b , c , d ):
    return a*X + b*Y + c*X*Y+d

def plane2(X , Y , a , b , c , d , e , f):
    return a*X + b*Y + c*X*Y+d*(x**2)+e*(y**2)+f

fig = plt.figure()
plt.title(("ax+by+cxy+d"))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, plane(X , Y , a , b , c , d ))

plt.show()

#correction:

numRows = len(light)
numCols = len(light[0])

# Create a new 2D array with double the values of the original array
corrected_light = [[light[i][j] - plane(i, j , a, b, c, d) for j in range(numCols)] for i in range(numRows)]
background = [[plane(i, j , a, b, c, d) for j in range(numCols)] for i in range(numRows)]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('raw vs corrected ( ax+by+cxy+d)')

fig = plt.figure(figsize=(16, 12))

# First subplot
ax1 = fig.add_subplot(121)
im1 = ax1.imshow(light, interpolation='None')
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')
ax1.set_title("raw light")

# Second subplot
ax2 = fig.add_subplot(122)
im2 = ax2.imshow(corrected_light, interpolation='None')
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical')
ax2.set_title("dark corrected light")

plt.show()

