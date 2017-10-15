import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean
import sys

import utils

filepath_1 = sys.argv[1]
filepath_2 = sys.argv[2]
row = int(sys.argv[3])
col = int(sys.argv[4])

filename_1 = filepath_1.split('/')[-1]
filename_2 = filepath_2.split('/')[-1]

crop_size = 400
lag = crop_size/2

crop = np.s_[row-lag:row+lag , col-lag:col+lag]

sst_1 = utils.read_var(filepath_1, 'sea_surface_temperature')[crop]
sza_1 = utils.read_var(filepath_1, 'satellite_zenith_angle')[crop]
dtime_1 = utils.read_var(filepath_1, 'sst_dtime')[crop]
dtimes_1 = float(filename_1[8:10])  + float(filename_1[10:12])/60 + (float(filename_1[12:14]) + dtime_1)/(60*60)
#sses_1 = utils.read_var(filepath_1, 'sses_standard_deviation')[crop]
l2p_flags_1 = utils.read_var(filepath_1,'l2p_flags')[crop]
land_mask = np.bitwise_and(l2p_flags_1,2).astype(bool)[crop]
land = np.zeros((crop_size,crop_size,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

sst_2 = utils.read_var(filepath_2, 'sea_surface_temperature')[crop]
sza_2 = utils.read_var(filepath_2, 'satellite_zenith_angle')[crop]
#sses_2 = utils.read_var(filepath_2, 'sses_standard_deviation')[crop]
l2p_flags_2 = utils.read_var(filepath_2,'l2p_flags')[crop]
dtime_2 = utils.read_var(filepath_2, 'sst_dtime')[crop]
dtimes_2 = float(filename_2[8:10])  + float(filename_2[10:12])/60 + (float(filename_2[12:14]) + dtime_2)/(60*60)

#Collation
mask_1 = l2p_flags_1 == 0
mask_2 = l2p_flags_2 == 0
mask_3 = ~mask_1 & ~mask_2 & ~sst_1.mask & ~sst_2.mask

sst_collated = np.full(mask_1.shape, np.nan)
time_collated = np.full(mask_1.shape, np.nan)

w1 = np.exp(-((sza_1/ 30.0)**2))
w2 = np.exp(-((sza_2/ 30.0)**2))

sst_collated[mask_3] = ((w1*sst_1 + w2*sst_2)/ (w1+w2))[mask_3]

sst_collated[mask_1] = sst_2[mask_1]
sst_collated[mask_2] = sst_1[mask_2]

time_collated[mask_3] = ((w1*dtimes_1 + w2*dtimes_2)/ (w1+w2))[mask_3]

time_collated[mask_1] = dtimes_2[mask_1]
time_collated[mask_2] = dtimes_1[mask_2]



plt.figure()
#cmap = cmocean.cm.thermal
vmin = 282
vmax = 305
ax1 = plt.subplot(231)
img1 = ax1.imshow(sst_1,vmin=vmin,vmax=vmax, cmap='jet')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(232, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_collated-sst_1, vmin=-0.5, vmax=0.5, cmap='jet')
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(233, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_collated, vmin=vmin, vmax=vmax, cmap='jet')
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(234, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_2,vmin=vmin,vmax=vmax, cmap='jet')
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(235, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_collated-sst_2, vmin=-0.5, vmax=0.5, cmap='jet')
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(236, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_1-sst_2, vmin=-.5, vmax=.5, cmap='jet')
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

plt.show()

plt.figure()
plt.plot(sst_1[200,:],label='sst 1')
plt.plot(sst_2[200,:],label='sst 2')
plt.plot(sst_collated[200,:], label='sst collated')
plt.legend()
plt.show()

plt.figure()
plt.plot(sst_collated[200,:] - sst_2[200,:], label='sst 2')
plt.plot(sst_collated[200,:] - sst_1[200,:], label='sst 1')
plt.legend()
plt.show()