import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean
import sys

import utils

filepath = sys.argv[1]
filename = filepath.split('/')[-1]
time = filename[8:10] + ':' + filename[10:12]

#crop = np.s_[1000:1250,4850:5150]
crop_1 = np.s_[1000:1250,4850:5150]
crop_2 = np.s_[2250:2750,5500:6000]
crops = [crop_1, crop_2]

vmins = [273, 287]
vmaxs = [278, 303]

sst = utils.read_var(filepath, 'sea_surface_temperature')
sses_bias = utils.read_var(filepath, 'sses_bias')
l2p_flags = utils.read_var(filepath,'l2p_flags')
lats = utils.read_var(filepath, 'lat')
lons = utils.read_var(filepath, 'lon')

land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((land_mask.shape[0],land_mask.shape[1],4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

sst_corrected = sst - sses_bias



step = 60
angle = -45
tick_font = { 'size': 15  }

for crop, vmin, vmax in zip(crops, vmins, vmaxs):
	x_min = crop[1].start 
	x_max = crop[1].stop 
	y_min = crop[0].start 
	y_max = crop[0].stop 

	ylen = y_max - y_min 
	xlen = x_max - x_min 


	plt.figure()

	plt.suptitle(time, fontsize=30)

	ax1 = plt.subplot(121)
	ax1.set_title("SST")
	ax1.set_xticks(np.arange(0,xlen,step))
	ax1.set_yticks(np.arange(0,ylen,step))
	ax1.set_xticklabels(['{0:.2f}'.format(x) for x in lons[x_min:x_max+1:step]],fontdict=tick_font, rotation=angle)
	ax1.set_yticklabels(['{0:.2f}'.format(x) for x in lats[y_min:y_max+1:step]],fontdict=tick_font)
	img1 = ax1.imshow(sst[crop],vmin=vmin,vmax=vmax, cmap='jet')
	ax1.imshow(land[crop],interpolation='nearest')
	div1 = make_axes_locatable(ax1)
	cax1 = div1.append_axes("right", size="5%", pad=0.05)
	cax1.set_visible(False)
	#cbar1 = plt.colorbar(img1, cax=cax1)
	#cbar1.ax.tick_params(labelsize=15)

	ax2 = plt.subplot(122, sharex=ax1,sharey=ax1)
	ax2.set_title("SST - SSES BIAS")
	ax2.set_xticks(np.arange(0,xlen,step))
	ax2.set_yticks(np.arange(0,ylen,step))
	ax2.set_xticklabels(['{0:.2f}'.format(x) for x in lons[x_min:x_max:step]],fontdict=tick_font, rotation=angle)
	ax2.set_yticklabels(['{0:.2f}'.format(x) for x in lats[y_min:y_max:step]],fontdict=tick_font)
	img2 = ax2.imshow(sst_corrected[crop], vmin=vmin, vmax=vmax, cmap='jet')
	ax2.imshow(land[crop],interpolation='nearest')
	div2 = make_axes_locatable(ax2)
	cax2 = div2.append_axes("right", size="5%", pad=0.05)
	cbar2 = plt.colorbar(img2, cax=cax2)
	cbar2.ax.tick_params(labelsize=15)


	plt.show()
