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
granule_name = sys.argv[2]
filename = filepath.split('/')[-1]
time = filename[8:10] + ':' + filename[10:12]

crops = []
crops.append(np.s_[2400:3600,2400:3600])
crops.append(np.s_[5100:5600,14600:15100])
crops.append(np.s_[850:1100,12050:12450])
crops.append(np.s_[950:1250,4800:5300])
crops.append(np.s_[1700:2600,800:1500])

vmins = [285,295,270,270,275]
vmaxs = [305,307,285,280,290]

sst = sio.loadmat(filepath)['frequencies_day']
sst[sst==0] = np.nan

l2p_flags = utils.read_var(granule_name,'l2p_flags')
lats = utils.read_var(granule_name, 'lat')
lons = utils.read_var(granule_name, 'lon')

land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((land_mask.shape[0],land_mask.shape[1],4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]


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

	ax1 = plt.subplot(111)
	ax1.set_title("SST")
	ax1.set_xticks(np.arange(0,xlen,step))
	ax1.set_yticks(np.arange(0,ylen,step))
	ax1.set_xticklabels(['{0:.2f}'.format(x) for x in lons[x_min:x_max+1:step]],fontdict=tick_font, rotation=angle)
	ax1.set_yticklabels(['{0:.2f}'.format(x) for x in lats[y_min:y_max+1:step]],fontdict=tick_font)
	img1 = ax1.imshow(sst[crop],vmin=vmin,vmax=vmax, cmap='jet')
	ax1.imshow(land[crop],interpolation='nearest')
	div1 = make_axes_locatable(ax1)
	cax1 = div1.append_axes("right", size="5%", pad=0.05)
	cbar1 = plt.colorbar(img1, cax=cax1)
	cbar1.ax.tick_params(labelsize=15)


	plt.show()
