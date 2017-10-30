"""
Compute the frequency of clear-sky samples for each pixel within a list
of L3U files.
"""
import glob
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sio
import netCDF4

import utils

class FrequencyMap:
    

    def __init__(self, data_files, land_file):
        self.data_file_1 = data_file_1
        self.data_file_2 = data_file_2        
        self.variable = "quality_frequency_night"
        self.land_file = land_file
        self.create_map_day_night('l2p_flags',2)

    def create_map_day_night(self, l2p_flag_variable_name, land_bit):
        filename = "/home/fitz/viirs/20171007235000-STAR-L3U_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.50B04-v02.0-fv01.0.nc"
        cdf = netCDF4.Dataset(filename)
        self.height = cdf.dimensions['lat'].size
        self.width = cdf.dimensions['lon'].size

    
        self.lats = np.squeeze(cdf['lat'][:])
        self.lons = np.squeeze(cdf['lon'][:])
        land_mask = sio.loadmat(self.land_file)['land_mask']
        self.land_layer = np.zeros((self.height, self.width, 4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0
        self.land_layer[land_mask==255] = [r,g,b,1]


    def compute_overlay_frequency_day_night(self):

        self.data_1 = sio.loadmat(self.data_file_1)[variable]
        self.data_2 = sio.loadmat(self.data_file_2)[variable]


    def plot_mat(self, mat, cmap, day_night):
        step = 1000
        angle = -45
        tick_font = { 'size': 15  }

        fig, ax1 = plt.subplots(1)

        #mat_min = np.nanmin(mat)
        #mat_max = np.nanmax(mat)
        #norm = mpc.BoundaryNorm(np.arange(-.5 , max_val + 1 ,1), cmap.N)   

        img1 = ax1.imshow(mat, cmap=cmap)#, norm=norm)
        ax1.imshow(self.land_layer, interpolation='nearest')

        ax1.set_xticks(np.arange(0,self.width,step))
        ax1.set_yticks(np.arange(0,self.height,step))
        ax1.set_xticklabels(['{0:.2f}'.format(x) for x in self.lons[::step]],fontdict=tick_font, rotation=angle)
        ax1.set_yticklabels(['{0:.2f}'.format(x) for x in self.lats[::step]],fontdict=tick_font)

        if day_night == 0:
            ax1.set_title("AQUA VIIRS Day")
        else:
            ax1.set_title("AQUA VIIRS Night")

        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)#,ticks=np.linspace(0,max_val,max_val+1,dtype=np.int8))
        plt.show()

    def display_map_day_night(self):       
    

        for i, variable in enumerate(self.variable_list):
            full_mat = self.data[variable]
            day_night = i%2
            cmap = plt.cm.jet
            #max_val = 16
            self.plot_mat(full_mat, cmap, day_night)



crop_file = sys.argv[1]

all_sensors_data = sio.loadmat('../data/overlap_data/frequency_all_sensors.mat')['quality_frequency_night']
viirs_data =  sio.loadmat('../data/overlap_data/frequency_viirs.mat')['quality_frequency_night']

filename = "/home/fitz/viirs/20171007235000-STAR-L3U_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.50B04-v02.0-fv01.0.nc"
cdf = netCDF4.Dataset(filename)
height = cdf.dimensions['lat'].size
width = cdf.dimensions['lon'].size

crops = utils.get_crops(crop_file)

lats = np.squeeze(cdf['lat'][:])
lons = np.squeeze(cdf['lon'][:])

land_mask = sio.loadmat("../data/land_mask.mat")['land_mask']
land_layer = np.zeros((height, width, 4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land_layer[land_mask==255] = [r,g,b,1]


angle = -45
tick_font = { 'size': 15  }

for crop in crops:
    x_min = crop[1].start
    x_max = crop[1].stop
    y_min = crop[0].start
    y_max = crop[0].stop

    xlen = x_max - x_min
    ylen = y_max - y_min

    step_x = xlen/10
    step_y = ylen/10

    fig, (ax1, ax2) = plt.subplots(1,2,sharex=True,sharey=True)
    cmap = plt.cm.jet
    norm = mpc.BoundaryNorm(np.arange(-.5 , 7 + 1 ,1), cmap.N)   
    img1 = ax1.imshow(viirs_data[crop], cmap=cmap, norm=norm)
    ax1.imshow(land_layer[crop], interpolation='nearest')

    ax1.set_xticks(np.arange(0,xlen,step_x))
    ax1.set_yticks(np.arange(0,ylen,step_y))
    ax1.set_xticklabels(['{0:.2f}'.format(x) for x in lons[x_min+1:x_max+1:step_x]],fontdict=tick_font, rotation=angle)
    ax1.set_yticklabels(['{0:.2f}'.format(x) for x in lats[y_min+1:y_max+1:step_y]],fontdict=tick_font)

    ax1.set_title("VIIRS")

    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1,ticks=np.linspace(0,7,7+1,dtype=np.int8))

    cmap = plt.cm.Spectral_r
    img1 = ax2.imshow(all_sensors_data[crop], cmap=cmap, vmin=1, vmax=14)
    ax2.imshow(land_layer[crop], interpolation='nearest')

    ax2.set_xticks(np.arange(0,xlen,step_x))
    ax2.set_yticks(np.arange(0,ylen,step_y))
    ax2.set_xticklabels(['{0:.2f}'.format(x) for x in lons[x_min+1:x_max+1:step_x]],fontdict=tick_font, rotation=angle)
    ax2.set_yticklabels(['{0:.2f}'.format(x) for x in lats[y_min+1:y_max+1:step_y]],fontdict=tick_font)
    ax2.set_title("All Sensors")

    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    plt.show()