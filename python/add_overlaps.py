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
    

    def __init__(self, data_file_1, data_file_2,land_file):
        self.data_file_1 = data_file_1
        self.data_file_2 = data_file_2
        self.variable_list = ["full_frequency_day", "full_frequency_night",
                              "masked_frequency_day", "masked_frequency_night",
                              "quality_frequency_day", "quality_frequency_night" ]
        self.satellite_name_1 = data_file_1.split('_')[-1]
        self.satellite_name_2 = data_file_2.split('_')[-1]
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

        data_1 = sio.loadmat(self.data_file_1)
        data_2 = sio.loadmat(self.data_file_2)

        self.data_list_1 = []
        self.data_list_2 = []
        for i,key in enumerate(self.variable_list):
            self.data_list_1.append(data_1[key])
            self.data_list_2.append(data_2[key])

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
            ax1.set_title(self.satellite_name_1 + " " + self.satellite_name_2 + " Day")
        else:
            ax1.set_title(self.satellite_name_1 + " " + self.satellite_name_2 + " Night")

        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)#,ticks=np.linspace(0,max_val,max_val+1,dtype=np.int8))
        plt.show()

    def display_map_day_night(self):       
    

        for i,(variable,data_1,data_2) in enumerate(zip(self.variable_list,self.data_list_1, self.data_list_2)):
            if data_1.dtype != np.uint8:
                data_1[np.isnan(data_1)] = 0
                data_2[np.isnan(data_2)] = 0
            full_mat = data_1 + data_2
            if data_1.dtype != np.uint8:
                full_mat[full_mat==0] = np.nan
            day_night = i%2
            cmap = plt.cm.jet
            #max_val = 16
            self.plot_mat(full_mat, cmap, day_night)



    



data_2_file = sys.argv[2]
data_1_file = sys.argv[1]
land_file = sys.argv[3]


frequency_map = FrequencyMap(data_1_file, data_2_file, land_file)
frequency_map.compute_overlay_frequency_day_night()
frequency_map.display_map_day_night()
