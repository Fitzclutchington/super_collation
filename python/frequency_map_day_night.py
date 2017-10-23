"""
Compute the frequency of clear-sky samples for each pixel within a list
of L3U files.
"""
import glob
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sio
import netCDF4

import utils

class FrequencyMap:
    

    def __init__(self, l3u_granule_list, nc_variable_name, sst):
        self.l3u_granule_list = l3u_granule_list
        self.num_files = len(l3u_granule_list)
        self.nc_variable_name = nc_variable_name
        self.sst = sst
        self.satellite_name = l3u_granule_list[0].split('/')[-1].split('-')[4]
        self.create_map_day_night('l2p_flags',2)

        

    def create_map_day_night(self, l2p_flag_variable_name, land_bit):
        filename = self.l3u_granule_list[0]
        cdf = netCDF4.Dataset(filename)
        self.height = cdf.dimensions['lat'].size
        self.width = cdf.dimensions['lon'].size

        self.frequency_map_day = np.zeros((self.height, self.width))
        self.frequency_map_night = np.zeros((self.height, self.width))
        self.frequency_map_day_masked = np.zeros((self.height, self.width))
        self.frequency_map_night_masked = np.zeros((self.height, self.width))
        self.frequency_map_day_quality = np.zeros((self.height, self.width))
        self.frequency_map_night_quality = np.zeros((self.height, self.width))
        if self.sst:
            self.sst_map_day = np.full((self.height, self.width), np.nan)
            self.sst_map_night = np.full((self.height, self.width), np.nan)
    
            
        self.land_layer = np.zeros((self.height, self.width, 4))
        

    def compute_overlay_frequency_day_night(self):

        r = 146/256.0
        g = 98/256.0
        b = 57/256.0

        for i,filename in enumerate(self.l3u_granule_list):
            ql_mat = utils.read_var(filename, self.nc_variable_name)
            l2p_mat = utils.read_var(filename, 'l2p_flags')

            land_mask = np.bitwise_and(l2p_mat, 2).astype(bool)
            day = np.bitwise_and(l2p_mat,512).astype(bool)

            quality_mask = ql_mat == 5
            valid_mask = ql_mat >= 0
            day_mask_valid = (valid_mask & day)
            night_mask_valid = (valid_mask & (~day))
            day_mask_quality = (quality_mask & day)
            night_mask_quality = (quality_mask & (~day))

            self.land_layer[land_mask] = [r,g,b,1]
            self.frequency_map_day[day_mask_valid] += 1
            self.frequency_map_night[night_mask_valid] += 1
            
            self.frequency_map_day_masked[day_mask_valid] += 1
            self.frequency_map_night_masked[night_mask_valid] += 1
            
            self.frequency_map_day_quality[day_mask_quality] += 1
            self.frequency_map_night_quality[night_mask_quality] += 1
            
            if self.sst:
                sst_mat = utils.read_var(filename, 'sea_surface_temperature')
                self.sst_map_day[day_mask] = sst_mat[day_mask]
                self.sst_map_night[night_mask] = sst_mat[night_mask]
    
           
            print i, '/', self.num_files

        day_mask = self.frequency_map_day_quality == 0
        night_mask = self.frequency_map_night_quality == 0
        self.frequency_map_day_masked[day_mask] = np.nan
        self.frequency_map_night_masked[night_mask] = np.nan
        self.frequency_map_day_quality[day_mask] = np.nan
        self.frequency_map_night_quality[night_mask] = np.nan

    def display_map_day_night(self):

        fig, ax1 = plt.subplots(1)        
        img1 = ax1.imshow(self.frequency_map_day, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        ax1.set_title(self.satellite_name + " Day")
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)
     
        img1 = ax1.imshow(self.frequency_map_night, cmap='jet')
        ax1.set_title(self.satellite_name + " Night")
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)        
        img1 = ax1.imshow(self.frequency_map_day_masked, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        ax1.set_title(self.satellite_name + " Day")
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)
     
        img1 = ax1.imshow(self.frequency_map_night_masked, cmap='jet')
        ax1.set_title(self.satellite_name + " Night")
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)        
        img1 = ax1.imshow(self.frequency_map_day_quality, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        ax1.set_title(self.satellite_name + " Day")
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)
     
        img1 = ax1.imshow(self.frequency_map_night_quality, cmap='jet')
        ax1.set_title(self.satellite_name + " Night")
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        if self.sst:
            fig, ax1 = plt.subplots(1)
            img1 = ax1.imshow(self.sst_map_day, vmin=270, vmax=307,  cmap='jet')       
            ax1.imshow(self.land_layer, interpolation='nearest')
            ax1.set_title("Day")
            div1 = make_axes_locatable(ax1)
            cax1 = div1.append_axes("right", size="5%", pad=0.05)
            cbar1 = plt.colorbar(img1, cax=cax1)
            plt.show()

            fig, ax1 = plt.subplots(1)
        
            img1 = ax1.imshow(self.sst_map_night, vmin=270, vmax=307, cmap='jet')
            ax1.set_title("Night")
            ax1.imshow(self.land_layer, interpolation='nearest')
            div1 = make_axes_locatable(ax1)
            cax1 = div1.append_axes("right", size="5%", pad=0.05)
            cbar1 = plt.colorbar(img1, cax=cax1)
            plt.show()


    def save_map_day_night(self, save_loc):
        data = { "frequencies_day": self.frequency_map_day,  "frequencies_night": self.frequency_map_night }
        sio.savemat(save_loc, data)


save_loc = sys.argv[2]
l3u_folder = sys.argv[1]
l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))

l3u_variable = 'sea_surface_temperature'

vmin = 271.15
vmax = 310

"""
frequency_map = FrequencyMap(l3u_file_list, 'sea_surface_temperature')
frequency_map.compute_frequencies(vmin,vmax)
"""

frequency_map = FrequencyMap(l3u_file_list, 'quality_level',False)
frequency_map.compute_overlay_frequency_day_night()
frequency_map.display_map_day_night()
#frequency_map.save_map_day_night(save_loc)
"""
frequency_map.display_map()
frequency_map.save_map(save_loc)
"""

"""
day_mask = np.isnan(frequency_map.sst_map_day)
night_mask = np.isnan(frequency_map.sst_map_night)
frequency_map.frequency_map_day[day_mask] = np.nan
frequency_map.frequency_map_night[night_mask] = np.nan
"""