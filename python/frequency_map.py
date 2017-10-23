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

import utils

class FrequencyMap:

    

    def __init__(self, l3u_granule_list, nc_variable_name, sat_type, day_night):
        self.l3u_granule_list = l3u_granule_list
        self.num_files = len(l3u_granule_list)
        self.nc_variable_name = nc_variable_name
        self.sat_type = sat_type
        if day_night:
            self.create_map_day_night('l2p_flags',2)
        else:
            self.create_map('l2p_flags', 2)
        

    def create_map(self, l2p_flag_variable_name, land_bit):
        filename = self.l3u_granule_list[0]
        land_mask_val = (1 << (land_bit - 1))
        land_mask = np.bitwise_and(utils.read_var(filename, l2p_flag_variable_name), land_mask_val).astype(bool)
        self.frequency_map = np.zeros(land_mask.shape)

        self.land_layer = np.zeros((land_mask.shape[0], land_mask.shape[1], 4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0

        if self.sat_type == 'polar':
            for filename in self.l3u_granule_list:
                mat = utils.read_var(filename, l2p_flag_variable_name)
                land_mask = np.bitwise_and(mat, land_mask_val).astype(bool)
                ice_mask = np.bitwise_and(mat, 4).astype(bool)
                self.land_layer[land_mask] = [r,g,b,1]
                self.frequency_map[ice_mask] = np.nan


    def create_map_day_night(self, l2p_flag_variable_name, land_bit):
        filename = self.l3u_granule_list[0]
        land_mask_val = (1 << (land_bit - 1))
        land_mask = np.bitwise_and(utils.read_var(filename, l2p_flag_variable_name), land_mask_val).astype(bool)
        self.frequency_map_day = np.zeros(land_mask.shape)
        self.frequency_map_night = np.zeros(land_mask.shape)

        self.land_layer = np.zeros((land_mask.shape[0], land_mask.shape[1], 4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0

        if self.sat_type == 'polar':
            for i,filename in enumerate(self.l3u_granule_list):
                mat = utils.read_var(filename, l2p_flag_variable_name)
                land_mask = np.bitwise_and(mat, land_mask_val).astype(bool)
                ice_mask = np.bitwise_and(mat, 4).astype(bool)
                self.land_layer[land_mask] = [r,g,b,1]
                self.frequency_map_day[ice_mask] = np.nan
                self.frequency_map_night[ice_mask] = np.nan
                print i, '/', self.num_files

    def compute_frequencies(self, vmin, vmax):

        for i,filename in enumerate(self.l3u_granule_list):
            mat = utils.read_var(filename, self.nc_variable_name)
            mask = (mat >= vmin) & (mat <= vmax)
            self.frequency_map[mask] += 1
            print i, '/', self.num_files

    def compute_overlay_frequency(self):

         for i,filename in enumerate(self.l3u_granule_list):
            mat = utils.read_var(filename, self.nc_variable_name)
            mask = (mat >= 0)
            self.frequency_map[mask] += 1
            print i, '/', self.num_files

    def compute_overlay_frequency_day_night(self):

         for i,filename in enumerate(self.l3u_granule_list):
            mat = utils.read_var(filename, self.nc_variable_name)
            l2p_mat = utils.read_var(filename, 'l2p_flags')
            day = np.bitwise_and(l2p_mat,512).astype(bool)
            day_mask = ((mat >= 0) & day)
            night_mask = ((mat >= 0) & (~day))
            self.frequency_map_day[day_mask] += 1
            self.frequency_map_night[night_mask] += 1
            print i, '/', self.num_files


    def display_map(self):

        fig, ax1 = plt.subplots(1)

        img1 = ax1.imshow(self.frequency_map, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

    def display_map_day_night(self):

        fig, ax1 = plt.subplots(1)
        img1 = ax1.imshow(self.frequency_map_day, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        ax1.set_title("Day")
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

        fig, ax1 = plt.subplots(1)
        img1 = ax1.imshow(self.frequency_map_night, cmap='jet')
        ax1.set_title("Night")
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

    def save_map(self, save_loc):
        data = { "frequencies": self.frequency_map  }
        sio.savemat(save_loc, data)

    def save_map_day_night(self, save_loc):
        data = { "frequencies_day": self.frequency_map_day,  "frequencies_night": self.frequency_map_night }
        sio.savemat(save_loc, data)

def main():

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
    frequency_map = FrequencyMap(l3u_file_list, 'quality_level', 'polar', True)
    frequency_map.compute_overlay_frequency_day_night()
    frequency_map.display_map_day_night()
    """
    frequency_map.display_map()
    frequency_map.save_map(save_loc)
    """
if __name__ == '__main__':
    main()