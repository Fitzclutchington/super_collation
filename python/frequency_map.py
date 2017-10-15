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

import utils

class FrequencyMap:

    

    def __init__(self, l3u_granule_list, nc_variable_name):
        self.l3u_granule_list = l3u_granule_list
        self.nc_variable_name = nc_variable_name
        self.create_map('l2p_flags', 2)

    def create_map(self, l2p_flag_variable_name, land_bit):
        filename = self.l3u_granule_list[0]
        land_mask_val = (1 << (land_bit - 1))
        land_mask = np.bitwise_and(utils.read_var(filename, l2p_flag_variable_name), land_mask_val).astype(bool)
        self.frequency_map = np.zeros(land_mask.shape, dtype=np.int16)

        self.land_layer = np.zeros((land_mask.shape[0], land_mask.shape[1], 4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0
        self.land_layer[land_mask] = [r,g,b,1]


    def compute_frequencies(self, vmin, vmax):

        num_files = len(self.l3u_granule_list)

        for i,filename in enumerate(self.l3u_granule_list):
            mat = utils.read_var(filename, self.nc_variable_name)
            mask = (mat >= vmin) & (mat <= vmax)
            self.frequency_map[mask] += 1
            print i, '/', num_files

    def display_map(self):

        fig, ax1 = plt.subplots(1)

        img1 = ax1.imshow(self.frequency_map, cmap='jet')
        ax1.imshow(self.land_layer, interpolation='nearest')
        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1)
        plt.show()

def main():

    l3u_folder = sys.argv[1]
    l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))

    l3u_variable = 'sea_surface_temperature'

    vmin = 271.15
    vmax = 310

    frequency_map = FrequencyMap(l3u_file_list, 'sea_surface_temperature')
    frequency_map.compute_frequencies(vmin,vmax)
    frequency_map.display_map()

if __name__ == '__main__':
    main()