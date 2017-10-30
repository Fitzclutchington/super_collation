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

class LandMap:
    
    def __init__(self, l3u_granule_list):
        self.l3u_granule_list = l3u_granule_list
        self.num_files = len(l3u_granule_list)
        self.create_map()

        

    def create_map(self):
        filename = self.l3u_granule_list[0]
        cdf = netCDF4.Dataset(filename)
        self.height = cdf.dimensions['lat'].size
        self.width = cdf.dimensions['lon'].size
         
        self.land_layer = np.zeros((self.height, self.width), dtype=np.uint8)
        

    def compute_land(self):

        for i,filename in enumerate(self.l3u_granule_list):
            l2p_mat = utils.read_var(filename, 'l2p_flags')
            land_mask = np.bitwise_and(l2p_mat, 2).astype(bool)
            self.land_layer[land_mask] = 255           
            print i, '/', self.num_files


    def display_land(self):

        fig, ax1 = plt.subplots(1)        
        ax1.imshow(self.land_layer, interpolation='nearest')
        ax1.set_title("Land Cover")
        plt.show()



    def save_land(self, save_loc):
        data = { "land_mask": self.land_layer }
        sio.savemat(save_loc, data)


save_loc = sys.argv[2]
l3u_folder = sys.argv[1]
l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))

l3u_variable = 'sea_surface_temperature'

land_map = LandMap(l3u_file_list)
land_map.compute_land()
land_map.display_land()
land_map.save_land(save_loc)
