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


    

def create_index_stack(l3u_file_list, day_flag=True):

    filename = l3u_file_list[0]
    height,width = utils.get_dimensions(filename)
    depth = 8

    file_ind_stack = np.full((height, width, depth),255, dtype=np.uint8)
    inds = np.zeros((height, width), dtype=np.uint8)
    quality_data = np.zeros((height, width), dtype=np.uint8)

    frequency = np.zeros((height, width), dtype=np.uint8)
    num_files = len(l3u_file_list)

    for i,filename in enumerate(l3u_file_list):

        l2p_mat = utils.read_var(filename, 'l2p_flags')
        if day_flag:
            time = np.bitwise_and(l2p_mat,512).astype(bool)
        else:
            time = ~np.bitwise_and(l2p_mat,512).astype(bool)
        if time.sum() == 0:
            continue
        else:
            ql_mat = utils.read_var(filename, 'quality_level')
            ice_mask = ~np.bitwise_and(l2p_mat,4).astype(bool)
            land_mask = ~np.bitwise_and(l2p_mat, 2).astype(bool)

            high_quality_mask = ql_mat == 5
            valid_mask = ql_mat >= 0
            full_mask = time & ice_mask & land_mask
            time_quality_mask = (high_quality_mask & full_mask)
            time_valid_mask = (valid_mask & full_mask)
            
            depth_slices = inds[time_valid_mask]
            file_ind_stack[time_valid_mask,depth_slices] = i
            inds[time_valid_mask] += 1
            quality_data[time_quality_mask] += 1
            
            print i+1, '/', num_files

    low_quality_mask = quality_data==0
    file_ind_stack[low_quality_mask,:] = 255
    return file_ind_stack

l3u_folder = sys.argv[1]
l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))


vmin = 271.15
vmax = 310

ind_stack = create_index_stack(l3u_file_list)

save_dict = {}

for i in range(7):
    name = "slice_"+str(i)
    save_dict[name] = ind_stack[:,:,i]
sio.savemat("index_stack",save_dict)