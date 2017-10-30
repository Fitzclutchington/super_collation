import glob
import sys
import os
import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sio
from scipy import interpolate
import netCDF4

import utils

save_folder = sys.argv[3]
folders = { 
            'l3u' : sys.argv[1],
            'l2p' : sys.argv[2]
        }

files = {
         'l3u' : sorted(glob.glob(folders['l3u']+'*.nc')),
         'l2p' : sorted(glob.glob(folders['l2p']+"*.nc")) 
        }

sets = { 
         'l3u' : set(files['l3u']),
         'l2p' : set(files['l2p'])
        }

endings = { 
            'l3u': utils.compute_ending(files['l3u'][0]),
            'l2p': utils.compute_ending(files['l2p'][0]) 
        }

dimensions = { 
               'l3u' : (utils.get_dimensions(files['l3u'][0],'lat','lon')),
               'l2p' : (utils.get_dimensions(files['l2p'][0],'nj','ni'))
            }

current_day_string = files['l3u'][0].split('/')[-1].split('-')[0][0:8]

current_time = datetime.datetime.strptime(current_day_string+'0000', '%Y%m%d%H%M%S')
time_delta = datetime.timedelta(minutes=10)

next_day = current_time + datetime.timedelta(days=1)

l3u_sst = np.full(dimensions['l3u'],np.nan)
l3u_sza = np.full(dimensions['l3u'],np.nan)
l3u_time = np.full(dimensions['l3u'],np.nan)

sza_interp = np.zeros(dimensions['l2p'])

#while current_time != next_day:
for i in range(2):
    date_time_string = current_time.strftime('%Y%m%d%H%M%S')
    l3u_filename = folders['l3u'] + date_time_string + endings['l3u']
    l2p_filename = folders['l2p'] + date_time_string + endings['l2p']

    # append sst
    if l3u_filename in sets['l3u']:
        current_sst = utils.read_var(l3u_filename,'sea_surface_temperature')
        QL_flags = utils.read_var(l3u_filename, 'quality_level')
        QL_mask = QL_flags >= 0
        l3u_sst[QL_mask] = current_sst[QL_mask]
        
        # compute time and sza
        if l2p_filename in sets['l2p']:
            sza = utils.read_var(l2p_filename,'satellite_zenith_angle')
            sza_diff = np.diff(sza)
            end = sza.shape[1]
            num_rows = sza.shape[0]
            x_inds = np.arange(end)
            for i in range(num_rows):
                edges = np.where(sza_diff[i,:] != 0)[0]
                mid_ind = np.zeros(edges.shape[0] + 1,dtype=np.uint16)
                mid_ind[-1] = end-1
                mid_ind[1:-1] = np.round((edges[0:-1]+edges[1:])/2)
                f = interpolate.interp1d(mid_ind, sza[i,mid_ind])
                sza_interp[i,:] = f(x_inds)


    else:
        current_time = current_time + time_delta
        continue

    if (current_time + time_delta).minute == 0:
        sio.savemat( save_folder + (current_time - datetime.timedelta(minutes=50)).strftime('%Y%m%d%H%M%S'),
                     {'sst':l3u_sst, 'sza': l3u_sza, 'time' : l3u_time})
        l3u_sst.fill(np.nan)
        l3u_time.fill(np.nan)
        l3u_sza.fill(np.nan)

    current_time = current_time + time_delta