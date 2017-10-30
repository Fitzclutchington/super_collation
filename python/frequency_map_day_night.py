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
    

    def __init__(self, l3u_granule_list, nc_variable_name,land_file, sst=False):
        self.l3u_granule_list = l3u_granule_list
        self.num_files = len(l3u_granule_list)
        self.nc_variable_name = nc_variable_name
        self.sst = sst
        self.satellite_name = l3u_granule_list[0].split('/')[-1].split('-')[4]
        self.land_file = land_file
        self.create_map_day_night('l2p_flags',2)
                

    def create_map_day_night(self, l2p_flag_variable_name, land_bit):
        filename = self.l3u_granule_list[0]
        cdf = netCDF4.Dataset(filename)
        self.height = cdf.dimensions['lat'].size
        self.width = cdf.dimensions['lon'].size

        self.frequency_map_day = np.zeros((self.height, self.width),dtype=np.uint8)
        self.frequency_map_night = np.zeros((self.height, self.width),dtype=np.uint8)
        self.frequency_map_day_masked = np.zeros((self.height, self.width))
        self.frequency_map_night_masked = np.zeros((self.height, self.width))
        self.frequency_map_day_quality = np.zeros((self.height, self.width))
        self.frequency_map_night_quality = np.zeros((self.height, self.width))
        if self.sst:
            self.sst_map_day = np.full((self.height, self.width), np.nan)
            self.sst_map_night = np.full((self.height, self.width), np.nan)
    
        self.lats = np.squeeze(cdf['lat'][:])
        self.lons = np.squeeze(cdf['lon'][:])
        land_mask = sio.loadmat(self.land_file)['land_mask']
        self.land_layer = np.zeros((self.height, self.width, 4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0
        self.land_layer[land_mask==255] = [r,g,b,1]

    def compute_overlay_frequency_day_night(self):

        for i,filename in enumerate(self.l3u_granule_list):
            ql_mat = utils.read_var(filename, self.nc_variable_name)
            l2p_mat = utils.read_var(filename, 'l2p_flags')

            day = np.bitwise_and(l2p_mat,512).astype(bool)

            quality_mask = ql_mat == 5
            valid_mask = ql_mat >= 0
            day_mask_valid = (valid_mask & day)
            night_mask_valid = (valid_mask & (~day))
            day_mask_quality = (quality_mask & day)
            night_mask_quality = (quality_mask & (~day))

            
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

    def plot_mat(self, mat, cmap, day_night, max_val, title):
        step = 1000
        angle = -45
        tick_font = { 'size': 15  }

        fig, ax1 = plt.subplots(1)
        fig.canvas.set_window_title(title)
        mat_min = np.nanmin(mat)
        mat_max = np.nanmax(mat)
        norm = mpc.BoundaryNorm(np.arange(-.5 , max_val + 1 ,1), cmap.N)   

        img1 = ax1.imshow(mat, cmap=cmap, norm=norm)
        ax1.imshow(self.land_layer, interpolation='nearest')

        ax1.set_xticks(np.arange(0,self.width,step))
        ax1.set_yticks(np.arange(0,self.height,step))
        ax1.set_xticklabels(['{0:.2f}'.format(x) for x in self.lons[::step]],fontdict=tick_font, rotation=angle)
        ax1.set_yticklabels(['{0:.2f}'.format(x) for x in self.lats[::step]],fontdict=tick_font)

        if day_night == 0:
            ax1.set_title(self.satellite_name + " Day")
        else:
            ax1.set_title(self.satellite_name + " Night")

        div1 = make_axes_locatable(ax1)
        cax1 = div1.append_axes("right", size="5%", pad=0.05)
        cbar1 = plt.colorbar(img1, cax=cax1,ticks=np.linspace(0,max_val,max_val+1,dtype=np.int8))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        fig.tight_layout(pad=0.05)
        plt.show()

    def display_map_day_night(self):
        
        full_mats = [ self.frequency_map_day, self.frequency_map_night ]
        full_titles = ["full_frequency_day", "full_frequency_night"]

        masked_mats = [ self.frequency_map_day_masked, self.frequency_map_night_masked,
                        self.frequency_map_day_quality, self.frequency_map_night_quality]
        masked_titles = ["masked_frequency_day", "masked_frequency_night",
                         "quality_frequency_day", "quality_frequency_night"]

        for i,(full_mat,title) in enumerate(zip(full_mats,full_titles)):
            day_night = i%2
            cmap = plt.cm.jet
            max_val = 16
            self.plot_mat(full_mat, cmap, day_night, max_val, title)

        for i,(masked_mat,title) in enumerate(zip(masked_mats,masked_titles)):
            day_night = i%2
            cmap = plt.cm.jet
            max_val = 7
            self.plot_mat(masked_mat, cmap, day_night, max_val, title)



    def open_map(self, data_file):
        data = sio.loadmat(data_file)
        self.frequency_map_day = data['full_frequency_day']
        self.frequency_map_night = data['full_frequency_night']
        self.frequency_map_day_masked = data['masked_frequency_day']
        self.frequency_map_night_masked = data['masked_frequency_night']
        self.frequency_map_day_quality = data['quality_frequency_day']
        self.frequency_map_night_quality = data['quality_frequency_night']

    def save_map_day_night(self, save_loc):
        data = { "full_frequency_day": self.frequency_map_day,  "full_frequency_night": self.frequency_map_night,
                 "masked_frequency_day": self.frequency_map_day_masked, "masked_frequency_night": self.frequency_map_night_masked,
                 "quality_frequency_day": self.frequency_map_day_quality, "quality_frequency_night":self.frequency_map_night_quality }
        sio.savemat(save_loc, data)


if len(sys.argv) == 4:
    save_loc = sys.argv[2]
    l3u_folder = sys.argv[1]
    l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))
    land_file = sys.argv[3]

    l3u_variable = 'sea_surface_temperature'

    vmin = 271.15
    vmax = 310


    frequency_map = FrequencyMap(l3u_file_list, 'quality_level', land_file, False)
    frequency_map.compute_overlay_frequency_day_night()
    #frequency_map.display_map_day_night()
    frequency_map.save_map_day_night(save_loc)

else:
    data_file = sys.argv[2]
    l3u_folder = sys.argv[1]
    l3u_file_list = sorted(glob.glob((os.path.join(l3u_folder,'*.nc'))))
    land_file = sys.argv[3]

    l3u_variable = 'sea_surface_temperature'

    vmin = 271.15
    vmax = 310


    frequency_map = FrequencyMap(l3u_file_list, 'quality_level',land_file,False )
    #frequency_map.compute_overlay_frequency_day_night()
    frequency_map.open_map(data_file)
    frequency_map.display_map_day_night()
    #frequency_map.save_map_day_night(save_loc)