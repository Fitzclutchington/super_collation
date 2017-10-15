"""
Compute the time_series of clear-sky samples for each pixel within a list
of L3U files.
"""
import glob
import sys
import os
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import utils

class TimeSeries:

    
    # TODO: preform all opeartions on ABI crop

    def __init__(self, l3u_geo_granule_list, l3u_polar_granule_list, pixel_list, nc_variable_name):
        self.l3u_geo_granule_list = l3u_geo_granule_list
        self.l3u_polar_granule_list = l3u_polar_granule_list
        self.nc_variable_name = nc_variable_name
        self.pixel_list = pixel_list
        self.satellite_list = self.create_satellite_list()
        self.data = self.create_data_dict()

    def create_satellite_list(self):
        
        satellite_list = set()
        for filepath in self.l3u_polar_granule_list:
            filename = filepath.split('/')[-1]
            satellite_name = filename.split('-')[4]
            satellite_list.add(satellite_name)
        return list(satellite_list)

    def create_data_dict(self):

        variables = ['sst', 'time', 'sses', 'sza']
        data = { 'geo' : {},
                 'polar' : {}}
        

        for category in data:
            for variable in variables:
                data[category][variable] = {}
                for pixel in map(str,self.pixel_list):
                    if category == 'geo':
                        data[category][variable][pixel] = []
                    else:
                        data[category][variable][pixel] = {}
        
        for variable in variables:
            for pixel in map(str,self.pixel_list):
                for satellite_name in self.satellite_list:
                    data['polar'][variable][pixel][satellite_name] = []
        return data

    def extract_time_series(self, vmin, vmax):

        num_files = len(self.l3u_geo_granule_list)

        
        for i,filepath in enumerate(self.l3u_geo_granule_list):

            filename = filepath.split('/')[-1]
            satellite_name = filename.split('-')[4]

            sst = utils.read_var(filepath, self.nc_variable_name)
            dtime = utils.read_var(filepath, 'sst_dtime')
            sza = utils.read_var(filepath, 'satellite_zenith_angle')
            sses = utils.read_var(filepath, 'sses_bias')

            for pixel in self.pixel_list:
                row, col = pixel
                p = str(pixel)
                sst_val = sst[row,col]
                time_val = float(filename[8:10])  + float(filename[10:12])/60 + (float(filename[12:14]) + dtime[row,col])/(60*60)
                sza_val = sza[row,col]
                sses_val = sza[row,col]

                if sst_val >=vmin and sst_val <= vmax:
                    self.data['geo']['sst'][p].append(sst_val)
                else:
                    self.data['geo']['sst'][p].append(np.nan)
                self.data['geo']['time'][p].append(time_val)
                self.data['geo']['sses'][p].append(sses_val)
                self.data['geo']['sza'][p].append(sza_val)

            print i, '/', num_files, 'geo'

        num_files = len(self.l3u_polar_granule_list)

        for i,filepath in enumerate(self.l3u_polar_granule_list):

            filename = filepath.split('/')[-1]
            satellite_name = filepath.split('-')[4]

            sst = utils.read_var(filepath, self.nc_variable_name)
            dtime = utils.read_var(filepath, 'sst_dtime')
            sza = utils.read_var(filepath, 'satellite_zenith_angle')
            sses = utils.read_var(filepath, 'sses_bias')

            for pixel in self.pixel_list:
                row = pixel[0]
                col = pixel[1]
                p = str(pixel)
                sst_val = sst[row,col]
                time_val = float(filename[8:10])  + float(filename[10:12])/60 + (float(filename[12:14]) + dtime[row,col])/(60*60)
                sza_val = sza[row,col]
                sses_val = sses[row,col]

                if sst_val >=vmin and sst_val <= vmax:
                    self.data['polar']['sst'][p][satellite_name].append(sst_val)
                else:
                    self.data['polar']['sst'][p][satellite_name].append(np.nan)
                self.data['polar']['time'][p][satellite_name].append(time_val)
                self.data['polar']['sza'][p][satellite_name].append(sza_val)
                self.data['polar']['sses'][p][satellite_name].append(sses_val)
            print i, '/', num_files, 'polar'

    def display_time_series(self):

        for pixel in map(str,self.pixel_list):
            title = ("SST time series at", str(pixel))
            fig = plt.figure()
            fig.canvas.set_window_title(pixel)
            plt.title(title)
            plt.plot(self.data['geo']['time'][pixel], self.data['geo']['sst'][pixel], c='k', label='SST ABI')
            plt.plot(self.data['geo']['time'][pixel], self.data['geo']['sst'][pixel], 'k.')
            colors = ['m','b','g','c','r','#ffa500', '#551a8b']
            for satellite_name,c in zip(self.satellite_list,colors):
                plt.plot(self.data['polar']['time'][pixel][satellite_name], self.data['polar']['sst'][pixel][satellite_name], '.', markersize=20, c=c, label=satellite_name)
            plt.legend()
            plt.show()

#def main():


def display_time_series(ts):

    for pixel in map(str,ts.pixel_list):
        title = ("SST time series at", str(pixel))
        fig = plt.figure()
        plt.title(title)
        fig.canvas.set_window_title(pixel)
        plt.plot(ts.data['geo']['time'][pixel], ts.data['geo']['sst'][pixel], c='k', label='SST ABI')
        plt.plot(ts.data['geo']['time'][pixel], ts.data['geo']['sst'][pixel], 'k.')
        colors = ['m','b','g','c','r','#ffa500', '#551a8b']
        for satellite_name,c in zip(ts.satellite_list,colors):
            plt.plot(ts.data['polar']['time'][pixel][satellite_name], ts.data['polar']['sst'][pixel][satellite_name], '.', markersize=20, c=c, label=satellite_name)
        plt.legend()
        plt.show()

l3u_geo_folder = sys.argv[1]
l3u_polar_folder = sys.argv[2]
pixel_list_file = sys.argv[3]

pixel_list = utils.read_pixels(pixel_list_file)
l3u_geo_file_list = sorted(glob.glob((os.path.join(l3u_geo_folder,'*.nc'))))
l3u_polar_file_list = sorted(glob.glob((os.path.join(l3u_polar_folder,'*.nc'))))

l3u_variable = 'sea_surface_temperature'

vmin = 271.15
vmax = 310

time_series = TimeSeries(l3u_geo_file_list, l3u_polar_file_list, pixel_list, 'sea_surface_temperature')
time_series.extract_time_series(vmin,vmax)
time_series.display_time_series()
sio.savemat("time_series_data.mat",time_series.data)

"""
if __name__ == '__main__':
    main()
"""