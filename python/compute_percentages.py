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
    

    def __init__(self, data_files):
        self.data_files = data_files

    def compute_percentage(self):

        self.data = {}

        for data_file in self.data_files:
            data = sio.loadmat(data_file)
            day = data['quality_frequency_day']
            night = data['quality_frequency_night']
            self.data[data_file] = {}
            self.data[data_file]['day'] = (np.isfinite(day)).sum()/float(day.shape[0]*day.shape[1])
            self.data[data_file]['night'] = (np.isfinite(night)).sum()/float(day.shape[0]*day.shape[1])

    def print_percentages(self):
        with open("counts.txt", "w") as f:
            for data_file in self.data_files:
                string = data_file.split('/')[-1].split('.')[0] + ": day = " + str(self.data[data_file]['day']) + " night = " +  str(self.data[data_file]['night'])+"\n"
                f.write(string)


data_files = ["../data/overlap_data/frequency_all_sensors.mat", "../data/overlap_data/frequency_aqua_viirs.mat", "../data/overlap_data/frequency_viirs.mat"]
frequency_map = FrequencyMap(data_files)
frequency_map.compute_percentage()
frequency_map.print_percentages()
