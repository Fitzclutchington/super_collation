import netCDF4
import numpy as np

def read_var(filepath, variable):
    cdf = netCDF4.Dataset(filepath)
    data = np.squeeze(cdf[variable][:])
    cdf.close()
    return data

def read_pixels(point_file):
    points = []
    with open(point_file, 'r') as f:
        for line in f:
            point = map(int,line.strip().strip('()').split(','))
            points.append((point[0],point[1]))
    return points

def get_dimensions(granule_file):
    cdf = netCDF4.Dataset(granule_file)
    height = cdf.dimensions['lat'].size
    width = cdf.dimensions['lon'].size
    return height, width