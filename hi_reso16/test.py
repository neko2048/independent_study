import netCDF4 
from matplotlib.pyplot import *
import numpy as np 
import sys
import glob
import os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar, buoyancy
expname = 'hi_reso16'

data = np.array(netCDF4.Dataset('./nc_files/data000000.nc')['buoyancy'])
x = np.array(netCDF4.Dataset('./nc_files/data000000.nc')['xc'])
y = np.array(netCDF4.Dataset('./nc_files/data000000.nc')['yc'])
z = np.array(netCDF4.Dataset('./nc_files/data000000.nc')['zc'])

#contourf(x, z, data[0, :, 63, :]);colorbar()
#savefig('test.jpg')
