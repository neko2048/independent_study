import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from nclcmaps import nclcmap
import netCDF4
import sys
import glob
##### retrieve file names & location #####
exp_name = 'walter_3k'
td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
td_files, dy_files = [], []
for td, dy in zip(glob.glob(td_loc +exp_name + '.L.Thermodynamic-??????.nc'), glob.glob(td_loc + exp_name + '.L.Dynamic-??????.nc')):
	td_files.append(td)
	dy_files.append(dy)
print('Appended nc files: '+ str(len(td_files)))

td = netCDF4.Dataset(td_files[0])
z, y, x = np.array(td['zc']), np.array(td['yc']), np.array(td['xc'])
x -= np.max(x) / 2
y -= np.max(y) / 2
xx, yy = np.meshgrid(x, y)

level = 15
distance = np.sqrt(xx ** 2 + yy ** 2)

contourf(distance)
colorbar()
savefig('distance.png', dpi=150)
