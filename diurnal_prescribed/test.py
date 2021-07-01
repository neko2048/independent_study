import netCDF4 
from matplotlib.pyplot import *
import numpy as np 
import sys
import glob
import os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar, buoyancy
expname = 'diurnal_prescribed'

data = netCDF4.Dataset('./nc_files/data000303.nc')
for i in [data['xc'], data['yc'], data['zc'], data['buoyancy']]:
    print(i.shape)


if __name__ == '__main__':
    pass
    ##### retrieve file names & location #####
    exp_name = expname
    td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
    td_files, dy_files = [], []
    for td, dy in zip(glob.glob(td_loc + 'exp.L.Thermodynamic-??????.nc'), glob.glob(td_loc + 'exp.L.Dynamic-??????.nc')):
            td_files.append(td)
            dy_files.append(dy)
    print('Appended nc files: '+ str(len(td_files)))
    # ========================================

    for i in range(0, 10):
        print('-'*10+str(i)+'-'*10)
        data = netCDF4.Dataset(dy_files[i])
        print(np.sum(np.array(data['w'])[0, 0, :, :]))
