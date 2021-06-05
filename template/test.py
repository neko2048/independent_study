import netCDF4 
from matplotlib.pyplot import *
import numpy as np 
import sys
import glob
import os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar, buoyancy
expname = 'hi_reso24'

dt = netCDF4.Dataset("/media/wk2/atmenu10246/VVM/DATA/"+expname+"/archive/"+expname+".L.Thermodynamic-000000.nc")
z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])
z = z[z<= 5000]
th = dt['th'][0, :len(z), :, :]
thbar = np.mean(dt['th'][0, :len(z), :, :], axis=(1, 2))

th = dt['th'][0, :len(z), :, :]
area = th[0, :, :] != th[0, 0, 0]

dt = netCDF4.Dataset("/media/wk2/atmenu10246/VVM/DATA/"+expname+"/archive/"+expname+".L.Thermodynamic-000008.nc")
th = dt['th'][0, :len(z), :, :]
thbar = np.mean(dt['th'][0, :len(z), :, :], axis=(1, 2))
qv = dt['qv'][0, :len(z), :, :]
def Draw():
    B = buoyancy(th, qv)
    for i in range(len(z)):
        figure(figsize=(10, 10))
        contourf(x, x, B[i] * area, cmap='coolwarm', vmin=-0.01, vmax=0.01)
        colorbar()
        title(str(z[i]))
        savefig('tharea%s.jpg'%z[i])
        clf()

Draw()
