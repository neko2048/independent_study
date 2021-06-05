import netCDF4 
from matplotlib.pyplot import *
import numpy as np 
expname = 'walter_small'

dt = netCDF4.Dataset("/media/wk2/atmenu10246/VVM/DATA/"+expname+"/archive/"+expname+".L.Thermodynamic-000000.nc")
z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])
print(np.max(x) - np.min(x))
z = z[z<= 9000]
th = np.array(dt['th'])[0]#, 5, 63, :]
level = 1
thbar = th[level, 5, 5]
th = th[level, 63, :]
cell = th - thbar
xlx = x[cell != 0.]
print(np.min(xlx), np.max(xlx))
print('radius: '+str(-np.min(xlx)+np.max(xlx)))
