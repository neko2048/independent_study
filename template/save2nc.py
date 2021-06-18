import xarray as xa
import numpy as np
import netCDF4
import glob, sys, os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar, buoyancy

def gen_dwdt(tidx):
    if tidx == 0:
        prevw = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0]
        nextw = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0]
        dwdt = (nextw - prevw) / (60)

    elif tidx == len(dy_files)-1:
        prevw = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0]
        nextw = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0]
        dwdt = (nextw - prevw) / (60)
    else:
        prevw = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0]
        nextw = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0]
        dwdt = (nextw - prevw) / (120)

    return dwdt

if __name__ == "__main__":
    ##### appending data location #####
    exp_name = str(input())
    td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
    td_files, dy_files = [], []
    for td, dy in zip(glob.glob(td_loc + exp_name + '.L.Thermodynamic-??????.nc'), glob.glob(td_loc + exp_name + '.L.Dynamic-??????.nc')):
        td_files.append(td)
        dy_files.append(dy)
    print('Appended nc files: '+ str(len(td_files)))
    # ========================================

    ##### universal variables #####
    dt = netCDF4.Dataset(td_files[0]) # take initial state
    z, y, x = dt['zc'], dt['yc'], dt['xc']
    
    for tidx in range(len(td_files)):
        progressbar(now=tidx, length=len(td_files), text='generating nc file')
        ##### buoyancy #####
        qv = np.array(netCDF4.Dataset(td_files[tidx])['qv'])[0]
        th = np.array(netCDF4.Dataset(td_files[tidx])['th'])[0]
        buo = buoyancy(th, qv)

        ##### dwdt #####
        dwdt = gen_dwdt(tidx)

        buo, dwdt = buo[np.newaxis, :], dwdt[np.newaxis]
        data = xa.Dataset({'buoyancy': (('Time', 'zc', 'yc', 'xc'), buo), 
                           'dwdt': (('Time', 'zc', 'yc', 'xc'), dwdt)}, 
                           coords={
                           'Time': [tidx], 
                           'zc': np.array(z),
                           'yc': np.array(y),
                           'xc': np.array(x)
                           })
        data.to_netcdf('./nc_files/data{0:06d}.nc'.format(tidx))
