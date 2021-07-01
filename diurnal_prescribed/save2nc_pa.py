import xarray as xa
import numpy as np
import netCDF4
import glob, sys, os
from CC_buoyancy import buoyancy
import concurrent.futures
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar, buoyancy

def buoyancy(data, save=True):
    """
    input: data
    retrun buoyancy(z, y, x) [m/s^2]
    """
    th = np.array(data['th'][0, :len(z), yslice, xslice])
    qv = np.array(data['qv'][0, :len(z), yslice, xslice])
    qc = np.array(data['qc'][0, :len(z), yslice, xslice])
    qi = np.array(data['qi'][0, :len(z), yslice, xslice])
    qr = np.array(data['qr'][0, :len(z), yslice, xslice])
    qs = np.array(data['qs'][0, :len(z), yslice, xslice])
    qg = np.array(data['qg'][0, :len(z), yslice, xslice])

    B = ((th - thbar) / thbar + 0.61 * qv - qc - qr - qs - qg) * 9.81
    if save:
        np.save('buoyancy{0:06d}.npy'.format(tidx), B)

    return B

def gen_dwdt(tidx, save=True):
    if tidx == 0:
        prevw1 = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        prevw2 = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        prevw = (prevw1 + prevw2) / 2

        nextw1 = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        nextw2 = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        nextw = (nextw1 + nextw2) / 2

        dwdt = (nextw - prevw) / (60)

    elif tidx == len(dy_files)-1:
        prevw1 = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        prevw2 = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        prevw = (prevw1 + prevw2) / 2

        nextw1 = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        nextw2 = np.array(netCDF4.Dataset(dy_files[tidx])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        nextw = (nextw1 + nextw2) / 2

        dwdt = (nextw - prevw) / (60)
    else:
        prevw1 = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        prevw2 = np.array(netCDF4.Dataset(dy_files[tidx-1])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        prevw = (prevw1 + prevw2) / 2

        nextw1 = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0, :-1, yargmin:yargmax, xargmin:xargmax]
        nextw2 = np.array(netCDF4.Dataset(dy_files[tidx+1])['w'])[0, 1:, yargmin:yargmax, xargmin:xargmax]
        nextw = (nextw1 + nextw2) / 2

        dwdt = (nextw - prevw) / (120)

    if save:
        np.save('dwdt.npy{0:06d}'.format(tidx), dwdt)
    return dwdt

def execution(case_idxs):
    for tidx in case_idxs:
        #progressbar(now=tidx-case_idxs[0], length=case_idxs[-1]-case_idxs[0], text='generating nc file')
        print(tidx)
        ##### buoyancy #####
        td = netCDF4.Dataset(td_files[tidx])
        buo = buoyancy(td)
        buo = buo[:-1]
        ##### dwdt #####
        dwdt = gen_dwdt(tidx)

        buo, dwdt = buo[np.newaxis, :], dwdt[np.newaxis]
        data = xa.Dataset({'buoyancy': (('Time', 'zc', 'yc', 'xc'), buo), 
                           'dwdt': (('Time', 'zc', 'yc', 'xc'), dwdt)}, 
                           coords={
                           'Time': [tidx], 
                           'zc': np.array(z[:-1]),
                           'yc': np.array(y),
                           'xc': np.array(x)
                           })
        data.to_netcdf('./nc_files/data{0:06d}.nc'.format(tidx))


def multiprocess(case_group):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(execution, case_group)

if __name__ == "__main__":
    ##### Input setting #####
    exp_name, Tinit, Tend, Ngroup = input().split()
    exp_name = str(exp_name)
    init, end = int(Tinit), int(Tend)
    Ngroup = int(Ngroup)
    # ========================================

    ##### appending data location #####
    td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
    td_files, dy_files = [], []
    for td, dy in zip(glob.glob(td_loc + 'exp.L.Thermodynamic-??????.nc'), glob.glob(td_loc + 'exp.L.Dynamic-??????.nc')):
        td_files.append(td)
        dy_files.append(dy)
    print('Appended nc files: '+ str(len(td_files)))
    # ========================================

    ##### universal variables #####
    dt = netCDF4.Dataset(td_files[0]) # take initial state
    z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])
    xargmin, xargmax = np.argmin(abs(x - x[0])), np.argmin(abs(x - y[-1]))
    yargmin, yargmax = np.argmin(abs(y - y[0])), np.argmin(abs(y - y[-1]))
    xslice, yslice = slice(xargmin, xargmax), slice(yargmin, yargmax)
    x = x[xargmin:xargmax]
    y = y[yargmin:yargmax]
    
    thbar = np.mean(dt['th'][0], axis=(1, 2))
    thbar = np.tile(thbar[:, np.newaxis, np.newaxis], (1, len(y), len(x)))
    # ========================================

    ##### parallel parameters #####
    case_group = np.array_split([x for x in range(init, end)], Ngroup)
    
    multiprocess(case_group)
    #execution(case_group)
