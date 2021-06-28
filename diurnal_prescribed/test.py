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

if __name__ == '__main__':
    ##### retrieve file names & location #####
    exp_name = expname
    td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
    td_files, dy_files = [], []
    for td, dy in zip(glob.glob(td_loc + 'exp.L.Thermodynamic-??????.nc'), glob.glob(td_loc + 'exp.L.Dynamic-??????.nc')):
            td_files.append(td)
            dy_files.append(dy)
    print('Appended nc files: '+ str(len(td_files)))
    # ========================================

    for i in range(250, 350):
        print('-'*10+str(i)+'-'*10)
        data = netCDF4.Dataset(td_files[i])
        qv, qc, qr = np.array(data['qv'])[0], np.array(data['qc'])[0], np.array(data['qr'])[0]
        qi, qs, qg = np.array(data['qi'])[0], np.array(data['qs'])[0], np.array(data['qg'])[0]
        th = np.array(data['th'])[0]
        thbar = np.mean(th, axis=(1, 2))
        thbar = np.tile(thbar[:, np.newaxis, np.newaxis], (1, th.shape[1], th.shape[2]))
        g = 9.8
        arakawa_b = ((th - thbar) / thbar + 0.61 * qv - qc - qi - qr - qs - qg)
        arakawa_b1 = ((th - thbar) / thbar + 0.61 * qv)# - qc - qi - qr - qs - qg)
        # ================================================================================
        thv = th * (1 + 0.61 * qv)
        thvbar = np.mean(thv, axis=(1, 2))
        thvbar = np.tile(thvbar[:, np.newaxis, np.newaxis], (1, th.shape[1], th.shape[2]))
        my_b = (thv - thvbar) / (thvbar)
        
        # ================================================================================
        print(np.max(arakawa_b))
        print(np.max(arakawa_b1))
        print(np.max(my_b))
