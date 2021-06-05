import numpy as np
import netCDF4
import sys
import glob
##### retrieve file names & location #####
exp_name = 'walter'
td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
td_files, dy_files = [], []
for td, dy in zip(glob.glob(td_loc +exp_name + '.L.Thermodynamic-??????.nc'), glob.glob(td_loc + exp_name + '.L.Dynamic-??????.nc')):
        td_files.append(td)
        dy_files.append(dy)
print('Appended nc files: '+ str(len(td_files)))
qc_max, qr_max = np.zeros((len(td_files))), np.zeros((len(td_files)))

for i in range(len(td_files)):
        qc = netCDF4.Dataset(td_files[i])['qc']
        qr = netCDF4.Dataset(td_files[i])['qr']
        qc_max[i] = np.max(qc)
        qr_max[i] = np.max(qr)

print('cloud max: %f'%np.max(qc_max))
print('rain max:  %f'%np.max(qr_max))
