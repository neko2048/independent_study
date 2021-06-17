import numpy as np
import netCDF4


td_loc = '/media/wk2/atmenu10246/VVM/DATA/hi_reso08/archive/'
td0 = netCDF4.Dataset(td_loc+'hi_reso08.L.Thermodynamic-000000.nc')
thbar = td0['th'][0, :, 5, 5]
x = td0['xc']
z = td0['zc']
th = td0['th'][0, :, :, 63] - np.tile(thbar, (len(x), 1)).transpose()
rmax = 0
for i in range(len(z)):
	if np.max(th[i]) == 0:
		continue	
	else:
		radius = np.max(x[th[i] != 0]) - np.min(x[th[i] != 0])
	print(radius)
	print(np.argmax(x[th[i] != 0]), np.argmin(x[th[i] != 0]))
