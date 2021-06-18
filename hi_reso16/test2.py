import netCDF4
import numpy as np
import matplotlib.pyplot as plt
data = netCDF4.Dataset('./nc_files/data000000.nc')
buo = data['buoyancy']

plt.contourf(buo[0, :, 63, :]);plt.colorbar()
plt.savefig('a.png', dpi=300)
