import numpy as np
import netCDF4
from matplotlib.pyplot import *
expname = 'hi_reso16'

dt = netCDF4.Dataset("/media/wk2/atmenu10246/VVM/DATA/"+expname+"/archive/"+expname+".L.Thermodynamic-000000.nc")

z, theta, vapor = dt['zc'], dt['th'], dt['qv']
print(np.array(z))

#fig, axs = subplots(1, 2, sharey=True)
#axs[0].plot(theta[0, :, 0, 0], z, color='#a22041')
#axs[0].set_title('Potential Temperature [K]')
#axs[0].grid(axis='y')
#axs[1].plot(vapor[0, :, 0, 0], z, color='#0f2350')
#axs[1].set_title('Vapor Mixing Ratio [kg/kg]')
#axs[1].grid(axis='y')
imshow([[-5.0, 5.0], [-3.0, 3.0]], cmap='coolwarm', vmin=-5, vmax=5)
colorbar()
#ylim(-100, 16000)
savefig('theta.jpg', dpi=400)
