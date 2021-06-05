import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import netCDF4
import sys
import glob
import os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar
from nclcmaps import nclcmap
# ===========================================


def draw_thprof(init, end):
    for n in range(init, end):
        ##### retrieve variables #####
        xx, zz = np.meshgrid(x, z)
        td = netCDF4.Dataset(td_files[n])
        #print(td_files[n])
        progressbar(now=n, length=end-init, text='draw theta')
        t = int(np.array(td['Time']))
        th = td['th'][0, :len(z), y_prof, xslice]
        thbar = np.mean(td['th'][0, :len(z), :, :], axis=(1, 2))
        thbar = np.tile(thbar, (len(x), 1)).transpose()

        ##### dw / dt #####
        prevw = np.array(netCDF4.Dataset(dy_files[n-1 + (n == init)])['w'])
        nexw  = np.array(netCDF4.Dataset(dy_files[n+1 - (n == end-1)])['w'])
        dwdt = (nexw - prevw) / 60
        dwdt = dwdt[0, :len(z), y_prof, xslice]
        ##### draw th-thbar #####
        contourf(x, z, th-thbar, levels=50, vmin=-5, vmax=5, cmap='coolwarm')
        colorbar(extend='both')
        wt = contour(x, z, dwdt, colors='black', linewidths=.5)
        #clabel(wt, inline=True, fontsize=5)
        ##### draw wind #####
        u = np.array(netCDF4.Dataset(dy_files[n])['u'])[0, :len(z), y_prof, xslice]
        w = np.array(netCDF4.Dataset(dy_files[n])['w'])[0, :len(z), y_prof, xslice]
        ws = np.sqrt(u ** 2 + w ** 2)
        interval = 2
        quiver(xx[::interval, ::interval], zz[::interval, ::interval],
               (u/ws)[::interval, ::interval], (w/ws)[::interval, ::interval],
                ws[::interval, ::interval], cmap='rainbow')
        
        title(' t = '+ str(t) + r" min [Shade: $\theta '$, Contour: $\frac{\partial w}{\partial t}$]")
        xlabel('X (m)')
        ylabel('Z (m)')
        savefig('th_dwdt'+td_files[n][-9:-3]+'.jpg', dpi=300)
        clf()

if __name__ == '__main__':

    exp_name = str(input())
    td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
    td_files, dy_files = [], []
    for td, dy in zip(glob.glob(td_loc +exp_name + '.L.Thermodynamic-??????.nc'), glob.glob(td_loc + exp_name + '.L.Dynamic-??????.nc')):
        td_files.append(td)
        dy_files.append(dy)
    print('Appended nc files: '+ str(len(td_files)))
    # ========================================

    ##### universal variables #####
    dt = netCDF4.Dataset(td_files[0]) # take initial state
    z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])
    x = x - max(x)/2
    xslice = slice(int(len(x)/2)-2, len(x)-35)
    x = x[xslice]
    z = z[z<=12000]
    thbar = np.mean(dt['th'][0, :len(z), :, :], axis=(1, 2))
    thbar = np.tile(thbar, (len(x), 1)).transpose()
    y_prof =  63# max thbar index
    # =============================

    init, end = 0, 30#len(td_files)
    draw_thprof(init, end)

