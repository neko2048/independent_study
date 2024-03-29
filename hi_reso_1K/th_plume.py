import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, DivergingNorm
import netCDF4
import sys
import glob
import os
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar
from nclcmaps import nclcmap
# ===========================================

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
print(np.max(z))
x = x - max(x)/2
xslice = slice(0, len(x))
z = z[z<=15000]
th = dt['th'][0, :len(z), :, :]
area = th[0, :, :] != th[0, 0, 0]
area = np.tile(area, (len(z), 1, 1))
# =============================
def average(init, end):
    total = end - init
    thmean = np.zeros((len(z), ))
    umean = np.zeros((len(z), ))
    wmean = np.zeros((len(z), ))
    wsmean = np.zeros((len(z), ))
    for n in range(init, end):
        ##### retrieve variables #####
        td = netCDF4.Dataset(td_files[n])
        progressbar(now=n-init, length=total, text='Averaging theta')
        th = td['th'][0, :len(z)] * area
        th = np.where(th == 0., np.nan, th)
        th = np.nanmean(th, axis=(1, 2))
        thbar = np.mean(td['th'][0, :len(z), :, :], axis=(1, 2))
        thmean += (th - thbar) / thbar / total
        dy = netCDF4.Dataset(dy_files[n])
        w = np.mean(dy['w'][0, :len(z)] * area, axis=(1, 2))
        wmean += w / total
    return thmean, wmean 

def draw_thp(init, end):
    ##### retrieve variables #####
    xx, zz = np.meshgrid(x, z)
    dt = netCDF4.Dataset(td_files[0])
    t = int(np.array(dt['Time']))
    thmean, wmean = average(init, end)

    ##### figure setting #####
    fig = figure(figsize=(16, 12))
    
    ##### draw th-thbar #####
    ax1 = subplot(111)
    wmean = np.expand_dims(wmean, axis=-1)
    bar1 = ax1.plot(wmean, z, color='black')
    ax1.set_title('Average', fontsize=30)
    if np.min(wmean) >= 0: 
        ax1.set_xlim(0, 0.07)
    else:
        ax1.set_xlim(np.min(wmean), 0.07)
    ax1.set_ylabel('Z (m)')
    ax1.set_xlabel('Average Vertical Velocity')
    ##### draw w #####
    ax2 = ax1.twiny()
    thmean = np.expand_dims(thmean, axis=-1)
    if np.min(wmean) >= 0:
        bar2 = ax1.pcolormesh([0, 0.07], z, thmean*1e5, cmap='coolwarm', vmin=-40, vmax=40, shading='nearest')
    else:
        bar2 = ax1.pcolormesh([np.min(wmean), 0.07], z, thmean*1e5, cmap='coolwarm', vmin=-40, vmax=40, shading='nearest')
    cbar1 = fig.colorbar(bar2, extend='both')
    cbar1.set_label('Buoyancy ($x 10^{-5}$)', rotation=270, fontsize=15)
    ax2.set_xticklabels([])
    savefig('th_average.png', dpi=300)
    clf()

if __name__ == '__main__':
    init, end = 11, 30 # from bubble up to dissipated
    draw_thp(init, end)
