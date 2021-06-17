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
from function import progressbar, buoyancy
from nclcmaps import nclcmap
# ===========================================

def average(init, end):
    total = end - init
    thmean = np.zeros((len(z), ))
    umean = np.zeros((len(z), ))
    wmean = np.zeros((len(z), ))
    buomean = np.zeros((len(z), ))
    for n in range(init, end):
        ##### retrieve variables #####
        td = netCDF4.Dataset(td_files[n])
        progressbar(now=n-init, length=total, text='Averaging theta')
        th = td['th'][0, :len(z)] * area
        th = np.where(th == 0., np.nan, th)
        th = np.nanmean(th, axis=(1, 2))
        thbar = np.mean(td['th'][0, :len(z), :, :], axis=(1, 2))
        thmean += (th - thbar) / total
        buomean += (th - thbar) / thbar / total
        #th = td['th'][0, :len(z)]
        #qv = td['qv'][0, :len(z)]
        #B = buoyancy(th=th, qv=qv) * area
        #B = np.where(B == 0., np.nan, B)
        #B = np.nanmean(B, axis=(1, 2))
        #buomean += B / total
        
        dy = netCDF4.Dataset(dy_files[n])
        w = np.mean(dy['w'][0, :len(z)] * area, axis=(1, 2))
        wmean += w / total
    return thmean, buomean, wmean 

def draw_snapshot(init, end):
    """draw in snap shot"""
    ##### retrieve variables #####
    xx, zz = np.meshgrid(x, z)
    dt = netCDF4.Dataset(td_files[0])
    t = int(np.array(dt['Time']))

    thss = np.zeros((end-init, len(z)))
    buoss = np.zeros((end-init, len(z)))
    wss = np.zeros((end-init, len(z)))
    dwdtss = np.zeros((end-init, len(z)))
    dwdzss = np.zeros((end-init, len(z)-1))
    for tidx in range(init, end):
        td = netCDF4.Dataset(td_files[tidx])
        progressbar(now=tidx-init, length=end-init, text='theta in snap shot')
        thbar = np.mean(td['th'][0, :len(z), :, :], axis=(1, 2))
        th = td['th'][0, :len(z)] * area
        th = np.where(th == 0., np.nan, th)
        th = np.nanmean(th, axis=(1, 2))
        thss[tidx-init, :] = th - thbar
        buoss[tidx-init, :] = (th - thbar) / thbar

        dy = netCDF4.Dataset(dy_files[tidx])
        w = np.mean(dy['w'][0, :len(z)] * area, axis=(1, 2))
        wss[tidx-init, :] = w
        
        if tidx - (end-1): # if not at last time snap
            nextw = netCDF4.Dataset(dy_files[tidx+1])['w']
            nextw = np.mean(nextw[0, :len(z)] * area, axis=(1, 2))
            dwdtss[tidx-init, :] = (nextw - w) / 60 # dw / dt
            dwdzss[tidx-init, :] = (w[1:] - w[:-1]) / (z[1] - z[0]) # ! note that this only fit in fixed interval of z
        
    # ============= draw th spatial mean ================
    #pcolormesh(np.arange(init, end), z, thss.transpose(), 
    #           shading='near', cmap='coolwarm', vmin=-1., vmax=1.)
    #colorbar()
    #xlabel('Time (min)', fontsize=12); ylabel('Height (m)', fontsize=12)
    #title(r"$\theta '$ in spatial average")
    #savefig('thss.jpg', dpi=300)
    #clf()
    # ============= draw w spatial mean ================
    #pcolormesh(np.arange(init, end), z, wss.transpose(), 
    #           shading='near', cmap='coolwarm', vmin=-0.08, vmax=0.08)
    #colorbar()
    #xlabel('Time (min)', fontsize=12); ylabel('Height (m)', fontsize=12)
    #title(r"W in spatial average")
    #savefig('wss.jpg', dpi=300)
    #clf()
    # ============= draw buoyancy spatial mean ================
    pcolormesh(np.arange(init, end), z, buoss.transpose(),
               shading='near', cmap='coolwarm', vmin=-0.004, vmax=0.004)
    colorbar()
    #contour(np.arange(init, end), z, wss.transpose(), colors='black', linewidths=0.5)
    contour(np.arange(init, end), z, dwdtss.transpose(), colors='black', linewidths=0.5)
    xlabel('Time (min)', fontsize=12); ylabel('Height (m)', fontsize=12)
    title(r"Buoyancy [Shaded] & W [Contour] in spatial average")
    savefig('buoss.jpg', dpi=300)
    clf()
    # ============= draw dw/dt spatial mean ============
    #pcolormesh(np.arange(init, end), z, dwdtss.transpose(),
    #           shading='near', cmap='coolwarm', vmin=-0.0004, vmax=0.0004)
    #colorbar()
    #xlabel('Time (min)', fontsize=12); ylabel('Height (m)', fontsize=12)
    #title(r"$\frac{\partial w}{\partial t}$ in spatial average")
    #savefig('dwdtss.jpg', dpi=300)
    #clf()
    # ============= draw dw/dz spatial mean ============
    #pcolormesh(np.arange(init, end), z, dwdzss.transpose(),
    #           shading='near', cmap='coolwarm', vmin=-0.00015, vmax=0.00015)
    #colorbar()
    #xlabel('Time (min)', fontsize=12); ylabel('Height (m)', fontsize=12)
    #title(r"$\frac{\partial w}{\partial z}$ in spatial average")
    #savefig('dwdzss.jpg', dpi=300)
    #clf()

def draw_thp(init, end):
    """draw th in average"""
    ##### retrieve variables #####
    xx, zz = np.meshgrid(x, z)
    dt = netCDF4.Dataset(td_files[0])
    t = int(np.array(dt['Time']))
    thmean, buomean, wmean = average(init, end)

    ##### figure setting #####
    fig = figure(figsize=(16, 12))
    
    ##### draw th-thbar #####
    ax1 = subplot(111)
    wmean = np.expand_dims(wmean, axis=-1)
    bar1 = ax1.plot(wmean, z, color='black')
    ax1.set_title('Average from '+str(init) + ' to ' + str(end), fontsize=30)
    if np.min(wmean) >= 0: 
        ax1.set_xlim(0, 0.07)
    else:
        ax1.set_xlim(np.min(wmean), 0.07)
    ax1.set_ylabel('Z (m)')
    ax1.set_xlabel('Average Vertical Velocity')
    ##### draw w #####
    ax2 = ax1.twiny()
    thmean = np.expand_dims(thmean, axis=-1)
    buomean = np.expand_dims(buomean, axis=-1)
    if np.min(wmean) >= 0:
        bar2 = ax1.pcolormesh([0, 0.07], z, buomean, cmap='coolwarm', 
        vmin=-np.max(buomean), vmax=np.max(buomean), shading='nearest')
    else:
        bar2 = ax1.pcolormesh([np.min(wmean), 0.07], z, buomean, cmap='coolwarm', 
        vmin=-np.max(buomean), vmax=np.max(buomean), shading='nearest')
    #bar2 = ax2.plot(buomean, z, color='red')
    #ax2.set_xlim(-0.0015, 0.0005)
    cbar1 = fig.colorbar(bar2, extend='both')
    cbar1.set_label('Buoyancy ($x 10^{-5}$)', rotation=270, fontsize=15)
    savefig('th_average.png', dpi=300)
    clf()

def massflux_ss(init, end):
    """draw th in average"""
    ##### retrieve variables #####
    xx, zz = np.meshgrid(x, z)
    dt = netCDF4.Dataset(dy_files[0])
    t = int(np.array(dt['Time']))
    density = np.loadtxt('../density.txt')[:len(z)]
    dwdzss = np.zeros((end-init, len(z)))
    for tidx in range(init, end):
        td = netCDF4.Dataset(td_files[tidx])
        progressbar(now=tidx-init, length=end-init, text='theta in snap shot')

        dy = netCDF4.Dataset(dy_files[tidx])
        w = np.mean(dy['w'][0, :len(z)] * area, axis=(1, 2))
        dwdzss[tidx-init, :] = w * density
    
    figure()
    for tidx in range(init, end):
        plot(dwdzss[tidx-init, :], z, c=cm.jet_r((tidx-init)/(end-init)), label='%s'%tidx, alpha=0.5)
    plot(np.mean(dwdzss, axis=0), z, color='black', label='average')
    legend(fontsize=5)
    title(r'Convective Mass Flux $\rho w$')
    savefig('dwdzss_ct.jpg', dpi=300)
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
    xslice = slice(0, len(x))
    z = z[z<=15000]
    th = dt['th'][0, :len(z), :, :]
    area = th[0, :, :] != th[0, 0, 0]
    area = np.tile(area, (len(z), 1, 1))
    # =============================

    init, end = 5, 25 # from bubble up to dissipated
    
    #draw_snapshot(0, len(td_files)) # from started to end
    draw_thp(init, end)
    massflux_ss(init, end)
