import numpy as np
from matplotlib.pyplot import * 
import netCDF4
import sys
import os
import glob
##### OS setting #####
if not os.path.isdir('./buoyancy/'):
    os.mkdir('./buoyancy/')
if not os.path.isdir('./W/'):
    os.mkdir('./W/')

##### retrieve file names & location #####
exp_name = 'walter_small'
td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
td_files, dy_files = [], []
for td, dy in zip(glob.glob(td_loc +exp_name + '.L.Thermodynamic-??????.nc'), glob.glob(td_loc + exp_name + '.L.Dynamic-??????.nc')):
    td_files.append(td)
    dy_files.append(dy)
print('Appended nc files: '+ str(len(td_files)))

##### universal variables #####
dt = netCDF4.Dataset(td_files[0]) # take initial state
z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])
x = x - max(x)/2
interval = x[1]-x[0]# == y[1] - y[0]
z = z[z<=17500]
xx, zz = np.meshgrid(x, z)
thbar = np.array(dt['th'][0, :len(z), 5, 5])
thbar = np.tile(thbar, (128, 1)).transpose()
y_prof =  63# max thbar index

def buoyancy(th, qv):
    thv = th * (1 + 0.622 * qv) # theta_v
    thvbar = np.mean(thv, (1, 2))
    thvbar = np.tile(thvbar[:, np.newaxis, np.newaxis], (1, 128, 128)) # horizontal averaged theta_v
    B = (thv - thvbar) / thvbar
    return B

def h_laplace(var3d, dx, dy):
    # first order gradient using central method
    x1 = np.gradient(var3d, axis=2) / dx # d var / dx
    y1 = np.gradient(var3d, axis=1) / dy # d var / dy

    # process the boundary to periodic boundary
    x1[:, :, 0] = (var3d[:, :, 1] - var3d[:, :, -1]) / (2 * dx)
    y1[:, :, 0] = (var3d[:, 1, :] - var3d[:, -1, :]) / (2 * dy)
    x1[:, :, -1] = (var3d[:, :, 0] - var3d[:, :, -2]) / (2 * dx)
    y1[:, :, -1] = (var3d[:, 0, :] - var3d[:, -2, :]) / (2 * dy)

    # second order gradient using central method
    x2 = np.gradient(x1, axis=2) / dx # d var / dx
    y2 = np.gradient(y1, axis=1) / dy # d var / dy

    # process the boundary to periodic boundary
    x2[:, :, 0] = (x1[:, :, 1] - x1[:, :, -1]) / (2 * dx)
    y2[:, :, 0] = (y1[:, 1, :] - y1[:, -1, :]) / (2 * dy)
    x2[:, :, -1] = (x1[:, :, 0] - x1[:, :, -2]) / (2 * dx)
    y2[:, :, -1] = (y1[:, 0, :] - y1[:, -2, :]) / (2 * dy)

    return (x2 + y2)

def B_laplace(t_idx, save=True):
    td, dy = netCDF4.Dataset(td_files[t_idx]), netCDF4.Dataset(dy_files[t_idx])
    print('buoyancy time: %s'%str(t_idx))
    t = int(np.array(td['Time']))
    # thermal dynamic variables
    th = td['th'][0, :, :, :]
    qv = td['qv'][0, :, :, :]
    qc = td['qc'][0, :, :, :]
    B = buoyancy(th=th, qv=qv)
    L_buoyancy = h_laplace(var3d=B, dx=interval, dy=interval)
    if save:
        L_buoyancy.dump('./buoyancy/L_buoyancy%06d.npy'%t_idx)
    return L_buoyancy

def W_laplace(t_idx, save=True):
    dt = 60 # Assert 1 min since we know 
    if t_idx == 0 or t_idx == len(dy_files): # boundary
        w = netCDF4.Dataset(dy_files[t_idx])['w'][0, :, :, :]
        w_after = netCDF4.Dataset(dy_files[t_idx+1])['w'][0, :, :, :]
        dwdt = (w_after - w) / dt
    else:
        w = netCDF4.Dataset(dy_files[t_idx-1])['w'][0, :, :, :]
        w_after = netCDF4.Dataset(dy_files[t_idx+1])['w'][0, :, :, :]
        dwdt = (w_after - w) / (2 * dt)
    dwdt = np.array(dwdt)
    print('W time: %s'%str(t_idx))
    L_W = h_laplace(var3d=dwdt, dx=interval, dy=interval)
    if save:
        L_W.dump('./W/L_W%06d.npy'%t_idx)
    return L_W

def draw_L_buoyancy(init, end):
    for t_idx in range(init, end):
        L_buoyancy = B_laplace(t_idx=t_idx, save=True)
        contourf(x, z, L_buoyancy[:len(z), y_prof, :]*10**8, vmin=-20, vmax=20)
        cbar = colorbar()
        cbar.ax.set_ylabel(r'$\nabla ^2$ Buoyancy ($\times 10^{-8}$/$s^2 m$)', rotation=270)
        cbar.ax.yaxis.set_label_coords(5,0.5)
        title('Time: %06d'%t_idx)
        xlabel('X Domain (m)')
        ylabel('Height (m)')
        savefig('Lbuoyancy%06d.jpg'%t_idx, dpi=300)
        clf()

def draw_L_W(init, end):
    for t_idx in range(init, end):
        L_W = W_laplace(t_idx=t_idx, save=True)
        contourf(x, z, L_W[:len(z), y_prof, :]*10**8, vmin=-20, vmax=20)
        cbar = colorbar()
        cbar.ax.set_ylabel(r'$\nabla ^2 \ \frac{dw}{dt}$ ($\times 10^{-8}$/$s^2 m$)', rotation=270)
        cbar.ax.yaxis.set_label_coords(5, 0.5)
        title('Time: %06d'%t_idx)
        xlabel('X Domain (m)')
        ylabel('Height (m)')
        savefig('LW%06d.jpg'%t_idx, dpi=300)
        clf()


if __name__ == "__main__":
    init, end = 0, 50#len(td_files)
    draw_L_buoyancy(init, end)
    draw_L_W(init, end)
