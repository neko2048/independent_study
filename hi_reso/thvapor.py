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
thbar = np.mean(dt['th'][0, :len(z), :, :], axis=(1, 2))
thbar = np.tile(thbar, (len(x), 1)).transpose()
y_prof =  63# max thbar index
# =============================

def draw_thp(init, end):
	for n in range(init, end):
		##### retrieve variables #####
		xx, zz = np.meshgrid(x, z)
		td = netCDF4.Dataset(td_files[n])
		#print(td_files[n])
		progressbar(now=n, length=len(td_files), text='draw theta')
		t = int(np.array(td['Time']))
		th = td['th'][0, :len(z), y_prof, xslice]
		thbar = np.mean(td['th'][0, :len(z), :, :], axis=(1, 2))
		thbar = np.tile(thbar, (len(x), 1)).transpose()
		qc = td['qc'][0, :len(z), y_prof, xslice]
		qr = td['qr'][0, :len(z), y_prof, xslice]
		dy = netCDF4.Dataset(dy_files[n])
		u, w = dy['u'][0, :len(z), y_prof, xslice], dy['w'][0, :len(z), y_prof, xslice]
		ws = np.sqrt(u**2 + w ** 2)
		##### draw th-thbar #####
		contourf(x, z, th-thbar, levels=50, vmin=-5, vmax=5, cmap='coolwarm')
		#print(np.max(th-thbar))
		colorbar(extend='both')
		##### draw cloud and rain #####
		#contour(x, z, qc-1e-6, colors='grey')
		#contour(x, z, qr-1e-6, colors='blue')

		##### wind blob #####
		interval = 3
		quiver(xx[::interval, ::interval], zz[::interval, ::interval], 
			(u/ws)[::interval, ::interval], (w/ws)[::interval, ::interval],
			ws[::interval, ::interval], cmap='rainbow')
		##### figure setting #####
		title('Thermal Bubble @ t = '+ str(t) + ' min')
		xlabel('X (m)')
		ylabel('Z (m)')
		savefig('th'+td_files[n][-9:-3]+'.png', dpi=300)
		clf()

def draw_vapor(init, end):
	n = 0
	for f in td_files[init:end]:
		#print(f)
		progressbar(now=n, length=len(td_files), text='draw vapor')
		n += 1
		##### retrieve variables #####
		dt = netCDF4.Dataset(f)
		t = int(np.array(dt['Time']))
		qc = dt['qc'][0, :len(z), y_prof, xslice]
		#qv = dt['qv'][0, :len(z), y_prof, :]
		qr = dt['qr'][0, :len(z), y_prof, xslice]
		##### draw qv #####
		#qv = np.ma.masked_array(qv, qv<1e-10)
		#contourf(x, z, qv, levels=50, vmin=0., vmax=0.02, cmap='Blues')
		#colorbar(extend='max')
                ##### draw qc #####
		colors = cm.Greys(np.hstack([np.array([0.]*5), np.linspace(0.5, 0.75, 95)]))
		cmap = LinearSegmentedColormap.from_list('name', colors)
		#cmap.set_bad('white')
		qc = np.ma.masked_array(qc, qc<0.0)
		contourf(x, z, qc, cmap=cmap, vmin=0, vmax=0.001, levels=100)
		colorbar()
		##### draw qr #####
		cmap = nclcmap('BrownBlue12')
		colors = cm.Blues(np.linspace(0.5, 1))
		cmap = LinearSegmentedColormap.from_list('name', colors)
		qr = np.ma.masked_array(qr, qr<=1e-6)#5e-6)
		contourf(x, z, qr, cmap=cmap, vmin=0., vmax=0.01, alpha=.5, levels=20)
		##### figure setting #####
		title('Vapor Distribution @ t = '+ str(t) + ' min')
		xlabel('X (m)')
		ylabel('Z (m)')
		savefig('vapor'+f[-9:-3]+'.png', dpi=300)
		clf()


if __name__ == '__main__':
	init, end = 0, len(td_files)
	draw_thp(init, end)
	draw_vapor(init, end)
