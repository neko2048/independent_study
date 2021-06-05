import numpy as np 
from matplotlib.pyplot import *
import netCDF4
import sys
import os
import glob
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar
# ===========================================

##### retrieve file names & location #####
cd_files, cr_files = [], []
for cd, cr in zip(glob.glob('./cloud/cloud??????.npy'), glob.glob('./core/core??????.npy')):
	cd_files.append(cd)
	cr_files.append(cr)
print('Appended npy files: '+ str(len(cd_files)))
# ========================================

##### retrieve global variables
exp_name = str(input())
z = np.array(netCDF4.Dataset('/media/wk2/atmenu10246/VVM/DATA/'+exp_name+'/archive/'+exp_name+'.L.Thermodynamic-000000.nc')['zc'])
z = z[z<=15000]
i = 0
# =============================

def draw_core(init, end):
	"""
	input: init(int), end(init)
	return: figure(.png)

	return time series figure containing 
	buoyancy average filtered by core 
	"""
	core_collect = np.zeros((len(z), end-init))
	for tidx in range(init, end):
		progressbar(now=tidx, length=end, text='core progress')
		core = np.load('./core/core%06d.npy'%tidx, allow_pickle=True) # (z, y, x)
		core_collect[:, tidx-init] = np.mean(core, axis=(1, 2))[:len(z)]
	contourf(np.arange(init, end), z, core_collect, vmin=-0.01, vmax=0.01, levels=50, cmap='coolwarm')
	colorbar()
	title('Core Bouyancy')
	xlabel('Time (min)')
	ylabel('Height (m)')
	savefig('Tcore.png', dpi=300)
	clf()

def draw_cloud(init, end):
	"""
	input: init(int), end(init)
	return: figure(.png)

	return time series figure containing 
	buoyancy average filtered by cloud 
	"""
	cloud_collect = np.zeros((len(z), end-init))
	for tidx in range(init, end):
		progressbar(now=tidx, length=end, text='cloud progress')
		cloud = np.load('./cloud/cloud%06d.npy'%tidx, allow_pickle=True) # (z, y, x)
		cloud_collect[:, tidx-init] = np.mean(cloud, axis=(1, 2))[:len(z)]
	contourf(np.arange(init, end), z, cloud_collect, vmin=-0.01, vmax=0.01, levels=50, cmap='coolwarm')
	colorbar()
	title('Cloud Buoyancy')
	xlabel('Time (min)')
	ylabel('Height (m)')
	savefig('Tcloud.png', dpi=300)
	clf()

if __name__ == '__main__':
	init, end = 0, len(cd_files)
	draw_core(init, end)
	draw_cloud(init, end)
