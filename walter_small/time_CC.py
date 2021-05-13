import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from nclcmaps import nclcmap
import netCDF4
import sys
import glob
from time import sleep
##### retrieve file names & location #####
cd_files, cr_files = [], []
for cd, cr in zip(glob.glob('./cloud/cloud??????.npy'), glob.glob('./core/core??????.npy')):
	cd_files.append(cd)
	cr_files.append(cr)
print('Appended npy files: '+ str(len(cd_files)))

z = np.array(netCDF4.Dataset('/media/wk2/atmenu10246/VVM/DATA/walter/archive/walter.L.Thermodynamic-000000.nc')['zc'])
z = z[z<=6000]
i = 0

def progressbar(now, length, text):
	nbar = 10
	j = now / length # precent
	sys.stdout.write('\r')
	sys.stdout.write(f"[{'=' * int(nbar * j):{nbar}s}] {int(100 * j)}%  {text}")
	sys.stdout.flush()
	if now == length-1:
		sys.stdout.write('\r')
		sys.stdout.write(f"[{'=' * int(nbar * 1):{nbar}s}] {int(100 * 1)}%  {text}")
		print('\nend %s'%text)
def draw_core(init, end):
	core_collect = np.zeros((len(z), end-init))
	for tidx in range(init, end):
		progressbar(tidx, end, text='core progress')
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
	cloud_collect = np.zeros((len(z), end-init))
	for tidx in range(init, end):
		progressbar(tidx, end, text='cloud progress')
		cloud = np.load('./cloud/cloud%06d.npy'%tidx, allow_pickle=True) # (z, y, x)
		cloud_collect[:, tidx-init] = np.mean(cloud, axis=(1, 2))[:len(z)]
	contourf(np.arange(init, end), z, cloud_collect, vmin=-0.01, vmax=0.01, levels=50, cmap='coolwarm')
	colorbar()
	title('Cloud Bouyancy')
	xlabel('Time (min)')
	ylabel('Height (m)')
	savefig('Tcloud.png', dpi=300)
	clf()

if __name__ == '__main__':
	init, end = 0, len(cd_files)
	draw_core(init, end)
	draw_cloud(init, end)
