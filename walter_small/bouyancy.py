import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from nclcmaps import nclcmap
import netCDF4
import sys
import glob
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
z = z[z<=17500]
xx, zz = np.meshgrid(x, z)
thbar = np.array(dt['th'][0, :len(z), 5, 5])
thbar = np.tile(thbar, (128, 1)).transpose()
y_prof =  63# max thbar index

def bouyancy(th, qv):
	thv = th * (1 + 0.622 * qv) # theta_v
	thvbar = np.mean(thv, (1, 2))
	thvbar = np.tile(thvbar[:, np.newaxis, np.newaxis], (1, 128, 128)) # horizontal averaged theta_v
	B = (thv - thvbar) / thvbar
	return B

##### core generation #####
def core_mask(var, bouyancy, w, qc):
	mask = np.logical_or(qc <= 0., bouyancy <= 0., w <= 0.) # positive bouyancy and W
	masked_var = np.ma.masked_array(var, mask)
	#masked_var = np.ma.masked_array(masked_var, masked_var <= 0.01*np.max(masked_var))
	return masked_var

def core_filter(t_idx, save=True):
	td, dy = netCDF4.Dataset(td_files[t_idx]), netCDF4.Dataset(dy_files[t_idx])
	print('core time: %s'%str(t_idx))
	t = int(np.array(td['Time']))
	# thermal dynamic variables
	th = td['th'][0, :, :, :]
	qv = td['qv'][0, :, :, :]
	qc = td['qc'][0, :, :, :]
	B = bouyancy(th=th, qv=qv)
	# dynamic variables
	w = dy['w'][0, :, :, :]
	# masked variable
	core = core_mask(var=B, bouyancy=B, w=w, qc=qc)
	if save:
		core.dump('./core/core%06d.npy'%t_idx)
	return core
##### core generation #####
##### core draw #####
def core_draw(init, end):
	for t_idx in range(init, end):
		core = core_filter(t_idx=t_idx, save=True)
		contourf(x, z, core[:len(z), y_prof, :], vmin=0., vmax=0.020)
		colorbar()
		title('Time: %06d'%t_idx)
		savefig('core%06d.png'%t_idx, dpi=300)
		clf()
##### core draw #####
##### cloud generation #####
def cloud_mask(var, qc):
	mask = qc <= 0. # saturated
	masked_var = np.ma.masked_array(var, mask)
	return masked_var

def cloud_filter(t_idx, save=True):
	td, dy = netCDF4.Dataset(td_files[t_idx]), netCDF4.Dataset(dy_files[t_idx])
	print('cloud time: %s'%str(t_idx))
	t = int(np.array(td['Time']))
	# thermal dynamic variables
	th = td['th'][0, :, :, :]
	qv = td['qv'][0, :, :, :]
	qc = td['qc'][0, :, :, :]
	B = bouyancy(th=th, qv=qv)
	# dynamic variables
	w = dy['w'][0, :, :, :]
	# masked variable
	cloud = cloud_mask(var=B, qc=qc)
	if save:
		cloud.dump('./cloud/cloud%06d.npy'%t_idx)
	return cloud
##### cloud generation #####
##### cloud draw #####
def cloud_draw(init, end):
	for t_idx in range(init, end):
		cloud = cloud_filter(t_idx=t_idx, save=True)
		contourf(x, z, cloud[:len(z), y_prof, :], vmin=-0.010, vmax=0.020)
		colorbar()
		title('Time: %06d'%t_idx)
		savefig('cloud%06d.png'%t_idx, dpi=300)
		clf()
##### cloud draw #####
if __name__ == '__main__':
	init, end = 0, 50#len(td_files)
	core_draw(init, end)
	cloud_draw(init, end)
