import numpy as np 
from matplotlib.pyplot import *
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import netCDF4
import sys
import glob
##### import lib from the parent folder #####
sys.path.insert(1, '../') # move 2 parent folder (figures/)
from function import progressbar
from nclcmaps import nclcmap
# ===========================================




def buoyancy(th, qv):
	"""
	input: th(z, y, x) [K], qv(z, y, x) [kg/kg]
	retrun buoyancy(z, y, x) [m/s^2]
	"""
	thv = th * (1 + 0.622 * qv) # theta_v
	thvbar = np.mean(thv, (1, 2))
	thvbar = np.tile(thvbar[:, np.newaxis, np.newaxis], (1, thv.shape[1], thv.shape[2])) # horizontal averaged theta_v
	B = (thv - thvbar) / thvbar
	return B

def buoyancy_draw(init, end):
	for t_idx in range(init, end):
		progressbar(now=t_idx-init, length=end-init, text='draw_buoyancy')
		td = netCDF4.Dataset(td_files[t_idx])
		th = np.array(td['th'][0, :len(z), yslice, xslice])
		qv = np.array(td['qv'][0, :len(z), yslice, xslice])
		qc = np.array(td['qc'][0, :len(z), yslice, xslice])
		time = int(np.array(td['Time']))
		buo = buoyancy(th=th, qv=qv)
		contourf(y, z, buo[:len(z), :, x_prof], vmin=-0.01, vmax=0.01, cmap='coolwarm', 
			     levels=20)
		colorbar()
		contour(y, z, qc[:len(z), :, x_prof]>0, colors='grey', alpha=.5, linewidths=1)
		#la = contour(x, z, buo[:len(z), y_prof, :], vmin=-0.01, vmax=0.01, levels=20)
		#clabel(la, inline=True)
		title('Time: %06d'%time)
		savefig('buoyancy%06d.jpg'%time, dpi=300)
		clf()

##### core generation #####
def core_mask(var, buoyancy, w, qc):
	"""
	input: var(z, y, x), buoyancy(z, y, x) [m/s^2], w(z, y, x) [m/s], qc(z, y, x) [kg/kg]
	return masked_var(z, y, x)

	filter out var with selected condition: 
	1. Saturated (qc>0.)
	2. Positive buoyancy
	3. Positive W
	"""
	mask = np.logical_or(qc <= 0., buoyancy <= 0., w <= 0.) # positive buoyancy and W
	masked_var = np.ma.masked_array(var, mask)
	#masked_var = np.ma.masked_array(masked_var, masked_var <= 0.01*np.max(masked_var))
	return masked_var

def core_filter(t_idx, save=True):
	"""
	input: t_idx(int)
	return core(z, y, x)

	with the help of core_mask(),
	return buoyancy filtered by core
	"""
	td, dy = netCDF4.Dataset(td_files[t_idx]), netCDF4.Dataset(dy_files[t_idx])
	pprogressbar(now=t_idx-init, length=end-init, 
		    text='draw_core')
	t = int(np.array(td['Time']))
	# thermal dynamic variables
	th = td['th'][0, :, :, :]
	qv = td['qv'][0, :, :, :]
	qc = td['qc'][0, :, :, :]
	B = buoyancy(th=th, qv=qv)
	# dynamic variables
	w = dy['w'][0, :, :, :]
	# masked variable
	core = core_mask(var=B, buoyancy=B, w=w, qc=qc)
	if save:
		core.dump('./core/core%06d.npy'%t_idx)
	return core
##### core generation #####

##### core draw #####
def core_draw(init, end):
	"""
	input: init(int), end(int)
	return figure(.png)
	"""
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
	"""
	input: var(z, y, x), qc(z, y, x) [kg/kg]
	return masked_var(z, y, x)

	filter out var with selected condition:
	1. Saturated (qc>0.)
	"""
	mask = qc <= 0. # saturated
	masked_var = np.ma.masked_array(var, mask)
	return masked_var

def cloud_filter(t_idx, save=True):
	"""
	input: t_idx(int)
	return cloud(z, y, x)

	with the help of cloud_mask(),
	return buoyancy filtered by core
	"""

	td, dy = netCDF4.Dataset(td_files[t_idx]), netCDF4.Dataset(dy_files[t_idx])
	progressbar(now=t_idx-init, length=end-init, 
                    text='draw_core')
	t = int(np.array(td['Time']))
	# thermal dynamic variables
	th = td['th'][0, :, :, :]
	qv = td['qv'][0, :, :, :]
	qc = td['qc'][0, :, :, :]
	B = buoyancy(th=th, qv=qv)
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
	"""
	input: init(int), end(int)
	return figure(.png)
	"""
	for t_idx in range(init, end):
		cloud = cloud_filter(t_idx=t_idx, save=True)
		contourf(x, z, cloud[:len(z), y_prof, :], vmin=-0.010, vmax=0.020)
		colorbar()
		title('Time: %06d'%t_idx)
		savefig('cloud%06d.png'%t_idx, dpi=300)
		clf()
##### cloud draw #####


if __name__ == '__main__':
	##### retrieve file names & location #####
	exp_name = str(input())
	td_loc = '/media/wk2/atmenu10246/VVM/DATA/' + exp_name + '/archive/'
	td_files, dy_files = [], []
	for td, dy in zip(glob.glob(td_loc + 'exp.L.Thermodynamic-??????.nc'), glob.glob(td_loc + 'exp.L.Dynamic-??????.nc')):
		td_files.append(td)
		dy_files.append(dy)
	print('Appended nc files: '+ str(len(td_files)))
	# ========================================
	
	
	##### universal variables #####
	dt = netCDF4.Dataset(td_files[0]) # take initial state
	z, y, x = np.array(dt['zc']), np.array(dt['yc']), np.array(dt['xc'])

	xargmin, xargmax = np.argmin(abs(x - 39000)), np.argmin(abs(x - 43000))
	x = x[xargmin:xargmax]
	xslice = slice(xargmin, xargmax)
	x_prof = np.argmin(abs(x - 41000))

	yargmin, yargmax = np.argmin(abs(y - 30000)), np.argmin(abs(y - 45000))
	y = y[yargmin:yargmax]
	yslice = slice(yargmin, yargmax)
	y_prof =  np.argmin(abs(y - 41000))# max thbar index
	
	z = z[z<=13000]

	thbar = np.mean(dt['th'][0, :len(z), :, :], axis=(1, 2))
	thbar = np.tile(thbar, (len(x), 1)).transpose()
	# =============================
	init, end = 250, 350#len(td_files)
	#core_draw(init, end)
	#cloud_draw(init, end)
	buoyancy_draw(init, end)
