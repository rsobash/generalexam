#!/usr/bin/env python

from matplotlib.pyplot import *
from mpl_toolkits.basemap import *
from netCDF4 import Dataset
import sys
from datetime import datetime
from scipy.interpolate import interp1d, griddata 
from numpy import *
from wrfhelper2 import *

fh = Dataset('./genexam/wrfout_d01_1981-05-17_18:00:00_wsm6', 'r')

# GET VARIABLE TO PLOT
time = int(sys.argv[1])
level = int(sys.argv[2]);
wndlevel = 0;

print "Retrieving data from netcdf..."
gattr = fh.__dict__
u = fh.variables['U'][time,:,:,:]*1.94
v = fh.variables['V'][time,:,:,:]*1.94
w = fh.variables['W'][time,:,:]*1.94
t = fh.variables['T'][time,:,:]
nx = len(fh.dimensions['west_east'])
ny = len(fh.dimensions['south_north'])
nz = len(fh.dimensions['bottom_top'])
phb = fh.variables['PHB'][time,:,:]
php = fh.variables['PH'][time,:,:]
dx, dy = 2, 2
xpts,ypts,zpts = arange(0,nx), arange(0,ny), arange(0,nz)

z = (phb+php)/9.8
zm = (z[1:,:] + z[:-1,:])/2.0
ref = computeDBZ(fh, time, copy)

um = (u[:,:,1:] + u[:,:,:-1])/2.0
vm = (v[:,1:,:] + v[:,:-1,:])/2.0
wm = (w[1:,:,:] + w[:-1,:,:])/2.0
uvm = (um**2 + vm**2)**0.5

maxrefl = ref.max(axis=0)
#mr = round(maxrefl.max(), 0)

### PLOT DATA ###
print "Plotting data..."
figure(figsize=(12,6))
if (level == -99):
	f = maxrefl
else:
	f = average(ref[:,level-10:level+10,:], axis=1) 
	#f = ref[:,level,:]
levs = arange(10,70,5)

umcross = average(um[:,level-10:level+10,:], axis=1)
wmcross = average(wm[:,level-10:level+10,:], axis=1)
zptshgt = zm[:,0,0]/1000.0

# PLAN VIEW REFLECTIVITY AND GROUND-RELATIVE WINDS
#skipx = 15
#cs = contourf(xpts, ypts, f, levels=levs) 
#contourf(xpts, ypts, uvm[level,:]-40, levels=[0,10,20,30])
#barbs(xpts[::skipx], ypts[::skipx], um[level,::skipx,::skipx]-40, vm[level,::skipx,::skipx], length=6, barbcolor='0.4', lw=0.5, sizes=dict(emptybarb=0.05))
#axis('scaled')
#xlim((50,300))
#ylim((0,200))

# CROSS-SECTION OF REFLECTIVITY
skipx = 2
skipz = 2
cs = contourf(xpts, zptshgt, f, levels=[30,100], colors='0.5', alpha=0.5)
cs2 = contourf(xpts, zptshgt, wmcross/1.94, levels=[1,100], colors='red', alpha=0.2)
cs2 = contourf(xpts, zptshgt, wmcross/1.94, levels=[-1,-100], colors='blue', alpha=0.2)
cs = contour(xpts, zptshgt, (umcross-40)/1.94, levels=[0,5,10,20], colors=['0.5', '0.3', '0.1'], linewidths='2')
barbs(xpts[::skipx], zptshgt[::skipz], umcross[::skipz,::skipx]-40, wmcross[::skipz,::skipx], length=6, barbcolor='0.4', lw=0.5, sizes=dict(emptybarb=0.05))
xlim((240,280))
ylim((0,12))

#title('Max Reflectivity and Lowest Model Level Wind (kts)', size=12) # add a title
#cbar = colorbar(pad=0,aspect=40,orientation='horizontal')
#cbar.set_label('Reflectivity (dBZ)')
#axis('auto')
savefig('wrf.png', bbox_inches='tight')
