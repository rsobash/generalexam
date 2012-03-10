#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
from matplotlib.pyplot import *
import cPickle as pickle
from wrfhelper2 import *

fh = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6enhg", "r")
fh2 = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6", "r")
tlin    = fh.variables['T'][:,:]+300
twsm6    = fh2.variables['T'][:,:]+300

phb = fh.variables['PHB'][:]
php = fh.variables['PH'][:]
z = (php + phb)/9.8

ulin = fh.variables['U'][:]
vlin = fh.variables['V'][:]
ulinm = (ulin[:,:,:,1:] + ulin[:,:,:,:-1])/2.0
vlinm = (vlin[:,:,1:,:] + vlin[:,:,:-1,:])/2.0
maglin = (ulinm**2 + vlinm**2)**(0.5)

uwsm6 = fh2.variables['U'][:]
vwsm6 = fh2.variables['V'][:]
uwsm6m = (uwsm6[:,:,:,1:] + uwsm6[:,:,:,:-1])/2.0
vwsm6m = (vwsm6[:,:,1:,:] + vwsm6[:,:,:-1,:])/2.0
magwsm6 = (uwsm6m**2 + vwsm6m**2)**(0.5)

reflin = tlin[:,:,0,0]
pertlin = tlin - reflin[:,:,np.newaxis,np.newaxis]
refwsm6 = twsm6[:,:,0,0]
pertwsm6 = twsm6 - refwsm6[:,:,np.newaxis,np.newaxis]

times = np.arange(2,13,2)
for t in times:
	#c1 = contour(pertlin[t,0,:], levels=[-1], colors='blue', linestyles='solid', linewidths=0.5)
	#c2 = contour(pertwsm6[t,0,:], levels=[-1], colors='gray', linestyles='solid', linewidths=0.5)
	#contourf(maglin[t,0,:], levels=[25,100], colors='blue', alpha=0.5)
	#contourf(magwsm6[t,0,:], levels=[0,25,30,35,40,100], colors=['1.0','0.85', '0.7', '0.5', '0.3'])
        contourf(magwsm6[t,15,:]-20.0, levels=[7.5,10,12.5,15,20], colors=['0.7', '0.5', '0.3','0.1'])
        print z[t,15,0,0]

#legend((c1, c2), ('test', 'test2'))	
cbar = colorbar(pad=0,aspect=40,orientation='horizontal')
axis('scaled')
xlim((50,300))
ylim((0,200))
savefig('windswath.png')
