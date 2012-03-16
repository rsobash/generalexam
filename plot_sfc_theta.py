#!/usr/bin/env python

# CODE TO PLOT SURFACE WIND AREAS
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
from matplotlib.pyplot import *
import cPickle as pickle

# READ IN NECESSARY 4D MODEL FIELDS
#fh = Dataset("wrfout_d01_1981-05-17_18:00:00", "r")
fh = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6", "r")
fh2 = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6enhg", "r")
tlin    = fh.variables['T'][:,:]+300
twsm6    = fh2.variables['T'][:,:]+300

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

times = np.arange(6,7,2)
for t in times:
	#contour(pertlin[t,0,:], levels=[-1], colors='blue', linestyles='solid', linewidths=0.5)
	#contour(pertwsm6[t,0,:], levels=[-1], colors='red', linestyles='solid', linewidths=0.5)
	contourf(pertwsm6[t,0,:])
	#contourf(maglin[t,0,:], levels=[25,100], colors='blue', alpha=0.5)
	#contourf(magwsm6[t,0,:], levels=[25,100], colors='red', alpha=0.5)
	
savefig('test.png')
