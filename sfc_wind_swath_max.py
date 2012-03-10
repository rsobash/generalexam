#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
from matplotlib.pyplot import *
import cPickle as pickle
from wrfhelper2 import *

fh = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6", "r")
fh2 = Dataset("auxhist2_d01_1981-05-17_18:00:00_wsm6", "r")
twsm6    = fh.variables['T'][:,:]+300
uwsm6 = fh2.variables['U'][:,:]
vwsm6 = fh2.variables['V'][:,:]
uwsm6m = (uwsm6[:,:,:,1:] + uwsm6[:,:,:,:-1])/2.0
vwsm6m = (vwsm6[:,:,1:,:] + vwsm6[:,:,:-1,:])/2.0
magwsm6 = (uwsm6m**2 + vwsm6m**2)**(0.5)

magmax = np.amax(magwsm6, axis=0)
#magmax = np.amax(uwsm6, axis=0)

refwsm6 = twsm6[:,:,0,0]
pertwsm6 = twsm6 - refwsm6[:,:,np.newaxis,np.newaxis]

#contourf(magmax[0,:], levels=[0,25,30,35,40,100], colors=['1.0','0.85', '0.7', '0.5', '0.3'])
contourf(magmax[15,:]-20, levels=[7.5,10,12.5,15,20], colors=['0.7', '0.5', '0.3','0.1'])

times = np.arange(2,13,2)
for t in times:
#	#c1 = contour(pertlin[t,0,:], levels=[-1], colors='blue', linestyles='solid', linewidths=0.5)
	c2 = contour(pertwsm6[t,0,:], levels=[-1], colors='gray', linestyles='solid', linewidths=0.5)
#	#contourf(maglin[t,0,:], levels=[25,100], colors='blue', alpha=0.5)
#	contourf(magwsm6[t,0,:], levels=[0,25,30,35,40,100], colors=['1.0','0.85', '0.7', '0.5', '0.3'])
#	#contourf(magwsm6[t,8,:]-20, levels=[0,10,20,30], colors=['0.85', '0.7', '0.5', '0.3'])

#legend((c1, c2), ('test', 'test2'))	
#cbar = colorbar(pad=0,aspect=40,orientation='horizontal')
axis('scaled')
xlim((50,300))
ylim((0,200))
savefig('windswath.png')
