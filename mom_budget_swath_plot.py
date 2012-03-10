#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from matplotlib import pyplot
import numpy as np
import cPickle as pickle
from sys import argv
import scipy.ndimage as ndimage
from netCDF4 import Dataset

fh = Dataset('wrfout_d01_1981-05-17_18:00:00_wsm6', 'r')
php = fh.variables['PH'][:]
phb = fh.variables['PHB'][:]

z = (phb+php)/9.8
z = (z[0,1:,0,0] + z[0,:-1,0,0])/2.0
z = z[1:-1]/1000.0

field, level = int(argv[1]), int(argv[2])

stimes = np.arange(6,72,6)
TENDxmax = pickle.load(open('tendxmax%d-%d_wsm6.pk'%(stimes[0],stimes[0]+5), 'rb'))
for t in stimes:
	fname = 'tendxmax%d-%d_wsm6.pk'%(t,t+5)
	print fname
	thismax = pickle.load(open(fname, 'rb'))
	TENDxmax = np.maximum(TENDxmax, thismax)

# compute average over y interval
hw = 10
TENDxmax = np.average(TENDxmax[field,:,level-hw:level+hw,:], axis=1)
#TENDxall = np.average(TENDxall[:,:,level-hw:level+hw,:], axis=2)

# gaussian smoother
TENDxmax = ndimage.gaussian_filter(TENDxmax, sigma=[0.25,1.00])

pyplot.figure(figsize=(12,6))
#cs = pyplot.contourf(xvals, z, f, levels=np.arange(-80,80,10), alpha=0.7)
cs = pyplot.contourf(np.arange(0,297), z, TENDxmax, alpha=0.7, levels=np.arange(40,300,20))

pyplot.ylim((0,12))
pyplot.xlim((50,300))
pyplot.colorbar(cs, pad=0, aspect=40)
pyplot.savefig('tend.png')
