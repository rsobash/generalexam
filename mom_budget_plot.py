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
print z[0,:,0,0]
z = (z[0,1:,0,0] + z[0,:-1,0,0])/2.0
z = z[1:-1]/1000.0

field, level = int(argv[1]), int(argv[2])
st, et = 28,33

TENDx = pickle.load(open('tendx%d-%d_wsm6.pk'%(st,et), 'rb'))
TENDy = pickle.load(open('tendy%d-%d_wsm6.pk'%(st,et), 'rb'))
prevp = pickle.load(open('prevp%d-%d_wsm6.pk'%(st,et), 'rb'))*2257000
psdep = pickle.load(open('psdep%d-%d_wsm6.pk'%(st,et), 'rb'))*2591000
pgmlt = pickle.load(open('pgmlt%d-%d_wsm6.pk'%(st,et), 'rb'))*334000
psmlt = pickle.load(open('psmlt%d-%d_wsm6.pk'%(st,et), 'rb'))*334000
u = pickle.load(open('u%d-%d_wsm6.pk'%(st,et), 'rb'))
#TENDxmax = pickle.load(open('tendxmax%d-%d_wsm6.pk'%(st,et), 'rb'))
#TENDxall = pickle.load(open('tendxall%d-%d_wsm6.pk'%(st,et), 'rb'))
print prevp.shape, pgmlt.shape, pgmlt.shape, psmlt.shape, TENDx.shape, TENDy.shape

# y-average for cross-section
hw = 10
f = np.average(TENDx[field,:,level-hw:level+hw,:], axis=1)
prevp = np.average(prevp[:,level-hw:level+hw,:], axis=1)
psdep = np.average(psdep[:,level-hw:level+hw,:], axis=1)
pgmlt = np.average(pgmlt[:,level-hw:level+hw,:], axis=1)
psmlt = np.average(psmlt[:,level-hw:level+hw,:], axis=1)
u = np.average(u[:,level-hw:level+hw,:], axis=1)
#TENDxmax = np.average(TENDxmax[field,:,level-hw:level+hw,:], axis=1)
#TENDxall = np.average(TENDxall[:,:,level-hw:level+hw,:], axis=2)

# gaussian smoother
f = ndimage.gaussian_filter(f, sigma=0.75)
#TENDxmax = ndimage.gaussian_filter(TENDxmax, sigma=0.75)
prevp = ndimage.gaussian_filter(prevp, sigma=0.75)
psdep = ndimage.gaussian_filter(psdep, sigma=0.75)

#TENDsum = np.sum(TENDx[0, axis=0)
#TENDfract = TENDx/TENDsum[np.newaxis,:]
#print TENDfract.shape, TENDfract.max()
#pyplot.barbs(TENDx[0,7,::5,::5], TENDy[0,7,::5,::5])
#pyplot.contourf(np.arange(0,100), np.arange(0,247), f, levels=np.arange(-60,60,10))

pyplot.figure(figsize=(12,6))
xvals = np.arange(0,100)
cs = pyplot.contourf(xvals, z, f, levels=np.arange(-80,80,10), alpha=0.7)
#cs = pyplot.contourf(np.arange(0,297), z, TENDxmax, axis=0, alpha=0.7)
#cs2 = pyplot.contour(xvals, z, prevp, linestyles='solid', colors='gray')
#cs3 = pyplot.contour(xvals, z, psdep, colors='gray', linestyles='dashed')
#cs4 = pyplot.contour(xvals, z, pgmlt, colors='gray', linestyles='dotted')
#cs5 = pyplot.contour(xvals, z, psmlt, linestyles='dotted', colors='gray')
#cs6 = pyplot.contourf(xvals, z, u, levels=[15,20,30,40], alpha=0.3, colors=['0.7','0.5','0.3','0.1'])
#pyplot.clabel(cs2, fontsize=9, inline=1, fmt='%1.0f')
#pyplot.clabel(cs3, fontsize=9, inline=1, fmt='%1.0f')
#pyplot.clabel(cs4, fontsize=9, inline=1, fmt='%1.0f')
#pyplot.clabel(cs5, fontsize=9, inline=1, fmt='%1.0f')

pyplot.ylim((0,12))
#pyplot.xlim((20,90))
pyplot.colorbar(cs, pad=0)
pyplot.savefig('tend.png')
