#!/usr/bin/env python

# CODE TO PLOT MOMENTUM BUDGET OUTPUT
# AUTHOR: RYAN SOBASH

from matplotlib import pyplot
import numpy as np
import cPickle as pickle
from sys import argv
import scipy.ndimage as ndimage
from netCDF4 import Dataset
from wrfhelper2 import *

Ls, Lf, Le, cp, s2h = 2591, 334, 2257, 1.005, 3600.0
field, level = int(argv[1]), int(argv[2])
#st,et,ept = 28,33,156 #156
#st,et,ept = 34,39,176 #178
#st,et,ept = 40,45,195 #200
st,et,ept = 46,51,223 #223

fh = Dataset('wrfout_d01_1981-05-17_18:00:00_wsm6', 'r')
php = fh.variables['PH'][:]
phb = fh.variables['PHB'][:]

z = (phb+php)/9.8
z = (z[0,1:,0,0] + z[0,:-1,0,0])/2.0
z = z[1:-1]/1000.0

# get reflectivity and try to adjust to match moving average box
ref = computeDBZ(fh, 8, 0)[1:-1,1:-1,ept-60:ept+40]

# read in moving-box averages
TENDx = pickle.load(open('tendx%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
TENDy = pickle.load(open('tendy%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
prevp = pickle.load(open('prevp%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
psdep = pickle.load(open('psdep%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
pgmlt = pickle.load(open('pgmlt%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
psmlt = pickle.load(open('psmlt%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
u = pickle.load(open('u%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))
w = pickle.load(open('w%d-%d_%d_wsm6.pk'%(st,et,ept), 'rb'))

# convert tendencies to cooling rates (Khr-1)
prevp = prevp*Le*s2h/cp
psdep = psdep*Ls*s2h/cp
pgmlt = pgmlt*Lf*s2h/cp
psmlt = psmlt*Lf*s2h/cp

# y-average for cross-section
hw = 10
f = np.average(TENDx[field,:,level-hw:level+hw,:], axis=1)
prevp = np.average(prevp[:,level-hw:level+hw,:], axis=1)
psdep = np.average(psdep[:,level-hw:level+hw,:], axis=1)
pgmlt = np.average(pgmlt[:,level-hw:level+hw,:], axis=1)
psmlt = np.average(psmlt[:,level-hw:level+hw,:], axis=1)
u = np.average(u[:,level-hw:level+hw,:], axis=1)
w = np.average(w[:,level-hw:level+hw,:], axis=1)
ref = np.average(ref[:,level-hw:level+hw,:], axis=1)

# gaussian smoother
f = ndimage.gaussian_filter(f, sigma=0.75)
prevp = ndimage.gaussian_filter(prevp, sigma=0.75)
psdep = ndimage.gaussian_filter(psdep, sigma=0.75)

# combined cooling rate
#cooling = prevp+psdep+pgmlt+psmlt

# plot some crap
pyplot.figure(figsize=(12,6))
xvals = np.arange(0,100)
skipx = 2
skipz = 2

# PLOT REFLECTIVITY OUTLINE
cs = pyplot.contourf(xvals, z, ref, levels=[30,100], colors='0.5', alpha=0.5)

# PLOT CONTOURS FOR WIND SPEEDS
cs2 = pyplot.contour(xvals, z, u-20, levels=[0,5,10,20], colors=['0.5', '0.3', '0.1'], linewidths='2')

# PLOT MOMENTUM TENDENCIES 
cs3 = pyplot.contourf(xvals, z, f, levels=np.arange(-120,-10,20), cmap=pyplot.get_cmap('Blues_r'), alpha=0.6)
cs5 = pyplot.contourf(xvals, z, f, levels=np.arange(20,300,30), cmap=pyplot.get_cmap('Reds'), alpha=0.6)

# PLOT MICROPHYSICS TENDENCIES 
#cs3 = pyplot.contourf(xvals, z, prevp, levels=np.arange(-25,0,4), cmap=pyplot.get_cmap('Blues_r'), extend='min', alpha=0.6)
#cs3 = pyplot.contourf(xvals, z, psdep, levels=np.arange(-25,0,4), cmap=pyplot.get_cmap('Blues_r'), extend='min', alpha=0.6)
#cs3 = pyplot.contourf(xvals, z, pgmlt, levels=np.arange(-200,0,25), cmap=pyplot.get_cmap('Blues_r'), extend='min', alpha=0.6)
#cs3 = pyplot.contourf(xvals, z, psmlt, cmap=pyplot.get_cmap('Blues_r'), alpha=0.6)
#pyplot.contour(xvals, z, prevp, levels=[-1], color='blue', alpha=0.6, linewidths=0.5)
#pyplot.contour(xvals, z, psdep, levels=[-1], color='blue', linestyles='dotted')

# PLOT VERTICAL MOTION >= 1,-1
#cs2 = contourf(xvals, z, wmcross/1.94, levels=[1,100], colors='red', alpha=0.2)
#cs2 = contourf(xvals, z, wmcross/1.94, levels=[-1,-100], colors='blue', alpha=0.2)

# PLOT WIND BARBS IN PLANE OF CROSS-SECTION
pyplot.barbs(xvals[::skipx], z[::skipz], u[::skipz,::skipx]-20, w[::skipz,::skipx], length=6, barbcolor='0.4', lw=0.5, sizes=dict(emptybarb=0.05))

#cs = pyplot.contourf(np.arange(0,297), z, TENDxmax, axis=0, alpha=0.7)
#pyplot.clabel(cs2, fontsize=9, inline=1, fmt='%1.0f')

pyplot.ylim((0,12))
pyplot.xlim((20,60))
pyplot.colorbar(cs5, orientation='horizontal')
pyplot.savefig('tend.png', bbox_inches='tight')
