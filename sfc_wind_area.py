#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
import cPickle as pickle
from matplotlib.pyplot import *

st,et = 6, 24
# READ IN NECESSARY 4D MODEL FIELDS
#fh = Dataset("wrfout_d01_1981-05-17_18:00:00", "r")
fh = Dataset("auxhist2_d01_1981-05-17_18:00:00_wsm6enhg", "r")
u    = fh.variables['U'][st:et,:]
v    = fh.variables['V'][st:et,:]
w    = fh.variables['W'][st:et,:]
ub   = fh.variables['U_BASE'][st:et,:]
dx, dy, dt = 2000.0, 2000.0, 300.0
nxstag, nystag = 300, 250
nx, ny = nxstag-1, nystag-1

u = (u[:,:,:,1:] + u[:,:,:,:-1])/2.0
v = (v[:,:,1:,:] + v[:,:,:-1,:])/2.0
windmag = (u**2 + v**2)**(0.5)

counts, all_counts = [], []
for thresh in [25,30,35]:
    windmagthresh = np.where(windmag>=thresh, 1, 0)
    for i in np.arange(0, windmag.shape[0]):
        counts.append(np.count_nonzero(windmagthresh[i,0,100:120,:]))
    all_counts.append(counts)
    counts = []

bar(np.arange(0, windmag.shape[0]), all_counts[0], color='b')
bar(np.arange(0, windmag.shape[0]), all_counts[1], color='y')
bar(np.arange(0, windmag.shape[0]), all_counts[2], color='r')
show()
