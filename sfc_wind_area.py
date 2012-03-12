#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
import cPickle as pickle
from matplotlib.pyplot import *

st,et = 10, 70
# READ IN NECESSARY 4D MODEL FIELDS
#fh = Dataset("wrfout_d01_1981-05-17_18:00:00", "r")
fh = Dataset("auxhist2_d01_1981-05-17_18:00:00_wsm6enhg", "r")
u    = fh.variables['U'][st:et,:]
v    = fh.variables['V'][st:et,:]
#w    = fh.variables['W'][st:et,:]
dx, dy, dt = 2000.0, 2000.0, 300.0
nxstag, nystag = 300, 250
nx, ny = nxstag-1, nystag-1

stimes = np.arange(10,70,6)
tmax = []
for t in stimes:
	fname = 'vaup%d-%d_wsm6enhg.pk'%(t,t+5)
	print fname
	this = pickle.load(open(fname, 'rb'))
	for i in np.arange(0,6):
		if t == 10: ind = 8
		else: ind = 13
		thismax = this[i,:ind,:].max()
		tmax.append(thismax)

u = (u[:,:,:,1:] + u[:,:,:,:-1])/2.0
v = (v[:,:,1:,:] + v[:,:,:-1,:])/2.0
windmag = (u**2 + v**2)**(0.5)
windmagsfc = windmag[:,0,:]

counts, all_counts, windmagmax = [], [], []
for thresh in [25,30,35]:
    windmagthresh = np.where(windmagsfc>=thresh, 1, 0)
    for i in np.arange(0, windmag.shape[0]):
        counts.append(np.count_nonzero(windmagthresh[i,100:120,:]))
    all_counts.append(counts)
    counts = []

for i in np.arange(0, windmag.shape[0]):
    windmagmax.append(windmagsfc[i,100:120,:].max())

f = figure(figsize=(12,9))
subplots_adjust(hspace=0)
ax1 = subplot(311)
ax1.plot(np.arange(0, windmag.shape[0]), all_counts[0], color='0.3', linewidth='2', label='>= 25 m s-1')
ax1.plot(np.arange(0, windmag.shape[0]), all_counts[1], color='0.6', linewidth='2', label='>= 30 m s-1')
ax1.plot(np.arange(0, windmag.shape[0]), all_counts[2], color='0.9', linewidth='2', label='>= 35 m s-1')
ax1.set_yticks(np.arange(100,600,100))
ax1.set_xticks(np.arange(2,63,12))
ax1.set_xticklabels([])
ax1.grid()

#print len(windmagmax), windmag.shape[0]
#ax3 = subplot(313)
#ax3.plot(np.arange(0, windmag.shape[0]), windmagmax, color='0.4', linewidth='2')
#ax3.set_xticks(np.arange(2,63,12))
#ax3.set_xticklabels([])
#ax3.grid()

ax2 = subplot(312)
ax2.plot(np.arange(0, windmag.shape[0]), tmax, color='0.4', linewidth='2')
ax2.set_yticks(np.arange(0,300,50))
ax2.set_xticks(np.arange(2,63,12))
ax2.set_xticklabels(np.arange(1,7))
ax2.set_ylim((0,260))
ax2.grid()


#setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
savefig('windarea.png')
