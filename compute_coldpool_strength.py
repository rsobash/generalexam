#!/usr/bin/env python

# CODE TO COMPUTE COLD POOL STRENGTH PARAMETER C
# AND C/delu RATIO


from netCDF4 import Dataset
import numpy as np
from matplotlib.pyplot import *
from sys import argv

hr = int(argv[1])

def computeC(fh):
    ep, g = 0.622, 9.8
    php = fh.variables['PH'][:]
    phb = fh.variables['PHB'][:]
    theta = fh.variables['T'][:] + 300.0
    ph = fh.variables['PH']
    qv = fh.variables['QVAPOR'][:]
    qc = fh.variables['QCLOUD'][:]
    qi = fh.variables['QICE'][:]
    qr = fh.variables['QRAIN'][:]
    qs = fh.variables['QSNOW'][:]
    qg = fh.variables['QGRAUP'][:]

    # compute heights at mass points
    z = (phb+php)/9.8
    zmass = (z[:,1:,:,:] + z[:,:-1,:,:])/2.0
    dz = z[:,1:,:,:] - z[:,:-1,:,:] 

    # compute density potential temperature (base+pert)
    qt = qv+qc+qi+qr+qs+qg
    dentheta = theta*((1+(qv/ep))/(1+qt))
    dentheta_base = dentheta[0,:,0,0]
    dentheta_pert = dentheta - dentheta_base[np.newaxis,:,np.newaxis,np.newaxis]

    # compute potential temperature (base+pert)
    theta_base = theta[0,:,0,0]
    theta_pert = theta - theta_base[np.newaxis,:,np.newaxis,np.newaxis]
    
    # create mask where theta_pert is < -1K, make sure these points are associated with cold pool 
    cpool_mask = np.logical_or(theta_pert>-1, zmass>3000.0) # pick the values to mask (True=masked)
    integrand = (dentheta_pert/dentheta_base[np.newaxis,:,np.newaxis,np.newaxis])*dz
    integrand = np.ma.array(integrand, mask=cpool_mask) #convert to masked array so we can sum over points in cold pool using cpool_mask
    csqr = -2*g*integrand.sum(axis=1)
    return csqr**0.5

fh = Dataset("wrfout_d01_1981-05-17_18:00:00_wsm6", "r")
c = computeC(fh)

shear = 12
ratio = c/shear

cs = contourf(ratio[hr,:], levels=[0.5,1,1.5,2.0], cmap=get_cmap('Blues'), extend='both')
#cs.cmap.set_over('blue')

axis('scaled')
colorbar(aspect=20)
savefig('cpool.png')
