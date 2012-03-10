#!/usr/bin/env python

# CODE TO COMPUTE MOMENTUM BUDGET FROM WRF OUTPUT
# AUTHOR: RYAN SOBASH

from netCDF4 import Dataset
import numpy as np
import cPickle as pickle
from sys import argv

st,et = int(argv[1])-1, int(argv[2])+2
# READ IN NECESSARY 4D MODEL FIELDS
#fh = Dataset("wrfout_d01_1981-05-17_18:00:00", "r")
fh = Dataset("auxhist2_d01_1981-05-17_18:00:00_wsm6enhg", "r")
u    = fh.variables['U'][st:et,:]
v    = fh.variables['V'][st:et,:]
w    = fh.variables['W'][st:et,:]
ub   = fh.variables['U_BASE'][st:et,:]
t   = fh.variables['T'][st:et,:]
pp   = fh.variables['P'][st:et,:]
pb   = fh.variables['PB'][st:et,:]
php  = fh.variables['PH'][st:et,:]
phb  = fh.variables['PHB'][st:et,:]
mub  = fh.variables['MUB'][st:et,:]
mup  = fh.variables['MU'][st:et,:]
dnw  = fh.variables['DNW'][st:et,:]
prevp = fh.variables['PREVP'][st:et,:]
psdep = fh.variables['PSDEP'][st:et,:]
pgmlt = fh.variables['PGMLT'][st:et,:]
psmlt = fh.variables['PSMLT'][st:et,:]
#dx = fh.variables['DX'][0]
#dy = fh.variables['DY'][0]
dx, dy, dt = 2000.0, 2000.0, 300.0
nxstag, nystag = 300, 250
nx, ny = nxstag-1, nystag-1
mcsspeedx = 20.0
mcsspeedy = 0.0

### PREREQUISITES ###
# TOTAL PRESSURE AT MASS POINTS
p = pp + pb
# HEIGHTS
ph = phb + php
dz = (ph[:,1:,:,:] - ph[:,:-1,:,:])/9.8
# DENSITY
mu = mub + mup
dph = ph[:,1:,:,:] - ph[:,:-1,:,:]
phi_eta = dph/dnw[:,:,np.newaxis,np.newaxis]
rho = -mu[:,np.newaxis,:,:]/phi_eta

# COMPUTE BASE/PERTURBATION U AT W-E STAGGER POINTS
ub = np.repeat(ub[:,:,np.newaxis], ny, axis=2)
ub = np.repeat(ub[:,:,:,np.newaxis], nxstag, axis=3)
up = u-ub

### DERIVATIVES FOR MOMENTUM TENDENCY TERMS###
# Du/Dy
u1 = (u[:,:,1:,:] + u[:,:,:-1,:])/2.0
u2 = (u1[:,:,:,1:] + u1[:,:,:,:-1])/2.0
dudy = (u2[:,:,1:,:] - u2[:,:,:-1,:])/dy
# Dv/Dy
dvdy = (v[:,:,1:,:] - v[:,:,:-1,:])/dy

# Du/Dx
dudx = (u[:,:,:,1:] - u[:,:,:,:-1])/dx
# Dv/Dx
v1 = (v[:,:,:,1:] + v[:,:,:,:-1])/2.0
v2 = (v1[:,:,1:,:] + v1[:,:,:-1,:])/2.0
dvdx = (v2[:,:,:,1:] - v2[:,:,:,:-1])/dx

# Du/Dt
dudt = (u[2:,:,:,:] - u[:-2,:,:,:])/(2*dt)
dudt = (dudt[:,:,:,1:] + dudt[:,:,:,:-1])/2.0
# Dv/Dt
dvdt = (v[2:,:,:,:] - v[:-2,:,:,:])/(2*dt)
dvdt = (dvdt[:,:,1:,:] + dvdt[:,:,:-1,:])/2.0

# Du'/Dz
up1 = (up[:,1:,:,:] + up[:,:-1,:,:])/2.0
up2 = (up1[:,:,:,1:] + up1[:,:,:,:-1])/2.0
dupdz = (up2[:,1:,:,:] - up2[:,:-1,:,:])/dz[:,1:-1,:,:]
# Dv'/Dz
vp1 = (v[:,1:,:,:] + v[:,:-1,:,:])/2.0
vp2 = (vp1[:,:,1:,:] + vp1[:,:,:-1,:])/2.0
dvpdz = (vp2[:,1:,:,:] - vp2[:,:-1,:,:])/dz[:,1:-1,:,:]

# Dubar/Dz
ub1 = (ub[:,1:,:,:] + ub[:,:-1,:,:])/2.0
ub2 = (ub1[:,:,:,1:] + ub1[:,:,:,:-1])/2.0
dubdz = (ub2[:,1:,:,:] - ub2[:,:-1,:,:])/dz[:,1:-1,:,:]

# dpdx
dpdx = (p[:,:,:,2:] - p[:,:,:,:-2])/(2*dx)
dpdy = (p[:,:,2:,:] - p[:,:,:-2,:])/(2*dy)

# u,v,w at tendency points
u = (u[:,:,:,1:] + u[:,:,:,:-1])/2.0
v = (v[:,:,1:,:] + v[:,:,:-1,:])/2.0
w = (w[:,1:,:,:] + w[:,:-1,:,:])/2.0

# chop off edges due to t,z,y centered differencing
u = u[1:-1,1:-1,1:-1,1:-1]
v = v[1:-1,1:-1,1:-1,1:-1]
w = w[1:-1,1:-1,1:-1,1:-1]
dudx = dudx[1:-1,1:-1,1:-1,1:-1]
dvdx = dvdx[1:-1,1:-1,1:-1,:]
dudy = dudy[1:-1,1:-1,:,1:-1]
dvdy = dvdy[1:-1,1:-1,1:-1,1:-1]
dudt = dudt[:,1:-1,1:-1,1:-1]
dvdt = dvdt[:,1:-1,1:-1,1:-1]
dupdz = dupdz[1:-1,:,1:-1,1:-1]
dvpdz = dvpdz[1:-1,:,1:-1,1:-1]
dubdz = dubdz[1:-1,:,1:-1,1:-1]
dpdx = dpdx[1:-1,1:-1,1:-1,:]
dpdy = dpdy[1:-1,1:-1,:,1:-1]
rho = rho[1:-1,1:-1,1:-1,1:-1]
print dvdt.shape, dvdx.shape, dvdy.shape, dvpdz.shape, dpdy.shape
prevp = prevp[1:-1,1:-1,1:-1,1:-1]
psdep = psdep[1:-1,1:-1,1:-1,1:-1]
pgmlt = pgmlt[1:-1,1:-1,1:-1,1:-1]
psmlt = psmlt[1:-1,1:-1,1:-1,1:-1]

print 'Computing tendency terms'
# COMPUTE U MOMENTUM TENDENCY TERMS
TENx = dudt + mcsspeedx*dudx + mcsspeedy*dudy
HAux = -(u-mcsspeedx)*dudx
HAuy = -(v-mcsspeedy)*dudy
VAup = -w*dupdz
VAub = -w*dubdz
CORx = 1e-4*v
PGAx = -(1.0/rho)*dpdx
RESx = TENx - HAux - HAuy - VAup - VAub - CORx - PGAx

# COMPUTE V MOMENTUM TENDENCY TERMS
TENy = dvdt + mcsspeedx*dvdx + mcsspeedy*dvdy
HAvx = -(u-mcsspeedx)*dvdx
HAvy = -(v-mcsspeedy)*dvdy
VAvp = -w*dvpdz
CORy = 1e-4*u
PGAy = -(1.0/rho)*dpdy
RESy = TENy - HAvx - HAvy - VAvp - CORy - PGAy

# COMBINE TENDENCIES INTO ONE BIG-ASS ARRAY
TENDx = np.array([TENx, HAux, HAuy, VAup, VAub, CORx, PGAx, RESx])
TENDy = np.array([TENy, HAvx, HAvy, VAvp, CORy, PGAy, RESy])
print 'big array has shape:', TENDx.shape

TENDx = TENDx*3600.0 #convert to ms-1 hr-1
TENDy = TENDy*3600.0 #convert to ms-1 hr-1

print 'Performing averaging'
# PERFORM AVERAGE WITH MOVING BOX
#TENDx = np.average(TENDx.reshape(TENDx.shape[0],5,7,TENDx.shape[2],TENDx.shape[3],TENDx.shape[4]), axis=1)
#TENDy = np.average(TENDx.reshape(TENDx.shape[0],5,7,TENDx.shape[2],TENDx.shape[3],TENDx.shape[4]), axis=1)
hw, cpt = 50, 125
times = np.arange(0, TENDx.shape[1])
for i in times:
        #FIND EDGE OF COLD POOL AT Y CROSS
        sfcthetapert = t[i,0,:,:] - t[i,0,0,0]
        cpt = np.max(np.nonzero(sfcthetapert[110,:]<-2)[0])
        print 'edge of sfc cold pool at index:',cpt
	try:
	  TENDx_box += TENDx[:,i,:,:,cpt-hw:cpt+hw]
	  TENDy_box += TENDy[:,i,:,:,cpt-hw:cpt+hw]
	  prevp_box += prevp[i,:,:,cpt-hw:cpt+hw]
	  psdep_box += psdep[i,:,:,cpt-hw:cpt+hw]
	  pgmlt_box += pgmlt[i,:,:,cpt-hw:cpt+hw]
	  psmlt_box += psmlt[i,:,:,cpt-hw:cpt+hw]
	  u_box += u[i,:,:,cpt-hw:cpt+hw]
	except:
	  TENDx_box = np.copy(TENDx[:,i,:,:,cpt-hw:cpt+hw])
	  TENDy_box = np.copy(TENDy[:,i,:,:,cpt-hw:cpt+hw])
	  prevp_box = np.copy(prevp[i,:,:,cpt-hw:cpt+hw])
	  psdep_box = np.copy(psdep[i,:,:,cpt-hw:cpt+hw])
	  pgmlt_box = np.copy(pgmlt[i,:,:,cpt-hw:cpt+hw])
	  psmlt_box = np.copy(psmlt[i,:,:,cpt-hw:cpt+hw])
	  u_box = np.copy(u[i,:,:,cpt-hw:cpt+hw])
	#cpt += 3

TENDx_box = TENDx_box/times.shape[0]
TENDy_box = TENDy_box/times.shape[0]
prevp_box = prevp_box/times.shape[0]
psdep_box = psdep_box/times.shape[0]
pgmlt_box = pgmlt_box/times.shape[0]
psmlt_box = psmlt_box/times.shape[0]
u_box = u_box/times.shape[0]
#TENDx = np.average(TENDx, axis=1)
#TENDy = np.average(TENDy, axis=1)

# also compute the max in this time period
print 'non-moving maximum'
TENDx = np.where(TENDx>0, 0.0, -TENDx)
TENDx_max = np.amax(TENDx[:,:,:,:,:], axis=1)
pickle.dump(TENDx_max, open('tendxmax%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
#pickle.dump(TENDx[3,:,:], open('tendxall%d-%d_wsm6enhg.pk'%(st,et-1), 'wb'))
#pickle.dump(VAup, open('VAup%d-%d_wsm6enhg.pk'%(st,et-1), 'wb'))

print 'Pickling'
pickle.dump(TENDx_box, open('tendx%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(TENDy_box, open('tendy%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(prevp_box, open('prevp%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(psdep_box, open('psdep%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(pgmlt_box, open('pgmlt%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(psmlt_box, open('psmlt%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))
pickle.dump(u_box, open('u%d-%d_wsm6enhg.pk'%(st+1,et-2), 'wb'))

