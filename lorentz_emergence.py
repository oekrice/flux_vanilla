#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Script for plotting the Lorentz Force at given heights in the emergence sims
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.ndimage import gaussian_filter1d
#from fltrace import trace_fieldlines

#matplotlib.rcParams['text.usetex'] = True


if len(sys.argv) > 1:
    plot_id = int(sys.argv[1])
else:
    plot_id = 0


fig_width = 15#1.0*(513.11743/72)

nx = 128
ny = 128
nz = 128

xs = np.linspace(-130,130, nx+1)
ys = np.linspace(-130,130, ny+1)
zs = np.linspace(-25,100, nz+1)

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

photo_height = 10.0
corona_height = 20.0

z_photo = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))
z_corona = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))

if not os.path.exists('./lorentz/'):
    os.mkdir('./lorentz/')

class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

data_sources = ['./Data_150/','./Data_15/']
data_titles = ['More stratified', 'Less stratified']

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(10, 10))
#Columns are the stratification, rows are the height
vmaxs = []; vmins = []

for di, data_source in enumerate(data_sources):

    plot_heights = [20.0, 30.0, 40.0, 50.0]
    data_directory = data_source

    slice_index = ny//2
    wait = 0
    fname = '%s%04d.nc' % (data_directory, plot_id)
    print('Analysing plot', plot_id, 'fname', fname)

    try:
        data = netcdf_file(fname, 'r', mmap=False)
        print('File', fname, 'found')

    except:
        print('File', fname, 'not found')
        continue

    bx = np.zeros((nx+1,ny+2,nz+2))
    by = np.zeros((nx+2,ny+1,nz+2))
    bz = np.zeros((nx+2,ny+2,nz+1))

    bx[:,1:-1,1:-1] = np.swapaxes(data.variables['bx'][:],0,2)
    by[1:-1,:,1:-1] = np.swapaxes(data.variables['by'][:],0,2)
    bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)

    data.close()

    #Calculate Lorentz Force everywhere. Let's average to grid CENTRES because why not?
    jx = (bz[1:-1,1:,:] - bz[1:-1,:-1,:])/dy - (by[1:-1,:,1:] - by[1:-1,:,:-1])/dz
    jy =  (bx[:,1:-1,1:] - bx[:,1:-1,:-1])/dz - (bz[1:,1:-1,:] - bz[:-1,1:-1,:])/dx
    jz =  (by[1:,:,1:-1] - by[:-1,:,1:-1])/dx - (bx[:,1:,1:-1] - bx[:,:-1,1:-1])/dy

    #Average magnetic field to centres
    bx1 = 0.5*(bx[1:,1:-1,1:-1] + bx[:-1,1:-1,1:-1])
    by1 = 0.5*(by[1:-1,1:,1:-1] + by[1:-1,:-1,1:-1])
    bz1 = 0.5*(bz[1:-1,1:-1,1:] + bz[1:-1,1:-1,:-1])

    #Average current to centres
    jx1 = 0.25*(jx[:,1:,1:] + jx[:,1:,:-1] + jx[:,:-1,1:] + jx[:,:-1,:-1])
    jy1 = 0.25*(jy[1:,:,1:] + jy[1:,:,:-1] + jy[:-1,:,1:] + jy[:-1,:,:-1])
    jz1 = 0.25*(jz[1:,1:,:] + jz[1:,:-1,:] + jz[:-1,1:,:] + jz[:-1,:-1,:])

    lx1 = jy1*bz1 - jz1*by1
    ly1 = jz1*bx1 - jx1*bz1
    lz1 = jx1*by1 - jy1*bx1

    l2 = lx1**2 + ly1**2 + lz1**2

    for hi, height in enumerate(plot_heights):
        h_index = int(nz*(height - zs[0])/(zs[-1] - zs[0]))
        toplot =  np.log(l2[:,:,h_index].T)
        if di == 0:
            vmins.append(np.min(toplot)); vmaxs.append(np.max(toplot))
        print(vmins)
        im = axes[hi, di].pcolormesh(xs, ys, toplot, vmin = vmins[hi], vmax = vmaxs[hi])
        fig.colorbar(im, ax = axes[hi,di],label = 'Lorentz force (log)')
        axes[hi, di].set_title('%s, height = %.1f' % (data_titles[di], height))

plt.suptitle('Lorentz Force Magnitude with Height, snap id = %d' % plot_id)
plt.tight_layout()
plt.savefig('lorentz/lorentz%03d.png' % plot_id)
plt.show()






