#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from fltrace import trace_fieldlines
from scipy.ndimage import gaussian_filter1d

#matplotlib.rcParams['text.usetex'] = True


if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0


fig_width = 15#1.0*(513.11743/72)


nx= 128; ny = 128; nz = 128

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


class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

#data_sources = ['./Data_unstratified_64/', './Data_stratified_64/']
#data_sources = ['./Data/', './Data/']

data_sources = ['./Data_150/','./Data_50/','./Data_15/']
i = 0

if len(sys.argv) > 2:
    i = int(sys.argv[2])

snap_num = 0

cs = ['green', 'red', 'blue', 'orange']

photo_height = 10.0
corona_height = 20.0

z_photo = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))
z_corona = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))

for plot_num in range(300, 301):
    ymax = 0.0; ymin = 1e6
    for source_number in range(3):

        data_directory = data_sources[source_number]

        slice_index = ny//2
        i = plot_num
        wait = 0
        fname = '%s%04d.nc' % (data_directory, i)
        print('Examining plot', i, 'fname', fname)

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

        en = np.zeros((nx+2,ny+2,nz+2))
        rho = np.zeros((nx+2,ny+2,nz+2))

        en[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['en'][:],0,2)
        rho[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['rho'][:],0,2)

        vx = np.zeros((nx+1,ny+1,nz+1))
        vy = np.zeros((nx+1,ny+1,nz+1))
        vz = np.zeros((nx+1,ny+1,nz+1))

        vx = np.swapaxes(data.variables['vx'][:],0,2)
        vy = np.swapaxes(data.variables['vy'][:],0,2)
        vz = np.swapaxes(data.variables['vz'][:],0,2)

        pr = rho*en*(2/3)

        data.close()

        def current_test(bx,by,bz):
            jx = (bz[1:-1,1:,:] - bz[1:-1,:-1,:])/dy - (by[1:-1,:,1:] - by[1:-1,:,:-1])/dz
            jy =  (bx[:,1:-1,1:] - bx[:,1:-1,:-1])/dz - (bz[1:,1:-1,:] - bz[:-1,1:-1,:])/dx
            jz =  (by[1:,:,1:-1] - by[:-1,:,1:-1])/dx - (bx[:,1:,1:-1] - bx[:,:-1,1:-1])/dy
            return jx, jy, jz

        def magfield(bx, by, bz):
            bx1 = 0.5*(bx[1:,slice_index,1:-1] + bx[:-1,slice_index,1:-1])
            by1 = 0.5*(by[1:-1,slice_index,1:-1] + by[1:-1,slice_index,1:-1])
            bz1 = 0.5*(bz[1:-1,slice_index,1:] + bz[1:-1,slice_index,:-1])
            return 0.5*(bx1**2 + by1**2+ bz1**2)

        if np.max(magfield(bx,by,bz)) > 1e-6:
            beta = 4*np.pi*pr[1:-1,slice_index,1:-1].T/magfield(bx,by,bz).T
        else:
            beta = 0.0*pr[1:-1,slice_index,1:-1].T

        current_test(bx, by, bz)

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

        by_reference_flux = np.max(np.abs(by[ny//2,slice_index,z_photo:-1]))
        bx_interior = bx[:,1:-1,1:-1]
        #Is it an arcade or not? Check by taking a 1D slice through the x component, with some smooothing to remove noise early on
        checkslice = gaussian_filter1d(bx_interior[nx//2,slice_index,z_photo:], 2)

        def categorise(checkslice):
            #Determine the topology based on this slice
            #Initially will be negative followed by positive
            #Check for rope
            rope_index = -1
            arcade = False
            rope = False
            rope_height = np.nan
            arcade_height = np.nan
            signs = np.sign(checkslice)
            flips = np.where(signs[1:]*signs[:-1] == -1)[0]   #where the magnetic field changes sign
            for flip in flips:
                if signs[flip] > 0.0 and np.max(checkslice[:flip] > by_reference_flux*0.25):   #positive to negative
                    rope_height = zc[z_photo:][flip]
                    rope = True
                else:
                    rope_height = np.nan
            #Check for arcade (before the rope forms and erupts)
            for flip in flips:
                if signs[flip] < 0.0 and np.min(checkslice[:flip]) < -by_reference_flux*0.25:
                    arcade = True
                    arcade_height = zc[z_photo:][flip]

            return arcade, arcade_height, rope, rope_height

        try:
            arcade, aheight, rope, rheight = categorise(checkslice)
        except:
            arcade, aheight, rope, rheight = False, np.nan, False, np.nan

        #Run through vertical slices of these and output average values.
        #Can do Fourier transform if this makes sense, but will do so less than the flux rope sims

        lf_slice_avgs = np.zeros((nz))
        for k in range(nz):
            lf_slice_avgs[k] = np.sum(np.abs(l2[1:-1,1:-1,k]))/(np.size(l2[1:-1,1:-1,k]))

        plt.plot(zc[1:-1], lf_slice_avgs[1:-1], label = data_directory, c = cs[source_number])
        ymax = max(ymax, np.max(lf_slice_avgs[z_photo:-1]))
        ymin = min(ymin, np.min(lf_slice_avgs[z_photo:-1]))

        plt.plot([rheight, rheight], [ymin,ymax], linestyle = 'dashed', c = cs[source_number])

    plt.plot([10, 10], [ymin,ymax], linestyle = 'dotted', c = 'black')
    plt.plot([20, 20], [ymin,ymax], linestyle = 'dotted', c = 'black')
    plt.ylim(0,ymax)
    plt.xlabel('Height')
    plt.ylabel('Absolute Lorentz Force')
    plt.title('Output number %d' % plot_num)
    plt.legend()
    #plt.show()
    plt.savefig('./lorentz/%04d.png' % plot_num)
    plt.close()























