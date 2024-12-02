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

#matplotlib.rcParams['text.usetex'] = True
#Scripts for plotting various quantities with height -- perhaps for putting together into an animation?

if len(sys.argv) > 1:
    plot = int(sys.argv[1])
else:
    plot = 0


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

z_photo = int((nz)*(10.0 - zs[0])/(zs[-1] - zs[0]))   #Photosphere height

class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

height_min = int((nz)*(10.0 - zc[0])/(zc[-1] - zc[0]))
height_max = height_min + 1

#PLOT QUANTITIES AT VARIOUS SLICES IN HEIGHT
directories = ['../Data_15/', '../Data_150/']
titles = ['Less Stratified', 'More Stratified']
titles2 = ['B_x', 'B_y', 'B_z', 'Density']

height_min = 0
height_max = 128

if False:
    if not os.path.exists('./heightslice/'):
        os.mkdir('./heightslice')


    for hi in range(height_min, height_max):

        diffs = [[] for _ in range(4)]
        #hi = int((nz)*(height - zc[0])/(zs[-1] - zs[0]))
        fig, axs = plt.subplots(3,4, figsize = (10,7.5))

        for source in range(2):
            data_directory = directories[source]

            fname = '%s%04d.nc' % (data_directory, plot)

            try:
                data = netcdf_file(fname, 'r', mmap=False)
            except:
                continue

            bx = np.zeros((nx+1,ny+2,nz+2))
            by = np.zeros((nx+2,ny+1,nz+2))
            bz = np.zeros((nx+2,ny+2,nz+1))

            bx[:,1:-1,1:-1] = np.swapaxes(data.variables['bx'][:],0,2)
            by[1:-1,:,1:-1] = np.swapaxes(data.variables['by'][:],0,2)
            bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)

            en = np.zeros((nx,ny,nz))
            rho = np.zeros((nx,ny,nz))

            en = np.swapaxes(data.variables['en'][:],0,2)
            rho = np.swapaxes(data.variables['rho'][:],0,2)

            vx = np.zeros((nx+1,ny+1,nz+1))
            vy = np.zeros((nx+1,ny+1,nz+1))
            vz = np.zeros((nx+1,ny+1,nz+1))

            vx = np.swapaxes(data.variables['vx'][:],0,2)
            vy = np.swapaxes(data.variables['vy'][:],0,2)
            vz = np.swapaxes(data.variables['vz'][:],0,2)

            pr = rho*en*(2/3)

            data.close()

            def b_to_centres(bx, by, bz):
                #Averages the magnetic field to grid centres, allowing for easier plotting
                bx1 = 0.5*(bx[1:,1:-1,1:-1] + bx[:-1,1:-1,1:-1])
                by1 = 0.5*(by[1:-1,1:,1:-1] + by[1:-1,:-1,1:-1])
                bz1 = 0.5*(bz[1:-1,1:-1,1:] + bz[1:-1,1:-1,:-1])
                return bx1, by1, bz1

            bx1, by1, bz1 = b_to_centres(bx, by, bz)

            h_slice = hi#int((nz)*(10.0 - zs[0])/(zs[-1] - zs[0])) #slice index of the photosphere

            column = 0
            ax = axs[source,column]
            toplot = bx1[:,:,h_slice]
            im = ax.pcolormesh(xc,yc,toplot.T,vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)), cmap ='seismic')
            ax.set_title('%s \n %s' % (titles[source], titles2[column]), fontsize = 10)
            plt.colorbar(im, ax=ax)

            if source == 0:
                diffs[column] = toplot
            else:
                diffs[column] = diffs[column] - toplot

            column = 1
            ax = axs[source,column]
            toplot = by1[:,:,h_slice]
            im = ax.pcolormesh(xc,yc,toplot.T,vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)), cmap ='seismic')
            ax.set_title('%s \n %s' % (titles[source], titles2[column]), fontsize = 10)
            plt.colorbar(im, ax=ax)

            if source == 0:
                diffs[column] = toplot
            else:
                diffs[column] = diffs[column] - toplot

            column = 2
            ax = axs[source,column]
            toplot = bz1[:,:,h_slice]
            im = ax.pcolormesh(xc,yc,toplot.T,vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)), cmap ='seismic')
            ax.set_title('%s \n %s' % (titles[source], titles2[column]), fontsize = 10)
            plt.colorbar(im, ax=ax)

            if source == 0:
                diffs[column] = toplot
            else:
                diffs[column] = diffs[column] - toplot

            column = 3
            ax = axs[source,column]
            toplot = rho[:,:,h_slice]
            im = ax.pcolormesh(xc,yc,toplot.T,vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)), cmap ='seismic')
            ax.set_title('%s \n %s' % (titles[source], titles2[column]), fontsize = 10)
            plt.colorbar(im, ax=ax)

            if source == 0:
                diffs[column] = toplot
            else:
                diffs[column] = diffs[column] - toplot

        for column in range(4):
            ax = axs[2,column]
            toplot = diffs[column]
            im = ax.pcolormesh(xc,yc,toplot.T,vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)), cmap ='seismic')
            ax.set_title('%s \n %s' % ('Difference', titles2[column]), fontsize = 10)
            plt.colorbar(im, ax=ax)

            if source == 0:
                diffs[column] = toplot
            else:
                diffs[column] = diffs[column] - toplot

        plt.suptitle('Height = %.1f' % zc[hi])
        plt.tight_layout()
        #plt.show()
        plt.savefig('heightslice/%04d' % hi)
        plt.close()

if False:
    fig = plt.figure(figsize = (10,7))
    densities = []
    for source in range(2):
        densities.append([])
        hs = []
        data_directory = directories[source]

        fname = '%s%04d.nc' % (data_directory, plot)

        try:
            data = netcdf_file(fname, 'r', mmap=False)
        except:
            continue

        rho = np.swapaxes(data.variables['rho'][:],0,2)

        for hi in range(height_min, height_max):
            densities[source].append(np.sum(rho[:,:,hi])*dx*dy)
            hs.append(zc[hi])

    densities = np.array(densities)
    for source in range(2):
        plt.plot(hs, densities[source,:], label = titles[source])
    plt.yscale('log')
    plt.legend()
    plt.title('Densities with height')
    plt.tight_layout()
    plt.savefig('Densities')
    plt.show()

if True:
    fig = plt.figure(figsize = (10,7))
    energies = []
    for source in range(2):
        energies.append([])
        hs = []
        data_directory = directories[source]

        fname = '%s%04d.nc' % (data_directory, plot)

        try:
            data = netcdf_file(fname, 'r', mmap=False)
        except:
            continue

        en = np.swapaxes(data.variables['en'][:],0,2)
        rho = np.swapaxes(data.variables['rho'][:],0,2)

        pr = rho*en*(2/3)

        for hi in range(height_min, height_max):
            energies[source].append(np.sum(en[:,:,hi])*dx*dy)
            hs.append(zc[hi])

    energies = np.array(energies)
    for source in range(2):
        plt.plot(hs, energies[source,:], label = titles[source])
    plt.yscale('log')
    plt.legend()
    plt.title('Energies with height')
    plt.tight_layout()
    plt.savefig('../energies/%04d.png' % plot)
    plt.close()








