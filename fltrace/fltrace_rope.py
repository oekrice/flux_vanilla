#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for fortran field line tracer

Uses FLH density on the photosphere to determine which lines to trace.

Needs routine from emergence_comparisons to do that. Bit of a mess...

Plots field lines based on the height of the rope -- gaussian distribution from there? Gah...
"""


import os
import shutil
import numpy as np
import sys
from numpy import random
import time
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
import pyvista as pv
from get_flh import FLH
from scipy.ndimage import gaussian_filter1d
#pv.start_xvfb()

class trace_fieldlines():
    def __init__(self, snap_min, snap_max):
        #Establish random seeds for field line plotting
        if not os.path.isfile('start_seeds.txt'):
            self.start_seeds = random.rand(1000**2)
            np.savetxt('start_seeds.txt', self.start_seeds, delimiter = ',')
        else:
            self.start_seeds = np.loadtxt('start_seeds.txt', delimiter = ',')

        #Establish grid parameters (can be read in from elsewhere of course)
        for snap_number in range(snap_min, snap_max):
            self.run = snap_number
            self.snap = snap_number
            self.print_flag = 1

            self.save_number = self.snap
            self.option = 3   #tracing a plotting options (1 for jet, 2 for emergence)

            self.data_source = 0.
            self.data_root = '../Data_15/'
            data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
            self.bx = np.swapaxes(data.variables['bx'][:],0,2)
            self.by = np.swapaxes(data.variables['by'][:],0,2)
            self.bz = np.swapaxes(data.variables['bz'][:],0,2)
            data.close()

            #Establish start points for the field line plotting
            self.max_line_length = 10000
            self.ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
            self.weakness_limit = 5e-3   #Minimum field strength to stop plotting
            self.line_plot_length = 100  #To save time while plotting, reduce the length of the plotted lines
            self.nlines = 500

            #Import bz as a test of the resolutions (and for the pyvista plot)
            self.nx = np.shape(self.bz)[0]
            self.ny = np.shape(self.bz)[1]
            self.nz = np.shape(self.bz)[2] - 1

            self.x0 = -130.; self.x1 = 130.
            self.y0 = -130.; self.y1 = 130.
            self.z0 = -25.; self.z1 = 100.

            self.xs = np.linspace(self.x0,self.x1,self.nx+1)
            self.ys = np.linspace(self.y0,self.y1,self.ny+1)
            self.zs = np.linspace(self.z0,self.z1,self.nz+1)

            self.xc = np.zeros(self.nx + 2)
            self.yc = np.zeros(self.ny + 2)
            self.zc = np.zeros(self.nz + 2)

            self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
            self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
            self.zc[1:-1] = 0.5*(self.zs[1:] + self.zs[:-1])

            self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
            self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])
            self.zc[0] = self.zc[0] - (self.zc[2] - self.zc[2])

            self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
            self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])
            self.zc[-1] = self.zc[-2] + (self.zc[-2] - self.zc[-3])

            self.dx = self.xs[1] - self.xs[0]
            self.dy = self.ys[1] - self.ys[0]
            self.dz = self.zs[1] - self.zs[0]

            #Folder admin
            if not os.path.exists('./fl_data/'):
                os.mkdir('fl_data')

            if not os.path.exists('./plots/'):
                os.mkdir('plots')

            self.plot_difference()

            os.system('rm ./fl_data/flines%03d.nc' % self.snap)
            os.system('rm ./fl_data/flparameters%03d.txt' % self.snap)
            os.system('rm ./fl_data/starts%03d.txt' % self.snap)
            #os.system('rm ./fl_data/plot_base%d.npy' % self.snap)
            #os.system('rm ./fl_data/photo1%d.npy' % self.snap)
            #os.system('rm ./fl_data/photo2%d.npy' % self.snap)
            os.system('rm ./fl_data/flines15_%03d.nc' % self.snap)
            os.system('rm ./fl_data/flines150_%03d.nc' % self.snap)



    def plot_difference(self):

        #Plots the two fields in subplots, hopefully. Need to do compiling and stoof in here.
        x, y = np.meshgrid(self.xs, self.ys)
        z = 10*np.ones((np.shape(x)))
        surface = pv.StructuredGrid(x, y, z)

        def do_subplot():


            data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
            self.bx = np.swapaxes(data.variables['bx'][:],0,2)
            self.by = np.swapaxes(data.variables['by'][:],0,2)
            self.bz = np.swapaxes(data.variables['bz'][:],0,2)
            data.close()

            self.set_starts()

            self.setup_tracer()
            #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
            self.trace_lines_fortran()

            print('Plotting in pyvista')
            for li, line in enumerate(self.lines):
                line_length = 0
                for k in range(len(line)):
                    if line[k,2] < 1e6 and line[k,2] >= 10.0:
                        line_length += 1
                    else:
                        break
                line = line[:line_length,:]
                #Thin out the lines (if required)
                if line_length > 1:
                    thin_fact = max(int(line_length/self.line_plot_length), 1)
                    thinned_line = line[:line_length:thin_fact].copy()
                    thinned_line[-1] = line[line_length-1].copy()
                else:
                    continue

                line = np.array(thinned_line).tolist()
                doplot = True

                if doplot:
                    p.add_mesh(pv.Spline(line, len(line)),color='white',line_width=0.1)

            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            p.camera.position = (200.0,100,150.0)
            p.camera.focal_point = (0,0,50)

        p = pv.Plotter(off_screen=True, shape = (2,1))
        p.background_color = "black"

        p.subplot(0, 0)
        self.data_source = 15.
        self.data_root = '../Data_15/'
        do_subplot()
        p.add_mesh(surface, scalars = np.abs(self.plot_base), show_edges=False,cmap = 'inferno', clim = [-0.25*np.max(self.plot_base), np.max(self.plot_base)])
        #p.add_mesh(surface, scalars = np.abs(self.flh_photo1), show_edges=False,cmap = 'inferno', clim = [-0.25*np.max(self.flh_photo1), np.max(self.flh_photo1)])
        p.remove_scalar_bar()

        p.add_title('Snap number %d' % self.snap, color = 'Red', font_size = 25)

        p.subplot(1, 0)
        self.data_source = 150.
        self.data_root = '../Data_150/'


        do_subplot()

        p.add_mesh(surface, scalars = np.abs(self.plot_base), show_edges=False,cmap = 'inferno', clim = [-0.25*np.max(self.plot_base), np.max(self.plot_base)])
        #p.add_mesh(surface, scalars = np.abs(self.flh_photo2), show_edges=False,cmap = 'inferno', clim = [-0.25*np.max(self.flh_photo1), np.max(self.flh_photo1)])
        p.remove_scalar_bar()

        p.show(screenshot='plots/b%04d.png' % self.save_number, window_size = (2000,2000))
        p.close()

    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.

        #Plot up from the surface, based on some threshold of how strong the magnetic field is... Could be fun.
        z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

        #Set rope centre in x, z plane

        z_slice = self.bx[self.nx//2,self.nz//2,:]

        z_slice = gaussian_filter1d(z_slice, 2)

        checkslice = z_slice[z_photo:]

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

            for i in range(len(flips)):
                flip = flips[i]
                if signs[flip] > 0.0:   #positive to negative
                    rope_height = self.zc[z_photo:][flip]
                    rope = True
                else:
                    rope_height = np.nan

            #Check for arcade (before the rope forms and erupts)
            for flip in flips:
                if signs[flip] > 0.0 and np.min(checkslice[:flip]) < -by_reference_flux*0.25:
                    arcade = True
                    arcade_height = self.zc[z_photo:][flip]

            return arcade, arcade_height, rope, rope_height

        try:
            arcade, aheight, rope, rheight = categorise(checkslice)
        except:
            arcade, aheight, rope, rheight = False, np.nan, False, np.nan

        print('Rope height', rheight)

        rope_centre = 0.0
        for k in range(self.nz-1, 0, -1):
            if z_slice[k] > 0.0:
                rope_centre = self.zs[k]
                break

        print('R centre', rope_centre)
        centre_coordinates = [0.0, 0.0, rope_centre]

        gaussian_array = np.zeros((len(self.xc), len(self.zc)))
        alpha = 200.0
        for i in range(len(self.xc)):
            for j in range(len(self.zc)):
                gaussian_array[i,j] = np.exp((-(self.xc[i]-centre_coordinates[0])**2 - (self.zc[j] - centre_coordinates[2])**2)/alpha)

        #surface_array = self.bz[:,:,z_photo]   #Distribution of surface magnetic field. This array is all that needs changing
        surface_array = gaussian_array   #Distribution of surface FLH
        #Want to plot the lines based on where the flh differences are highest

        self.plot_base = surface_array
        max_surface = np.max(np.abs(surface_array)) + 1e-6

        nlines = self.nlines

        alpha = 2.0
        alphasum = np.sum(np.abs(surface_array)**alpha)
        pb = max_surface**alpha*nlines/alphasum

        print('prob', pb, nlines, alphasum)

        self.starts = []

        #Add capcbility for higher-resolution base
        xscale = np.shape(self.plot_base)[0]/self.nx
        zscale = np.shape(self.plot_base)[1]/self.nz

        nxl = np.shape(self.plot_base)[0]
        nzl = np.shape(self.plot_base)[1]

        self.xsl = np.linspace(self.xs[0], self.xs[-1], nxl)
        self.zsl = np.linspace(self.zs[0], self.zs[-1], nzl)

        xcl = 0.5*(self.xsl[1:] + self.xsl[:-1])
        zcl = 0.5*(self.zsl[1:] + self.zsl[:-1])


        cellcount = 0
        for i in range(nxl-1):  #run through lower surface cells
            for j in range(nzl-1):
                prop = np.abs(surface_array[i,j])/max_surface
                #if prop > 0.9:
                if self.start_seeds[cellcount] < pb*prop**alpha:
                    self.starts.append([xcl[i],0.0,zcl[j]])
                cellcount += 1

        print('Tracing', len(self.starts), 'lines')

        self.nstarts = len(self.starts)
        self.starts = np.array(self.starts).reshape(self.nstarts*3)

        self.plot_base = self.bz[:,:,z_photo]

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        variables = np.zeros((30))

        variables[0] = self.run
        variables[1] = self.nx
        variables[2] = self.ny
        variables[3] = self.nz
        variables[4] = self.x0
        variables[5] = self.x1
        variables[6] = self.y0
        variables[7] = self.y1
        variables[8] = self.z0
        variables[9] = self.z1
        variables[10] = self.snap
        variables[11] = self.nstarts
        variables[12] = self.print_flag
        variables[13] = self.max_line_length
        variables[14] = self.ds
        variables[15] = self.weakness_limit
        variables[16] = self.data_source

        np.savetxt('./fl_data/flparameters%03d.txt' % self.snap, variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('./fl_data/starts%03d.txt' % self.snap, self.starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):


        if False:#os.path.isfile('./fl_data/flines%d_%03d.nc' % (self.data_source, self.snap)):
            data = netcdf_file('./fl_data/flines%d_%03d.nc' % (self.data_source, self.snap), 'r', mmap=False)
            #Import existing field lines
        else:
            os.system('make')
            if os.uname()[1] == 'brillouin.dur.ac.uk':
                os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace %d' % self.snap)
            elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
                os.system('mpiexec -np 1 ./bin/fltrace %d' % self.snap)
            else:
                os.system('mpirun -np 1 ./bin/fltrace %d' % self.snap)
            os.system('scp -r ./fl_data/flines%03d.nc ./fl_data/flines%d_%03d.nc' % (self.snap, self.data_source, self.snap))

            data = netcdf_file('./fl_data/flines%d_%03d.nc' % (self.data_source, self.snap), 'r', mmap=False)
            try:
                data = netcdf_file('./fl_data/flines%03d.nc' % self.snap, 'r', mmap=False)
                print('Field lines found')

            except:
                print('File not found')

        self.lines = np.swapaxes(data.variables['lines'][:],0,2)

nset = 1000 #Number of concurrent runs. Receives input 0-(nset-1)
set_num = int(sys.argv[1])
snap_min = 0 + set_num
while snap_min < 400:
    print('Run number', snap_min)
    trace_fieldlines(snap_min = snap_min, snap_max = snap_min+1)
    snap_min = snap_min + nset
