#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for fortran field line tracer

Uses FLH density on the photosphere to determine which lines to trace.

Needs routine from emergence_comparisons to do that. Bit of a mess...
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

            self.data_root = '../Data_15/'
            data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
            self.bx = np.swapaxes(data.variables['bx'][:],0,2)
            self.by = np.swapaxes(data.variables['by'][:],0,2)
            self.bz = np.swapaxes(data.variables['bz'][:],0,2)
            data.close()

            #Establish start points for the field line plotting
            self.max_line_length = 10000
            self.ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
            self.weakness_limit = 1e-2   #Minimum field strength to stop plotting
            self.line_plot_length = 100  #To save time while plotting, reduce the length of the plotted lines

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

            if False:   #Option to plot based on the FLH density along
                flh = FLH(self)    #Do field-line helicity things
                self.flh_photo = flh.flh_photo #FLH density on the photosphere
                self.plot_base = self.flh_photo
            else:   #Option to plot based on the differences in FLH density
                self.data_root = '../Data_15/'
                data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
                self.bx = np.swapaxes(data.variables['bx'][:],0,2)
                self.by = np.swapaxes(data.variables['by'][:],0,2)
                self.bz = np.swapaxes(data.variables['bz'][:],0,2)
                data.close()

                flh1 = FLH(self)    #Do field-line helicity things
                self.flh_photo1 = flh1.flh_photo #FLH density on the photosphere

                self.data_root = '../Data_150/'
                data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
                self.bx = np.swapaxes(data.variables['bx'][:],0,2)
                self.by = np.swapaxes(data.variables['by'][:],0,2)
                self.bz = np.swapaxes(data.variables['bz'][:],0,2)
                data.close()

                flh2 = FLH(self)    #Do field-line helicity things
                self.flh_photo2 = flh2.flh_photo #FLH density on the photosphere

                self.plot_base = self.flh_photo1 - self.flh_photo2 #FLH difference
            #self.flh_photo = np.ones((self.nx, self.ny))
            #Find start points
            self.set_starts()
            #Create runtime variables for fortran
            self.setup_tracer()
            #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
            self.trace_lines_fortran()
            #Plot the field lines (using pyvista)
            if True:
                if not os.path.exists('./plots/'):
                    os.mkdir('plots')
                self.plot_vista()

            os.system('rm ./fl_data/flines%03d.nc' % self.snap)
            os.system('rm ./fl_data/flparameters%03d.txt' % self.snap)
            os.system('rm ./fl_data/starts%03d.txt' % self.snap)

    def plot_vista(self):
        print('Plotting...')
        x, y = np.meshgrid(self.xs, self.ys)
        z = 10*np.ones((np.shape(x)))
        surface = pv.StructuredGrid(x, y, z)
        p = pv.Plotter(off_screen=False)
        p.background_color = "black"

        for li, line in enumerate(self.lines):
            line = np.array(line)
            line_length = len(line[line[:,2]<1e6])
            #Thin out the lines (if required)
            if line_length > 0:
                thin_fact = max(int(line_length/self.line_plot_length), 1)
                thinned_line = line[:line_length:thin_fact].copy()
                thinned_line[-1] = line[line_length-1].copy()
            else:
                continue

            line = np.array(thinned_line).tolist()
            doplot = True
            if line_length == 0:
                doplot = False

            if line[-1][2] > 11.0 and line[-1][2] < 99.0:
                doplot = False
            if doplot:
                p.add_mesh(pv.Spline(line, len(line)),color='white',line_width=0.25)

        if self.option == 1:
            p.add_mesh(surface, scalars= self.bz[:,:,0], show_edges=True,cmap = 'plasma')
            p.camera.position = (20.0,40,20.0)
            p.camera.focal_point = (0,0,4)
        if self.option > 1:
            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))
            p.add_mesh(surface, scalars= self.plot_base, show_edges=False,cmap = 'plasma')
            p.camera.position = (400.0,200,250.0)
            p.camera.focal_point = (0,0,0)


        #p.show(screenshot='plots/b%04d.png' % self.save_number, window_size = (1000,1000))
        #print('Plot saved to file plots/b%04d.png' % self.save_number)
        p.show(screenshot='difftestb%04d.png' % self.save_number, window_size = (1000,1000))


    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.
        if self.option == 1:
            self.starts = []
            #Trace from the top
            nrs = 10; nthetas = 2
            ris = np.linspace(3.0/nrs,self.x0-3.0/nrs,nrs)
            tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
            for i in range(nrs):
                for j in range(nthetas):
                    self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),self.z1-1e-6])
            #And from the bottom (for the interior ones only)
            nrs = 20; nthetas = 10
            ris = np.linspace(0.5*self.x0/nrs,self.x0-0.5*self.x0/nrs,nrs)
            tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
            for i in range(nrs):
                for j in range(nthetas):
                    self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),1e-6])

        if self.option  == 2:
            #Plot outwards from the symmetry line
            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            self.starts = []

            nxs = 80; nys = 80; nzs = 80

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)
            zks = np.linspace(10.0,self.z1-1e-6,nzs)

            for i in range(nxs):
                j = nys//2
                for k in range(nzs):

                    xp = int((self.nx)*(xis[i] - self.x0)/(self.x1 - self.x0))
                    yp = int((self.ny)*(yjs[j] - self.y0)/(self.y1 - self.y0))
                    zp = int((self.nz)*(zks[k] - self.z0)/(self.z1 - self.z0))

                    if np.abs(self.by[xp,yp,zp]) > 0.01:
                        self.starts.append([xis[i],yjs[j],zks[k]])
            print('Tracing', len(self.starts), 'lines')

        if self.option  == 3:
            #Plot up from the surface, based on some threshold of how strong the magnetic field is... Could be fun.
            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))


            #surface_array = self.bz[:,:,z_photo]   #Distribution of surface magnetic field. This array is all that needs changing
            surface_array = self.plot_base   #Distribution of surface FLH
            #Want to plot the lines based on where the flh differences are highest


            max_surface = np.max(np.abs(surface_array)) + 1e-6

            nlines = 1000

            alpha = 2.0
            alphasum = np.sum(np.abs(surface_array)**alpha)
            pb = max_surface**alpha*nlines/alphasum

            print('prob', pb, nlines, alphasum)

            self.starts = []

            cellcount = 0
            for i in range(self.nx):  #run through grid cells
                for j in range(self.ny):
                    prop = np.abs(surface_array[i,j])/max_surface
                    if self.start_seeds[cellcount] < pb*prop**alpha:
                        self.starts.append([self.xc[i+1],self.yc[j+1],10.0])
                    cellcount += 1

            print('Tracing', len(self.starts), 'lines')


        self.nstarts = len(self.starts)
        self.starts = np.array(self.starts).reshape(self.nstarts*3)

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.05 #Tracing 'timestep' as a proportion of the grid size
        weakness_limit = 1e-3   #Minimum field strength to stop plotting
        print_flag = 1  #Print some things as the tracing happens

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

        np.savetxt('./fl_data/flparameters%03d.txt' % self.snap, variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('./fl_data/starts%03d.txt' % self.snap, self.starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):
        os.system('make')
        if os.uname()[1] == 'brillouin.dur.ac.uk':
            os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
            os.system('mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        else:
            os.system('mpirun -np 1 ./bin/fltrace %d' % self.snap)

        try:
            data = netcdf_file('./fl_data/flines%03d.nc' % self.snap, 'r', mmap=False)
            print('Field lines found')

        except:
            print('File not found')

        self.lines = np.swapaxes(data.variables['lines'][:],0,2)


if len(sys.argv) > 1:
    snap_min = int(sys.argv[1])
else:
    snap_min = 0

trace_fieldlines(snap_min = snap_min, snap_max = snap_min+1)
