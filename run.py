#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for running multiple instances of lare2d, varying certain parameters

"""


import os
import shutil
import numpy as np
import sys
from numpy import random
import time
from scipy.io import netcdf_file
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

if len(sys.argv) > 2:
    ncores = int(sys.argv[2])
else:
    ncores = 1

if os.uname()[1] == 'brillouin.dur.ac.uk':
    machine_flag = 0
elif os.uname()[1][-14:] == 'ham8.dur.ac.uk':
    machine_flag = 1
else:
    machine_flag = 2

if machine_flag == 0:
    data_directory = '/home/grads/trcn27/rdata/flux_emergence/%03d/' % run
if machine_flag == 1:
    data_directory = '/nobackup/trcn27/flux_emergence/%03d/' % run
if machine_flag == 2:
    data_directory = ('./Data%03d/' % run)

if os.path.isdir(data_directory):
    for i in range(1000):
        if os.path.isfile('%s%04d.nc' % (data_directory, i)):
            os.remove('%s%04d.nc' % (data_directory, i))
else:
    os.mkdir(data_directory)

if not os.path.isdir('./inits'):
    os.mkdir('./inits')

if not os.path.isdir('./parameters'):
    os.mkdir('./parameters')

if not os.path.isdir('./diagnostics'):
    os.mkdir('./diagnostics')

nx = 64
ny = 64
nz = 64

x0 = -130.; x1 = 130.0
y0 = -130.; y1 = 130.0
z0 = -25.0; z1 = 100.0

shearfact = 0.0
bfact = 5.0

density_init = 1.0
energy_init = 1.5e-2

nplots = 250
ndiags = 000
tmax = 1.0

eta = 1e-6

nu0_decay = 0.0

#factors = [0.0, 0.01, 0.1, 0.2, 0.5, 1.0, 2.0]
#energy_factor = factors[run]
energy_factor = 0.0#1e-3

zstar = 24.0/50.0
stratification_factor = 1.0

variables = np.zeros((30))

variables[0] = run
variables[1] = nx
variables[2] = tmax
variables[3] = nplots
variables[4] = ndiags
variables[5] = bfact
variables[6] = shearfact
variables[7] = x0
variables[8] = x1
variables[9] = y0
variables[10] = y1
variables[11] = z0
variables[12] = z1
variables[13] = nx
variables[14] = ny
variables[15] = nz
variables[16] = machine_flag
variables[17] = density_init
variables[18] = stratification_factor
variables[19] = energy_factor
variables[20] = zstar

if True:
    class Grid():
        def __init__(self):
            self.x0 = x0; self.x1 = x1
            self.y0 = y0; self.y1 = y1
            self.z0 = z0; self.z1 = z1
            self.nx = nx ; self.ny = ny; self.nz = nz
            self.dx = (self.x1 - self.x0)/nx
            self.dy = (self.y1 - self.y0)/ny
            self.dz = (self.z1 - self.z0)/nz
            self.xs = np.linspace(self.x0,self.x1,self.nx+1)
            self.ys = np.linspace(self.y0,self.y1,self.ny+1)
            self.zs = np.linspace(self.z0,self.z1,self.nz+1)

    def lbound_pariat(x,y):
        #Outputs lower boundary radial magnetic field as a function of position x
        sf = x1/12.0
        lbound_fn = np.zeros((len(x),len(y)))
        dipole_mag = 25.0; zstar = z1*1.5/24.0
        for i, xi in enumerate(x):
            for j,yj in enumerate(y):
                lbound_fn[i,j] = sf**3*dipole_mag*(2*(zstar)**2 - ((xi)**2 + (yj)**2))/(((xi)**2 + (yj)**2 + (zstar)**2)**2.5)
        print('Max lbound', np.max(lbound_fn))
        print('Lbound flux', np.sum(np.abs(lbound_fn)))

        return lbound_fn

    #Make initial jet condition
    grid = Grid()
    #compute_initial_condition(grid, lbound_pariat,run, background_strength = 1.0, boundary_error_limit = 1e-4, init_filename = './inits/init%03d.nc' % run)


os.system('make')

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
#compute_initial_condition(Grid(), lbound_fn, run)

if machine_flag < 0.5:
    os.system('/usr/lib64/openmpi/bin/mpiexec -np %d ./bin/lare3d %d' % (ncores, run))
