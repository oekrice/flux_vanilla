#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np 
from matplotlib import pyplot as plt
import fieldLineTopology as flt
from streamtracer import StreamTracer, VectorGrid
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simps
#import waveletRoutines as wr

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import netcdf_file
import os
import sys
# Script called from the field line tracer to get a distribution of the FLH on the lower surface.

class FLH():
    def __init__(self, grid):

        high_resolution = True
        #Get grid and stuff in here so don't have to do it twice
        bxOg = grid.bx
        byOg = grid.by
        bzOg = grid.bz

        bx = 0.5*(bxOg[1:,:,:] + bxOg[:-1,:,:])
        by = 0.5*(byOg[:,1:,:] + byOg[:,:-1,:])
        bz = 0.5*(bzOg[:,:,1:] + bzOg[:,:,:-1])

        print('Net flux', np.sum(bz[:,:,0]))
        bz = bz - np.sum(bz[:,:,0])/np.size(bz[:,:,0])
        print('Net flux', np.sum(bz[:,:,0]))

        xv = grid.xc[1:-1]
        yv = grid.yc[1:-1]
        zv = grid.zc[1:-1]

        X, Y = np.meshgrid(xv, yv, indexing='ij')
        # Flatten the grid arrays to form the input to the interpolator
        points = np.vstack([X.ravel(), Y.ravel()]).T
        dx = xv[1]-xv[0]
        dy = yv[1]-yv[0]
        dz = zv[1]-zv[0]
        dA = dx*dy
        grid_spacing = [dx,dy,dz]
        grid_ncells = [grid.nx, grid.ny, grid.nz]

        bField = flt.createSingleField(bx,by,bz)

        z_photo = int((grid.nz)*(10.0 - zv[0])/(zv[-1] - zv[0]))

        AField = flt.getAFastSingle(bField,grid_ncells,grid_spacing)

        bField_test = flt.curl(AField,grid_spacing)

        usf = flt.unitSpeedField(bField.copy(),0.01)  #transforms bz to be 'unit speed' in the z direction
        BUnit = flt.addDivergenceCleaningTerm(usf,grid_ncells,grid_spacing)   #Returns unit speed field.
        AWind = flt.getAFastSingle(BUnit,grid_ncells,grid_spacing)

        curlField = flt.curl(bField,grid_spacing)

        bzConst = np.sum(bz[:,:,0])/(dA*(grid.nx-1)*(grid.ny-1))
        AConst = flt.AConst(bzConst,points,dA,[grid.nx, grid.ny, grid.nz])
        AField = AField + AConst

        domain_size = 20
        print('Determining FLH')
        grid_spacing = [domain_size/grid.nx,domain_size/grid.ny,domain_size/grid.nz]   #dx, dy, dz

        #set domain lengths

        lx = grid_spacing[0]*grid.nx
        ly = grid_spacing[1]*grid.ny
        lz = grid_spacing[2]*grid.nz

        #set number of points used to calculate the distribution from
        if high_resolution:
            nxl = 2*grid.nx
            nyl = 2*grid.ny
        else:
            nxl = grid.nx
            nyl = grid.ny

        # set a minimum strength of field line cut off
        bCut = 0.01

        testFLHDen = flt.getFLHDenSingle(bField,AField)

        # set z value for anchoring plane (use photosphere z index = 116 here)
        z_photo = 20.0*(10.0 - zv[0])/(zv[-1] - zv[0])

        BxInterp,ByInterp,BzInterp = flt.getInterpolatedFieldSingle(bField,grid_spacing[0],grid_spacing[1],grid_spacing[2])

        #This interpolates onto a grid with coordinates of bx. Running from 0 to nx-1 etc.
        #Is a function


        #Shift everything so it matches up more nicely after the interpolation
        xscale = domain_size/(grid.x1- grid.x0)
        yscale = domain_size/(grid.y1- grid.y0)
        zscale = domain_size/(grid.z1- grid.z0)

        if high_resolution:
            fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[0+grid.dx*xscale/2,20-grid.dx*xscale/2],[0+grid.dy*yscale/2,20-grid.dy*yscale/2],z_photo,nxl,nyl,bCut)   #number ranges are the bounds to do the plotting (0-20)
        else:
            fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[0+grid.dx*xscale,20-grid.dx*xscale],[0+grid.dy*yscale,20-grid.dy*yscale],z_photo,nxl,nyl,bCut)   #number ranges are the bounds to do the plotting (0-20)

        flhInterp = flt.getInterpolatedQuantity(testFLHDen,grid_spacing)
        #flhBzInterp = flt.getInterpolatedQuantity(testFLHDen*bz[:,:,z_photo],grid_spacing)
        #flwindInterp = flt.getInterpolatedQuantity(testWindDen,grid_spacing)
        #twistInterp = flt.getInterpolatedQuantity(twistDensity,grid_spacing)
        # In[ ]:

        flh,indexesflh = flt.fieldLineIntegratedQuantity(flhInterp,goodSeeds,fieldLinesList,seeds,nxl,nyl)
        #flhBz,indexesflh = flt.fieldLineIntegratedQuantity(flhBzInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
        #flw,indexesflw = flt.fieldLineIntegratedQuantity(flwindInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
        #twistF,indexestwist = flt.fieldLineIntegratedQuantity(twistInterp,goodSeeds,fieldLinesList,seeds,nx,ny)

        self.flh_photo = flh
