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
# Code to compare field line quanitites (twist, helicity) etc. between the stratified and unstratified jet simulations

if len(sys.argv) > 1:
    snap_id = int(sys.argv[1])
else:
    snap_id = 300

flhs = []; flws = []; flhBzs = []; twistFs = []

# In[65]:
if False:
    if not os.path.exists('./Data_stratified/'):
        os.mkdir('./Data_stratified/')
    if not os.path.exists('./Data_unstratified/'):
        os.mkdir('./Data_unstratified/')

    os.system('scp -r pjb20205@login.hpc.strath.ac.uk:/users/pjb20205/lare3d_jet/Data_unstratified_192/%04d.nc ./Data_unstratified_192/' % snap_id)
    os.system('scp -r pjb20205@login.hpc.strath.ac.uk:/users/pjb20205/lare3d_jet/Data_verystratified_192/%04d.nc ./Data_verystratified_192/' % snap_id)

# read in file

paths = ['./Data_15/', './Data_150/']
titles = ['Strat ratio 15', 'Strat ratio 150']

nstrats = 2

for strat_flag in range(nstrats):   #do unstratified (top) and stratified (bottom)
    print('Dealing with field number', strat_flag)
    print('_______________________________________')
    path = paths[strat_flag]
    file2read = netcdf_file(path+'%04d.nc'% snap_id,'r')

    # In[67]:

    bxOg = np.swapaxes(file2read.variables['bx'][:], 0, 2).copy()
    byOg = np.swapaxes(file2read.variables['by'][:], 0, 2).copy()
    bzOg = np.swapaxes(file2read.variables['bz'][:], 0, 2).copy()

    file2read.close()
    # In[69]:

    bx = 0.5*(bxOg[1:,:,:] + bxOg[:-1,:,:])
    by = 0.5*(byOg[:,1:,:] + byOg[:,:-1,:])
    bz = 0.5*(bzOg[:,:,1:] + bzOg[:,:,:-1])

    ncells = 128  #Cut a bit off to test if I'm right about the scaling. I'm not...
    prop = 0
    if prop > 0:
        bx = bx[prop:-prop,prop:-prop,prop:-prop]
        by = by[prop:-prop,prop:-prop,prop:-prop]
        bz = bz[prop:-prop,prop:-prop,prop:-prop]

    print('Net flux', np.sum(bz[:,:,0]))
    bz = bz - np.sum(bz[:,:,0])/np.size(bz[:,:,0])
    print('Net flux', np.sum(bz[:,:,0]))

    ncells = int(np.size(bx)**(1/3)) + 1
    print(ncells)
    #Reshape so its 0= x comp 1 =y comp 2 = z comp

    xv = np.linspace(-130,130,ncells)
    yv = np.linspace(-130,130,ncells)
    zv = np.linspace(-25,100,ncells)

    X, Y = np.meshgrid(xv, yv, indexing='ij')
    # Flatten the grid arrays to form the input to the interpolator
    points = np.vstack([X.ravel(), Y.ravel()]).T
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]
    dz = zv[1]-zv[0]
    dA = dx*dy
    grid_spacing = [dx,dy,dz]
    grid_ncells = [ncells,ncells,ncells]


    # In[70]:

    bField = flt.createSingleField(bx,by,bz)
    #bFieldTest = flt.createSingleField(bxRot,byRot,bzRot)
    # In[71]:
    #get curl

    #calculate the winding gauge. Figure out if this actually works...
    #z component seems a bit wrong
    z_photo = int((ncells)*(10.0 - zv[0])/(zv[-1] - zv[0]))

    AField = flt.getAFastSingle(bField,grid_ncells,grid_spacing)

    bField_test = flt.curl(AField,grid_spacing)

    #calculate the winding gauge for the unit speed field
    usf = flt.unitSpeedField(bField.copy(),0.01)  #transforms bz to be 'unit speed' in the z direction
    BUnit = flt.addDivergenceCleaningTerm(usf,grid_ncells,grid_spacing)   #Returns unit speed field. Makes some sense... But why has BField changed?
    AWind = flt.getAFastSingle(BUnit,grid_ncells,grid_spacing)

    curlField= flt.curl(bField,grid_spacing)

    #flt.testfourier(bField)

    '''
    for i in range(3):
        im = axes[2,i].pcolormesh(BUnit[:,:,z_photo, i].T)
        fig.colorbar(im, ax=axes[2,i])
        axes[2,i].set_title('Unit Speed Field')
        im = axes[3,i].pcolormesh(curlField[:,:,z_photo, i].T)
        fig.colorbar(im, ax=axes[3,i])
        axes[3,i].set_title('Curlfield')

    fig.tight_layout()
    plt.savefig('extra_plots/inputs.png')
    plt.show()
    plt.close()
    '''

    # add constant component of A as there be net flux
    bzConst = np.sum(bz[:,:,0])/(dA*(ncells-1)*(ncells-1))
    AConst = flt.AConst(bzConst,points,dA,[ncells, ncells, ncells])
    AField = AField + AConst

    if False:
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 7))
        for i in range(3):
            im = axes[0,i].pcolormesh(bField[:,:,z_photo, i].T)
            fig.colorbar(im, ax=axes[0,i])
            axes[0,i].set_title('Original BField')
            im = axes[1,i].pcolormesh(bField_test[:,:,z_photo, i].T)
            fig.colorbar(im, ax=axes[1,i])
            axes[1,i].set_title('Curl of A')
            im = axes[2,i].pcolormesh(bField[:,:,z_photo, i].T - bField_test[:,:,z_photo, i].T)
            fig.colorbar(im, ax=axes[2,i])
            axes[2,i].set_title('Difference')


        plt.tight_layout()
        plt.savefig('comp.png')
        plt.show()

        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 7))

        for i in range(3):
            vmin = 1e6; vmax =  -1e6
            comp = bField[:,:,z_photo, i].T/bField_test[:,:,z_photo, i].T
            vmin = min(vmin, np.percentile(comp, 10))
            vmax = max(vmin, np.percentile(comp, 90))

            im = axes[i].pcolormesh(comp[:,:].T,vmin=vmin,vmax=vmax)
            fig.colorbar(im, ax = axes[i])
        plt.tight_layout()
        plt.show()
    # In[72]:

    # calculate the field line helcity and winding densities

    testFLHDen = flt.getFLHDenSingle(bField,AField)

    testWindDen = flt.getFLHDenSingle(bField,AWind)   #IS THIS RIGHT?
    # twist density

    twistDensity = flt.twistDen(bField,curlField,1e-6)

    # In[73]:

    z_photo = int((ncells)*(10.0 - zv[0])/(zv[-1] - zv[0]))
    #Check the field comps

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(bField[:,:,z_photo,i], origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.suptitle('Magnetic field')
    plt.savefig('extra_plots/0_%d.png' % strat_flag)

    plt.close()

    # In[1]:

    #Check the curl comps

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(twistDensity[:,:,z_photo*(i+1)], origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.suptitle('Twist Density')
    plt.savefig('extra_plots/1_%d.png' % strat_flag)
    plt.close()


    # In[61]:

    # FLH density, pretty !

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(testFLHDen[:,:,z_photo*(i+1)].T, origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.suptitle('FLH Density')
    plt.savefig('extra_plots/2_%d.png' % strat_flag)

    plt.close()


    # In[62]:

    # FLwind density, pretty !

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(testWindDen[:,:,z_photo*(i+1)].T, origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.suptitle('Winding Density')

    plt.savefig('extra_plots/3_%d.png' % strat_flag)

    plt.close()


    # Set up the field line tracer

    # In[63]:
    if True:
        grid_spacing = [20/ncells,20/ncells,20/ncells]   #dx, dy, dz

        #set domain lengths

        lx = grid_spacing[0]*ncells
        ly = grid_spacing[1]*ncells
        lz = grid_spacing[2]*ncells

        #set number of points used to calculate the distribution from
        nx = 200
        ny = 200

        # set a minimum strength of field line cut off
        bCut = 0.01

        # set z value for anchoring plane (use photosphere z index = 116 here)
        zv = 20.0*(10.0 - zv[0])/(zv[-1] - zv[0])

        print('Field read-in and tested')
        # In[ ]:

        print('Interpolating Field')
        BxInterp,ByInterp,BzInterp = flt.getInterpolatedFieldSingle(bField,grid_spacing[0],grid_spacing[1],grid_spacing[2])

        #This interpolates onto a grid with coordinates of bx. Running from 0 to nx-1 etc.
        #Is a function

        # In[ ]:

        print('Tracing field lines')

        fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[1,19],[1,19],zv,nx,ny,bCut)   #number ranges are the bounds to do the plotting (0-20)

        # In[ ]:

        print('Calculating quantities')
        flhInterp = flt.getInterpolatedQuantity(testFLHDen,grid_spacing)
        flhBzInterp = flt.getInterpolatedQuantity(testFLHDen*bz[:,:,z_photo],grid_spacing)
        flwindInterp = flt.getInterpolatedQuantity(testWindDen,grid_spacing)
        twistInterp = flt.getInterpolatedQuantity(twistDensity,grid_spacing)


        # In[ ]:

        flh,indexesflh = flt.fieldLineIntegratedQuantity(flhInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
        flhBz,indexesflh = flt.fieldLineIntegratedQuantity(flhBzInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
        flw,indexesflw = flt.fieldLineIntegratedQuantity(flwindInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
        twistF,indexestwist = flt.fieldLineIntegratedQuantity(twistInterp,goodSeeds,fieldLinesList,seeds,nx,ny)

        flhs.append(flh)
        flhBzs.append(flhBz)
        flws.append(flw)
        twistFs.append(twistF)

# In[44]:
if True:
    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(10, 10))

    vmin = 1e6; vmax =  -1e6
    for strat_flag in range(nstrats):
        vmin = min(vmin, np.percentile(flhs[strat_flag], 1))
        vmax = max(vmin, np.percentile(flhs[strat_flag], 99))
    for strat_flag in range(nstrats):
        im = axes[0,strat_flag].imshow(flhs[strat_flag].T,origin='lower',cmap='seismic', vmin=vmin, vmax =vmax)
        fig.colorbar(im, ax=axes[0,strat_flag])
        axes[0,strat_flag].set_title('FLH' + paths[strat_flag])
    im = axes[0,2].imshow(flhs[1].T - flhs[0].T,origin='lower',cmap='seismic')
    fig.colorbar(im, ax=axes[0,2])
    axes[0,strat_flag].set_title('FLH Difference')

    # In[45]:
    vmin = 1e6; vmax =  -1e6
    for strat_flag in range(nstrats):
        vmin = min(vmin, np.percentile(flws[strat_flag], 1))
        vmax = max(vmin, np.percentile(flws[strat_flag], 99))

    for strat_flag in range(nstrats):
        im = axes[1,strat_flag].imshow(flws[strat_flag].T,origin='lower',cmap='seismic', vmin=vmin, vmax =vmax)
        fig.colorbar(im, ax=axes[1,strat_flag])
        axes[1,strat_flag].set_title('FLW' + paths[strat_flag])

    im = axes[1,2].imshow(flws[1].T - flws[0].T,origin='lower',cmap='seismic')
    fig.colorbar(im, ax=axes[1,2])
    axes[1,strat_flag].set_title('FLW Difference')

    # In[46]:

    vmin = 1e6; vmax =  -1e6
    for strat_flag in range(nstrats):
        vmin = min(vmin, np.percentile(flhBzs[strat_flag], 1))
        vmax = max(vmin, np.percentile(flhBzs[strat_flag], 99))
    for strat_flag in range(nstrats):
        im = axes[2,strat_flag].imshow(flhBzs[strat_flag].T,origin='lower',cmap='seismic', vmin=vmin, vmax =vmax)
        fig.colorbar(im, ax=axes[2,strat_flag])
        axes[2,strat_flag].set_title('WFLH' + paths[strat_flag])

    im = axes[2,2].imshow(flhBzs[1].T - flhBzs[0].T,origin='lower',cmap='seismic')
    fig.colorbar(im, ax=axes[2,2])
    axes[2,strat_flag].set_title('WFLW Difference')

    # In[ ]:

    vmin = 1e6; vmax =  -1e6
    for strat_flag in range(nstrats):
        vmin = min(vmin, np.percentile(twistFs[strat_flag], 1))
        vmax = max(vmin, np.percentile(twistFs[strat_flag], 99))
    for strat_flag in range(nstrats):
        im = axes[3,strat_flag].imshow(twistFs[strat_flag].T,origin='lower',cmap='seismic', vmin=vmin, vmax =vmax)
        fig.colorbar(im, ax=axes[3,strat_flag])
        axes[3,strat_flag].set_title('Twists' + paths[strat_flag])

    im = axes[3,2].imshow(twistFs[1].T - twistFs[0].T,origin='lower',cmap='seismic')
    fig.colorbar(im, ax=axes[3,2])
    axes[3,strat_flag].set_title('Twists Difference')



    plt.suptitle('Snap id %d' % snap_id)

    plt.tight_layout()
    #plt.show()
    plt.savefig('./quantity_plots/%d.png' % snap_id)
    plt.close()

    # In[ ]:

