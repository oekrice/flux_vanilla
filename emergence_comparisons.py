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
# Code to compare field line quanitites (twist, helicity) etc. between the stratified and unstratified jet simulations

snap_id = 300   #number of the snap to compare

copy = True #copy over from archie-west

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

paths = ['./Data_15/', './Data_50/', './Data_150/']
titles = ['Strat ratio 15', 'Strat ratio 50', 'Strat ratio 150']

nstrats = 3

for strat_flag in range(nstrats):   #do unstratified (top) and stratified (bottom)
    print('Dealing with field number', strat_flag)
    print('_______________________________________')
    path = paths[strat_flag]
    file2read = netcdf_file(path+'%04d.nc'% snap_id,'r')

    # In[67]:

    bxOg = np.swapaxes(file2read.variables['bx'][:], 0, 2)
    byOg = np.swapaxes(file2read.variables['by'][:], 0, 2)
    bzOg = np.swapaxes(file2read.variables['bz'][:], 0, 2)

    file2read.close()
    # In[69]:

    bx = 0.5*(bxOg[1:,:,:] + bxOg[:-1,:,:])
    by = 0.5*(byOg[:,1:,:] + byOg[:,:-1,:])
    bz = 0.5*(bzOg[:,:,1:] + bzOg[:,:,:-1])

    ncells = int(np.size(bx)**(1/3)) + 1
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

    # In[70]:

    bField = flt.createSingleField(bx,by,bz)
    #bFieldTest = flt.createSingleField(bxRot,byRot,bzRot)
    # In[71]:
    #get curl

    #calculate the winding gauge
    AField = flt.getAFastSingle(bField)

    #calculate the winding gauge for the unit speed field
    usf = flt.unitSpeedField(bField,0.01)
    AUnit= flt.addDivergenceCleaningTerm(usf,grid_spacing)

    AWind = flt.getAFastSingle(AUnit)

    bField = flt.createSingleField(bx,by,bz)
    curlField= flt.curl(bField,grid_spacing)

    # add constant component of A as there be net flux

    bzConst = np.sum(bz[:,:,0])/(dA*(ncells-1)*(ncells-1))
    AConst = flt.AConst(bzConst,points,dA,[ncells, ncells, ncells])
    AField = AField + AConst

    # In[72]:

    # calculate the field line helcity and winding densities

    testFLHDen = flt.getFLHDenSingle(bField,AField)

    testWindDen = flt.getFLHDenSingle(bField,AWind)
    # twist density

    twistDensity = flt.twistDen(bField,curlField,0.1)

    # In[73]:

    z_photo = int((ncells)*(10.0 - zv[0])/(zv[-1] - zv[0]))
    #Check the field comps

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(bField[:,:,z_photo,i], origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.savefig('extra_plots/0_%d.png' % strat_flag)

    plt.close()

    # In[1]:

    #Check the curl comps

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(twistDensity[:,:,z_photo*(i+1)], origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.savefig('extra_plots/1_%d.png' % strat_flag)
    plt.close()


    # In[61]:

    # FLH density, pretty !

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(testFLHDen[:,:,z_photo*(i+1)].T, origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
    plt.savefig('extra_plots/2_%d.png' % strat_flag)

    plt.close()


    # In[62]:

    # FLwind density, pretty !

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
    for i in range(3):
        im = axes[i].imshow(testWindDen[:,:,z_photo*(i+1)].T, origin='lower')
        fig.colorbar(im, ax=axes[i])
    fig.tight_layout()
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
        nx = 400
        ny = 400

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


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
for strat_flag in range(nstrats):
    im = axes[strat_flag].imshow(flhs[strat_flag].T,origin='lower',cmap='seismic')
    axes[strat_flag].set_title(titles[strat_flag])
    fig.colorbar(im, ax=axes[strat_flag])
plt.suptitle('Field-line helicities')
plt.tight_layout()
plt.savefig('extra_plots/flh.png')
plt.close()

# In[45]:

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
for strat_flag in range(nstrats):
    im = axes[strat_flag].imshow(flws[strat_flag].T,origin='lower',cmap='seismic')
    axes[strat_flag].set_title(titles[strat_flag])
    fig.colorbar(im, ax=axes[strat_flag])
plt.suptitle('Field-line windings')
plt.tight_layout()
plt.savefig('extra_plots/flw.png')
plt.close()


# In[46]:


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
for strat_flag in range(nstrats):
    im = axes[strat_flag].imshow(flhBzs[strat_flag].T,origin='lower',cmap='seismic')
    axes[strat_flag].set_title(titles[strat_flag])
    fig.colorbar(im, ax=axes[strat_flag])
plt.suptitle('Weighted Field-line helicities')
plt.tight_layout()
plt.savefig('extra_plots/flhw.png')
plt.close()


# In[ ]:


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
for strat_flag in range(nstrats):
    im = axes[strat_flag].imshow(twistFs[strat_flag].T,origin='lower',cmap='seismic')
    axes[strat_flag].set_title(titles[strat_flag])
    fig.colorbar(im, ax=axes[strat_flag])
plt.suptitle('Twists')
plt.tight_layout()
plt.savefig('extra_plots/tw.png')
plt.close()


# In[ ]:


def getFieldLineMetrics(fileLoc,dimensions,grid_spacing):
    #open the file
    file2read = netcdf.NetCDFFile(fileLoc,'r')
    # read the field variables
    bxOg =file2read.variables['bx'][:]*1
    byOg =file2read.variables['by'][:]*1
    bzOg =file2read.variables['bz'][:]*1
    # smush to centred grid
    bx =0.5*(bxOg[:,:,:(bxOg.shape[2]-1)]+bxOg[:,:,1::bxOg.shape[2]])
    by =0.5*(byOg[:,:(byOg.shape[1]-1)]+byOg[:,1::byOg.shape[1],:])
    bz =0.5*(bzOg[:(bzOg.shape[0]-1),:,:]+bzOg[1::bzOg.shape[0],:,:])

    #Reshape so its 0= x comp 1 =y comp 2 = z comp 
    bx = bx.reshape(dimensions).transpose(2,1,0)
    by = by.reshape(dimensions).transpose(2,1,0)
    bz = bz.reshape(dimensions).transpose(2,1,0)

    #spatial grids    
    xv = np.linspace(-10,10,192)
    yv = np.linspace(-10,10,192)
    zv = np.linspace(0,20,192)
    X, Y = np.meshgrid(xv, yv, indexing='ij')
    # Flatten the grid arrays to form the input to the interpolator
    points = np.vstack([X.ravel(), Y.ravel()]).T
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]
    dz = zv[1]-zv[0]
    dA = dx*dy
    grid_spacing = [dx,dy,dz]
    
    bField = flt.createSingleField(bx,by,bz)
    
    #get curl
    curlField= flt.curl(bField,grid_spacing)

    #calculate the winding gauge 
    AField = flt.getAFastSingle(bField)

    #calculate the winding gauge for the unit speed field 
    usf = flt.unitSpeedField(bField,0.01)
    AUnit= flt.addDivergenceCleaningTerm(usf,grid_spacing)

    AWind = flt.getAFastSingle(AUnit)

    #add constant component of A as there be net flux

    bzConst = np.sum(bz[:,:,0])/(dA*191*191)
    AConst = flt.AConst(bzConst,points,dA,[192,192,192])
    AField = AField + AConst 
    
    # calculate the field line helcity and winding densities

    testFLHDen =flt.getFLHDenSingle(bField,AField)
    testWindDen =flt.getFLHDenSingle(bField,AWind)
    twistDensity = flt.twistDen(bField,curlField,0.1)
    
    # set up the field line tracer
    
    grid_spacing = [20/192,20/192,20/192]

    #set domain lengths

    lx = grid_spacing[0]*192
    ly = grid_spacing[1]*192
    lz = grid_spacing[2]*192

    #set number of points used to calculate the distribution from
    nx = 200
    ny = 200

    # set a minimum strength of field line cut off 
    bCut = 0.01

    # set z value for anchoring plane (use photosphere z index = 116 here)
    zv = 0
    
    #interpolate the field for tracer
    BxInterp,ByInterp,BzInterp = flt.getInterpolatedFieldSingle(bField,grid_spacing[0],grid_spacing[1],grid_spacing[2])
    
    #calculate and prepare curves
    fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[1,19],[1,19],zv,nx,ny,bCut)
    
    bz_surface_interp = BzInterp()
    #interpolate the quantities of interest
    flhInterp = flt.getInterpolatedQuantity(testFLHDen,grid_spacing)
    flhBzInterp = flt.getInterpolatedQuantity(testFLHDen*bz[:,:,0],grid_spacing)
    flwindInterp = flt.getInterpolatedQuantity(testWindDen,grid_spacing)
    twistInterp = flt.getInterpolatedQuantity(twistDensity,grid_spacing)
    
    flh,indexesflh = flt.fieldLineIntegratedQuantity(flhInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    flhBz,indexesflh = flt.fieldLineIntegratedQuantity(flhBzInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    flw,indexesflw = flt.fieldLineIntegratedQuantity(flwindInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    twistF,indexestwist = flt.fieldLineIntegratedQuantity(twistInterp,goodSeeds,fieldLinesList,seeds,nx,ny)

