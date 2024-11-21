#!/usr/bin/env python
# coding: utf-8

import numpy as np 
from matplotlib import pyplot as plt
from scipy.fft import ifft2,fft2,fftfreq, fftshift, ifftshift,ifft,fft
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simps
from streamtracer import StreamTracer, VectorGrid
#from numba import jit




def norm2d(vec):
    mag = np.linalg.norm(vec)
    if (mag > 0.0):
        v = vec/mag
    else:
        v = np.array([0, 0])
    return np.array([v[0],v[1],0.0])


def getFrequencyMatrix(ncells,spacing):
    freqlist1da =np.roll(np.linspace(-ncells[0]/2,ncells[0]/2-1,ncells[0]),round(ncells[0]/2))/(ncells[0]*spacing[0])
    freqlist1db =np.roll(np.linspace(-ncells[1]/2,ncells[1]/2-1,ncells[1]),round(ncells[1]/2))/(ncells[1]*spacing[1])
    return np.array([np.array([np.array([2.0*np.pi*freqlist1da[i],2.0*np.pi*freqlist1db[j]]) for j in range(len(freqlist1db))]) for i  in range(len(freqlist1da))]);       

def getFrequencyMatrixVert(ncells,spacing):
    freqlist1 = np.roll(np.linspace(-ncells[2]/2,ncells[2]/2-1,ncells[2]),round(ncells[2]/2))/(ncells[2]*spacing[2])
    return np.array([2.0*np.pi*freqlist1[i] for i  in range(len(freqlist1))]);       


def getAFastSingle(b,ncells,spacing):
    fm = getFrequencyMatrix(ncells,spacing);
    # make the basis
    kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
    kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])
    # note in the k matrix below the k=0 element is set to one so we can divide by it.
    k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T    
    A = np.zeros(b.shape)
    B0fx = np.zeros(b.shape[2],dtype=np.complex128)
    B0fy = np.zeros(b.shape[2],dtype=np.complex128)
    aftx = np.zeros([b.shape[0],b.shape[1],b.shape[2]],dtype=np.complex128)
    afty  =np.zeros([b.shape[0],b.shape[1],b.shape[2]],dtype=np.complex128)
    aftz = np.zeros([b.shape[0],b.shape[1],b.shape[2]],dtype=np.complex128)
    for i in range(b.shape[2]):
        fbx = fft2(b[:,:,i,0]); fby =fft2(b[:,:,i,1]); fbz = fft2(b[:,:,i,2])
        #set the zero element of the transform for the final vertical part
        B0fx[i] = fbx[0][0]
        B0fy[i] = fby[0][0]
        akperp = -1j*fbz/k
        ## fix i =j  element
        #print(akperp[0][0],fbz[0][0])
        akperp[0][0]=0.0
        akw = 1j*(-(kparr[:,:,1])*fbx + (kparr[:,:,0])*fby)/k
        ## fix i =j  element
        akw[0][0]=0.0    
        aftx[:,:,i] = akperp*kperp[:,:,0]
        afty[:,:,i] = akperp*kperp[:,:,1]
        aftz[:,:,i] = akperp*kperp[:,:,2]+akw
    # now do the vertical transforms
    kz= getFrequencyMatrixVert(ncells,spacing)
    fB0fx = fft(B0fx)
    fB0fy = fft(B0fy)
    kz[0]=1.0
    ax00 = -1j*fB0fy/kz
    ay00 =  1j*fB0fx/kz
    ax00[0]=0.0
    ay00[0]=0.0
    ax00 = ifft(ax00)
    ay00 = ifft(ay00)
    #finally transform back to real space
    for i in range(b.shape[2]):
        aftx[0][0] = ax00[i]
        afty[0][0] = ay00[i]
        ax = ifft2(aftx[:,:,i])
        ay = ifft2(afty[:,:,i])
        az = ifft2(aftz[:,:,i])
        A[:,:,i,0] = np.real(ax)
        A[:,:,i,1] = np.real(ay)
        A[:,:,i,2] = np.real(az)
    return A


def getHelicity(bx,by,bz,ax,ay,az):
    return np.sum(ax*bx + ay*by + az*bz)

def getFLHDen(bx,by,bz,ax,ay,az):
    bmag =np.sqrt(bx*bx + by*by + bz*bz)
    return (ax*bx + ay*by + az*bz)/bmag

def getFLHDenSingle(b,a):
    bmag = vector_norm(b)
    multProd = a*b
    return vector_norm(multProd)/bmag

def getFLHDenSingle(b,a):
    bmag = vector_norm(b)
    multProd = a*b
    return vector_norm(multProd)/bmag

def getAWindFast(bx,by,bz):
    fm = getFrequencyMatrix(ncells,spacing);
    # make the basis
    normFunc = np.vectorize(norm2d)
    kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
    kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])
    # note in the k matrix below the k=0 element is set to one so we can divide by it.
    k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T    
    Ax = np.zeros(bx.shape)
    Ay = np.zeros(by.shape)
    Az = np.zeros(bz.shape)
    for i in range(bz.shape[2]):
        # locate the zeros in 
        zerolocs = np.argwhere(np.abs(bz[:,:,i])<0.00000000001)
        bzSafe= bz[:,:,i]
        for j in range(len(zerolocs)):
            bzSafe[zerolocs[j][0],zerolocs[j][1]]=1.0
        tx = bx[:,:,i]/bzSafe
        ty = by[:,:,i]/bzSafe
        tz = np.ones(tx.shape)
        for j in range(len(zerolocs)):
            tx[zerolocs[j][0],zerolocs[j][1]]=0.0
            ty[zerolocs[j][0],zerolocs[j][1]]=0.0
        fbx = fft2(tx); fby =fft2(ty); fbz = fft2(tz)
        akperp = -1j*fbz/k
        ## fix i =j  element
        akperp[0][0]=0.0
        akw = 1j*(-kparr[:,:,1].T*fbx + kparr[:,:,0].T*fby)/k
        ## fix i =j  element
        akw[0][0]=0.0    
        aftx = akperp*kperp[:,:,0]
        afty = akperp*kperp[:,:,1]
        aftz = akperp*kperp[:,:,2]+akw
        ax = ifft2(aftx)
        ay = ifft2(afty)
        az = ifft2(aftz)
        Ax[:,:,i] = np.real(ax)
        Ay[:,:,i] = np.real(ay)
        Az[:,:,i] = np.real(az)
    return Ax,Ay,Az

def getCurve(tracer,index,zv):
    cv =tracer.xs[index]
    mask = (cv[:, 2]> zv)
    cv =cv[mask, :]
    return cv.T[0],cv.T[1],cv.T[2]
    

def getInterpolatedField(bx,by,bz,dx,dy,dz):
    x = np.linspace(0,dx*(bx.shape[0]-1),bx.shape[0])
    y = np.linspace(0,dy*(bx.shape[1]-1),bx.shape[1])
    z = np.linspace(0,dz*(bx.shape[2]-1),bx.shape[2])
    xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)
    bxInterp = RegularGridInterpolator((x, y, z), bx)
    byInterp = RegularGridInterpolator((x, y, z), by)
    bzInterp = RegularGridInterpolator((x, y, z), bz)
    return bxInterp,byInterp,bzInterp

def getInterpolatedFieldSingle(b,dx,dy,dz):
    x = np.linspace(0,dx*(b.shape[0]-1),b.shape[0])
    y = np.linspace(0,dy*(b.shape[1]-1),b.shape[1])
    z = np.linspace(0,dz*(b.shape[2]-1),b.shape[2])
    xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)
    bxInterp = RegularGridInterpolator((x, y, z), b[:,:,:,0], bounds_error=False, fill_value=0.0)
    byInterp = RegularGridInterpolator((x, y, z), b[:,:,:,1], bounds_error=False, fill_value=0.0)
    bzInterp = RegularGridInterpolator((x, y, z), b[:,:,:,2], bounds_error=False, fill_value=0.0)
    return bxInterp,byInterp,bzInterp

def makePoints(x,y,z):
    return [np.array([i,j,k]) for k in z for j in y for i in x]

def getSlice(interpField,xlims,ylims,nx,ny,zval):
    x = np.linspace(xlims[0],xlims[1],nx)
    y = np.linspace(ylims[0],ylims[1],ny)
    pts =makePoints(x,y,np.linspace(zval,zval,1))
    return interpField(pts).reshape(ny,nx).T
    
def createSeedGridFlat(nx,ny,xrange,yrange,zv):
    xv = np.linspace(xrange[0],xrange[1],nx)
    yv = np.linspace(yrange[0],yrange[1],ny)
    return np.array([[x,y,zv] for x in xv for y in yv])
    
def getFieldLinesStrCheck(seeds,cutoff,bxInterp,byInterp,bzInterp):
    bxset = bxInterp(seeds)  
    byset = byInterp(seeds) 
    bzset = bzInterp(seeds)
    bset = np.sqrt(np.square(bxset)+np.square(byset)+np.square(bzset))    
    seedsout=[]
    # weed out the weakling field
    for i in range(len(seeds)):
        if bset[i]>cutoff:
            seedsout.append(seeds[i])
    return np.array(seedsout)

def getCurve(tracer,index,zv):
    cv =tracer.xs[index]
    mask = (cv[:, 2]> zv)
    cv =cv[mask, :]
    return cv.T[0],cv.T[1],cv.T[2]

def getCurveVec(tracer,index,zv):
    cv =tracer.xs[index]
    mask = (cv[:, 2]> zv)
    cv =cv[mask, :]
    return cv
    
def getInterpolatedQuantity(data,spacings):
    x = np.linspace(0,spacings[0]*(data.shape[0]-1),data.shape[0])
    y = np.linspace(0,spacings[1]*(data.shape[1]-1),data.shape[1])
    z = np.linspace(0,spacings[2]*(data.shape[2]-1),data.shape[2])
    xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)
    return RegularGridInterpolator((x, y, z), data)

def getInterpolatedQuantity2D(quantity, dx, dy):
    x = np.linspace(0, dx * (quantity.shape[0] - 1), quantity.shape[0])
    y = np.linspace(0, dy * (quantity.shape[1] - 1), quantity.shape[1])
    return RegularGridInterpolator((x, y), quantity)

def fieldLineIntegration(fl,interpolatedQuantity):
    if (len(fl) > 3):
        difSet= np.delete(np.roll(fl,-1,axis=0)-fl,-1,axis=0)
        #dsSet=np.array([np.linalg.norm(v) for v in difSet])
        dsSet= np.linalg.norm(difSet,axis=1)
        integratedQuantity = interpolatedQuantity(fl)
        return simps(dsSet*np.vstack([integratedQuantity[1:], integratedQuantity[:-1]]).mean(axis=0))
    else:
        return 0.0

def fieldLineIntegratedQuantity(interpolatedQuantity,bx,by,bz,nx,ny,lx,ly,goodSeeds,fieldLines):
    outputDensity = np.zeros((nx,ny))
    dx = lx/(nx-1)
    dy = ly/(ny)
    for i in range(len(goodSeeds)):
        fl = fieldLines[i]
        val = fieldLineIntegration(fl,interpolatedQuantity)
        xind = round(goodSeeds[i][0]/dx)
        yind = round(goodSeeds[i][1]/dy)
        outputDensity[xind][yind] = val
    return outputDensity

        
def createSingleField(bx,by,bz):
        fieldOut = np.zeros((bx.shape[0],bx.shape[1],bx.shape[2],3));
        fieldOut[:,:,:,0]=bx
        fieldOut[:,:,:,1]=by
        fieldOut[:,:,:,2]=bz
        return fieldOut

def curl(b_field,steps):  # (x, y, z)
    _, dFx_dy, dFx_dz = np.gradient(b_field[..., 0],axis=[0, 1, 2], edge_order=2)
    dFy_dx, _, dFy_dz = np.gradient(b_field[..., 1],axis=[0, 1, 2], edge_order=2)
    dFz_dx, dFz_dy, _ = np.gradient(b_field[..., 2],axis=[0, 1, 2], edge_order=2)

    rot_x = dFz_dy/steps[1] - dFy_dz/steps[2]
    rot_y = dFx_dz/steps[2] - dFz_dx/steps[0]
    rot_z = dFy_dx/steps[0] - dFx_dy/steps[1]

    return createSingleField(rot_x,rot_y,rot_z)

def curlScale(b_field):  # (x, y, z)
    _, dFx_dy, dFx_dz = np.gradient(b_field[..., 0], axis=[0, 1, 2], edge_order=2)
    dFy_dx, _, dFy_dz = np.gradient(b_field[..., 1], axis=[0, 1, 2], edge_order=2)
    dFz_dx, dFz_dy, _ = np.gradient(b_field[..., 2], axis=[0, 1, 2], edge_order=2)

    rot_x = dFz_dy - dFy_dz
    rot_y = dFx_dz - dFz_dx
    rot_z = dFy_dx - dFx_dy

    return rot_x, rot_y, rot_z

#@jit(nopython=True)
def GetKXKY(kx,ky,nx,ny,Kx,Ky): 
    for i in range(0,nx):
        for r in range(0,ny):
            if ((kx[i]**2+ky[r]**2) !=0):
                Kx[r,i]= kx[i]/np.sqrt(kx[i]**2+ky[r]**2)
                Ky[r,i]= ky[r]/np.sqrt(kx[i]**2+ky[r]**2)

def getPotential(bz):
    nx = bz.shape[0]
    ny = bz.shape[1]
    freqx = fftfreq(nx)
    freqy = fftfreq(ny)
    kx= fftshift(freqx)
    ky= fftshift(freqy)
    Kx=np.zeros((ny,nx))
    Ky=np.zeros((ny,nx))   
    GetKXKY(kx,ky,nx,ny,Kx,Ky)
    Bxp=np.zeros((ny,nx), dtype=float)
    Byp=np.zeros((ny,nx), dtype=float)
    L=nx*ny         
    Bxp = np.zeros(bz.shape)
    Byp = np.zeros(bz.shape)
    for i in range(bz.shape[2]):
        FFTB=fftshift(fft2(bz[:,:,i].T))
        Mx=-FFTB*Kx*1j
        My=-FFTB*Ky*1j
        Bxp_slice=ifft2(ifftshift(Mx)).real
        Byp_slice=ifft2(ifftshift(My)).real
        Bxp[:,:,i] = Bxp_slice.T
        Byp[:,:,i] = Byp_slice.T
    return Bxp,Byp

def lorentz_force(b_field, j_field=None):
    j_field = j_field if j_field is not None else curl(b_field)
    l = np.cross(j_field, b_field, axis=-1)
    return l

def vector_norm(vector):
    return np.sqrt((vector ** 2).sum(-1))

def divergenceForWind(b_field,steps):
    # get in plane components
    dFx_dx = np.gradient(b_field[..., 0], axis=[0], edge_order=2)/steps[0]
    dFy_dy = np.gradient(b_field[..., 1], axis=[1], edge_order=2)/steps[1]
    vertDirecs = np.sign(b_field[..., 2])
    dFz_dz = np.gradient(vertDirecs, axis=[2], edge_order=2)/steps[2]
    return  dFx_dx+dFy_dy+dFz_dz

def divergence(b_field,steps):
    # get in plane components
    
    dFx_dx = np.gradient(b_field[..., 0],steps[0],axis=0, edge_order=2)
    dFy_dy = np.gradient(b_field[..., 1],steps[1],axis=1, edge_order=2)
    dFz_dz = np.gradient(b_field[..., 2],steps[2],axis=2, edge_order=2)
    return  dFx_dx+dFy_dy+dFz_dz


def magForFt(vec):
    mag = np.linalg.norm(vec)
    return mag

def addDivergenceCleaningTerm(Bfield,ncells,spacing):
    divField = divergence(Bfield,spacing)
    fm = getFrequencyMatrix(ncells,spacing);
    # note in the k matrix below the k=0 element is set to one so we can divide by it.
    k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T    
    divFreeField = np.zeros(Bfield.shape)
    for i in range(Bfield.shape[2]):
        # get the transformed divergence 
        divF = fft2(divField[:,:,i])
        divergenceFactor = -1j*divF/k
        divergenceTerm  = np.real(ifft2((fm*divergenceFactor[:,:,np.newaxis]).transpose(2,0,1)))
        ## fix i =j  element
        divergenceTerm[0][0][0] = 0.0
        divergenceTerm[1][0][0] = 0.0
        #divFreeField[:,:,i,0] = Bfield[:,:,i,0]+divergenceTerm[:,:,0]
        #divFreeField[:,:,i,1] = Bfield[:,:,i,1]+divergenceTerm[:,:,1]
        divFreeField[:,:,i,0] = Bfield[:,:,i,0]+divergenceTerm[0]
        divFreeField[:,:,i,1] = Bfield[:,:,i,1]+divergenceTerm[1]
        divFreeField[:,:,i,2] = np.sign(Bfield[:,:,i,2])
    return divFreeField

def addDivergenceCleaningTermTest(Bfield,spacing):
    divField = divergenceForWind(Bfield,spacing)
    fm = getFrequencyMatrix(ncells,spacing);
    # note in the k matrix below the k=0 element is set to one so we can divide by it.
    k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T    
        # get the transformed divergence 
    divF = fft2(divField[:,:,0])
    divergenceFactor = -1j*divF/k
        #divergenceTerm  = np.real(ifft2(np.array([np.array([divergenceFactor[j][k]*fm[j][k] for j in range(len(fm))]) for k  in range(len(fm[0]))]).T))
    return fm,divergenceFactor

def unitSpeedField(bField,tol):
    bz = bField[:,:,:,2]
    bz[abs(bz) < tol] = 1.0
    bFieldOut = np.zeros(bField.shape)
    bFieldOut[:,:,:,0]= bField[:,:,:,0]/bz
    bFieldOut[:,:,:,1]= bField[:,:,:,1]/bz
    bFieldOut[:,:,:,2]= np.sign(bField[..., 2])
    return bFieldOut


def fieldMags(BxInterp,ByInterp,BzInterp,curve):
    bxvals =BxInterp(curve)
    byvals =ByInterp(curve)
    bzvals =BzInterp(curve)
    return np.sqrt(bxvals*bxvals +  byvals*byvals +bzvals*bzvals)

def chopWeakFieldLine(BxInterp,ByInterp,BzInterp,curve,cutoff):
    bxvals =BxInterp(curve)
    byvals =ByInterp(curve)
    bzvals =BzInterp(curve)
    mags = np.sqrt(bxvals*bxvals +  byvals*byvals +bzvals*bzvals)
    if(len(curve)>1):
        cut =np.argmax(mags<cutoff)
        if cut != 0:
            return curve[:cut]
        else: 
            if mags[0]<cutoff:
                return []
            else:
                return curve
    else:
        return curve
    
    
    
def prepareCurves(field,bxInterp,byInterp,bzinterp,gridSpacing,xrange,yrange,zv,nx,ny,bcut):
        nsteps = 10000
        step_size = 0.05
        grid = VectorGrid(field,gridSpacing)
        tracer = StreamTracer(nsteps,step_size)
        seeds =createSeedGridFlat(nx,ny,xrange,yrange,zv)
        goodSeeds= getFieldLinesStrCheck(seeds,bcut,bxInterp,byInterp,bzinterp)
        #Get the field lines, needs to be done in batches else sometimes program runs out of memory
        noFieldLinePerSec= int(len(goodSeeds)/8)
        fieldLinesList=[]
        for j in range(1,8):
            if j==7:
                subFieldSet =goodSeeds[(6)*noFieldLinePerSec:]
            else:
                subFieldSet =goodSeeds[(j-1)*noFieldLinePerSec:j*noFieldLinePerSec]
            tracer.trace(subFieldSet,grid)
            for i in range(len(subFieldSet)):
                curve= getCurveVec(tracer,i,zv)
                fieldLinesList.append(np.array(curve))
        # sometimes you get odd field lines due to weak field, we cut this off at 10% of bcut
        fieldLinesList = [chopWeakFieldLine(bxInterp,byInterp,bzinterp,fl,bcut/10.0) for fl in fieldLinesList] 
        return fieldLinesList,goodSeeds,seeds 
    
def fieldLineIntegratedQuantity(interpolatedQuantity,goodSeeds,fieldLines,seeds,nx,ny):
    outputDensity = np.zeros((nx,ny))
    indexes = np.zeros((nx,ny))
    xdifs = np.abs((np.roll(seeds,1,0) - seeds).T[0])
    dx = np.unique(xdifs)[1]
    dy = (seeds[1]-seeds[0])[1]
    for i in range(len(fieldLines)):
        fl = fieldLines[i]
        val = fieldLineIntegration(fl,interpolatedQuantity)
        xind = round((goodSeeds[i][0]-seeds[0][0])/dx)
        yind = round((goodSeeds[i][1]-seeds[0][1])/dy)
        outputDensity[xind][yind] = val
        indexes[xind][yind] = i
    return outputDensity,indexes

def aConst(bvalue,xymeshIn,pt,dA,indexIgnore):
    xymesh = np.delete(xymeshIn,indexIgnore,axis=0)
    difs = pt-xymesh
    denominators = np.einsum('ij,ij -> i', difs, difs)
    numerators = np.flip(difs*np.array([1.0,-1.0]),axis=1)
    return np.array([np.sum(bvalue*dA*((numerators.T)[0])/denominators),np.sum(bvalue*dA*((numerators.T)[1])/denominators),0.0])

def AConst(bvalue,xymesh,dA,grid):
    aslice =  (1.0/(2.0*np.pi))*np.array([aConst(bvalue,xymesh,xymesh[i],dA,i) for i in range(len(xymesh))])
    AconsGrid = np.zeros((grid[0],grid[1],grid[2],3))
    print(grid[0], grid[1], grid[2])
    for i in range(grid[2]):
        AconsGrid[:,:,i] = aslice.reshape((grid[0],grid[1],3))
    return AconsGrid
    
def twistDen(bField,curlField,tol):
    bmag = vector_norm(bField*bField)
    currDen = vector_norm(bField*curlField)
    currDen[abs(bmag) < tol] = 0.0
    twist = currDen/bmag
    return twist                           
