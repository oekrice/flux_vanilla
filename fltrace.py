import numpy as np
import matplotlib.pyplot as plt
import random
import os

import pyvista as pv
pv.start_xvfb()


class trace_fieldlines():
    def __init__(self, grid, bx,by,bz,save = -1,plot_vista = True,plot_notvista=False):
        print('Tracing field lines...')
        self.xs = np.linspace(grid.x0,grid.x1,grid.nx+1)
        self.ys = np.linspace(grid.y0,grid.y1,grid.ny+1)
        self.zs = np.linspace(grid.z0,grid.z1,grid.nz+1)

        self.nx = grid.nx
        self.ny = grid.ny
        self.nz = grid.nz

        self.x0 = grid.x0; self.x1 = grid.x1
        self.y0 = grid.y0; self.y1 = grid.y1
        self.z0 = grid.z0; self.z1 = grid.z1

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

        self.dx = np.sum(self.xs[1:] - self.xs[:-1])/len(self.xs[1:])
        self.dy = np.sum(self.ys[1:] - self.ys[:-1])/len(self.ys[1:])
        self.dz = np.sum(self.zs[1:] - self.zs[:-1])/len(self.zs[1:])

        self.bx = bx; self.by = by; self.bz = bz
        self.ds = 0.1*min(self.dx,self.dy,self.dz)

        z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

        #Find start positions of each field line. For now a regular grid on the lower boundary
        if False: #Cartesian tracing grid
            nxs = 10; nys = 12
            self.starts = []

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)

            for i in range(nxs):
                for j in range(nys):
                    self.starts.append([xis[i],yjs[j],1e-6])

        elif False:   #Polar plotting grid
            self.starts = []
            #Trace from the top
            nrs = 10; nthetas = 2
            ris = np.linspace(3.0/nrs,self.x0-3.0/nrs,nrs)
            tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
            for i in range(nrs):
                for j in range(nthetas):
                    self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),self.z1-1e-6])
            #And from the bottom (for the interior ones only)
            nrs = 20; nthetas = 2
            ris = np.linspace(0.5*self.x0/nrs,self.x0-0.5*self.x0/nrs,nrs)
            tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
            for i in range(nrs):
                for j in range(nthetas):
                    self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),1e-6])

        elif False:   #grid for plotting the emerging flux. Flux rope is in y direction so try that.
            #Background field isn't much help - only start when field strength is decent'
            self.starts = []

            if False:
                nys = 5; nzs = 5
                self.starts = []

                yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)
                zks = np.linspace(self.z0+1e-6,self.z1-1e-6,nzs)

                for j in range(nys):
                    for k in range(nzs):
                        self.starts.append([self.x0+1e-6,yjs[j],zks[k]])

            nxs = 100; nzs = 100

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            zks = np.linspace(self.z0+1e-6,self.z1-1e-6,nzs)


            for i in range(nxs):
                for k in range(nzs):

                    xp = int((self.nx)*(xis[i] - self.x0)/(self.x1 - self.x0))
                    yp = int((self.ny)*(self.y0+1e-6 - self.y0)/(self.y1 - self.y0))
                    zp = int((self.nz)*(zks[k] - self.z0)/(self.z1 - self.z0))


                    if np.abs(by[xp,yp,zp]) > 1e-2:
                        self.starts.append([xis[i],self.y0+1e-6,zks[k]])


        elif False:   #plot field lines from the photosphere, where the magnetic field strength is above a certain threshold (0.025 maybe?)
            nlines = 10   #line distribution based on absolute magnetic field strength. More field lines where the field is stronger

            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            self.starts = []

            nxs = 40; nys = 40

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)

            for i in range(nxs):
                for j in range(nys):

                    xp = int((self.nx)*(xis[i] - self.x0)/(self.x1 - self.x0))
                    yp = int((self.ny)*(yjs[j] - self.y0)/(self.y1 - self.y0))
                    zp = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

                    if np.abs(bz[xp,yp,zp]) > 0.025:
                        self.starts.append([xis[i],yjs[j],10.0])
            print('Tracing', len(self.starts), 'lines')

        else:
            """Plot field lines from a slice well above the photosphere, to see if a rope is forming"""

            nlines = 10   #line distribution based on absolute magnetic field strength. More field lines where the field is stronger

            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            self.starts = []

            nxs = 40; nys = 40; nzs = 40

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)
            zks = np.linspace(10.0,self.z1-1e-6,nzs)

            for i in range(nxs):
                j = nys//2
                for k in range(nzs):

                    xp = int((self.nx)*(xis[i] - self.x0)/(self.x1 - self.x0))
                    yp = int((self.ny)*(yjs[j] - self.y0)/(self.y1 - self.y0))
                    zp = int((self.nz)*(zks[k] - self.z0)/(self.z1 - self.z0))

                    if np.abs(by[xp,yp,zp]) > 0.01:
                        self.starts.append([xis[i],yjs[j],zks[k]])
            print('Tracing', len(self.starts), 'lines')

        self.lines = []
        for start in self.starts:
            doplot = True
            for updown in [-1,1]:
                line = self.trace_line(start,updown=updown)
                if len(line) < 2:
                    continue

                if line[-2][1] < self.y0+1e-3:
                    continue



                if doplot:
                    self.lines.append(line)


        self.lines = []
        for start in self.starts:
            doplot = True
            for updown in [-1,1]:
                line = self.trace_line(start,updown=updown)
                if len(line) < 2:
                    continue

                if line[-2][1] < self.y0+1e-3:
                    continue

                if doplot:
                    self.lines.append(line)


        if plot_notvista:
            #Don't use pyvista, just the normal thing
            plt.figure().add_subplot(projection='3d')
            for line in self.lines:
                line = np.array(line)
                plt.gca().plot(line[:,0],line[:,1],line[:,2],c='black',linewidth=0.1)
            plt.gca().set_xlim(self.x0, self.x1)
            plt.gca().set_ylim(self.y0, self.y1)
            plt.gca().set_zlim(10.0, self.z1)

            if save >= 0:
                plt.savefig('plots/b%04d.png' % save, dpi = 400)
                plt.close()

            else:
                plt.show()

        elif plot_vista:
            #Do use pyvista
            x, y = np.meshgrid(self.xs, self.ys)
            z = 10.0*np.ones((len(self.xs),len(self.ys)))
            surface = pv.StructuredGrid(x, y, z)
            p = pv.Plotter(off_screen=True)
            p.background_color = "black"
            p.add_mesh(surface, scalars= bz[1:-1,1:-1,z_photo], show_edges=False,cmap = 'plasma')
            for line in self.lines:
                line = np.array(line)
                p.add_mesh(pv.Spline(line, len(line)),color='white',line_width=1)


            #p.camera.position = (0,0.0,1.0)
            print(p.camera.position)
            p.camera.position = (420.0,0.0,300.0)

            p.camera.focal_point = (0,0,0.25)

            p.show(screenshot='plots/b%04d.png' % save, window_size = (1000,1000))

    def interpolate_field(self,pt):
        #Outputs the magnetic field vector at point pt, using the individual magnetic field vectors (no unecessary averaging)
        #Find position of point on the grid
        b1 = np.zeros((3))
        xp = (self.nx)*(pt[0] - self.x0)/(self.x1 - self.x0)
        yp = (self.ny)*(pt[1] - self.y0)/(self.y1 - self.y0)
        zp = (self.nz)*(pt[2] - self.z0)/(self.z1 - self.z0)

        #Interpolate bx
        xi = int(xp); yi = int(yp + 0.5); zi = int(zp + 0.5)
        xf = xp - xi; yf = yp+0.5 - yi; zf = zp + 0.5 - zi


        b1[0] = b1[0] + self.bx[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.bx[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.bx[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.bx[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[0] = b1[0] + self.bx[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.bx[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.bx[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.bx[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        #Interpolate by
        xi = int(xp+0.5); yi = int(yp); zi = int(zp + 0.5)
        xf = xp + 0.5 - xi; yf = yp - yi; zf = zp + 0.5 - zi
        b1[1] = b1[1] + self.by[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.by[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.by[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.by[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[1] = b1[1] + self.by[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.by[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.by[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.by[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        #Interpolate bz
        xi = int(xp+0.5); yi = int(yp + 0.5); zi = int(zp)
        xf = xp + 0.5 - xi; yf = yp + 0.5 - yi; zf = zp - zi
        b1[2] = b1[2] + self.bz[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.bz[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.bz[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.bz[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[2] = b1[2] + self.bz[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.bz[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.bz[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.bz[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        return (b1/np.sqrt(np.sum(b1**2)), np.sqrt(np.sum(b1**2)))

    def trace_line(self, start, updown = 0):
        #Traces an individual field line
        go = True
        pts = [[start[0],start[1],start[2]]]
        pt = start
        if not self.inbounds(pt):
            return []
        grad, mag = self.interpolate_field(pts[-1])
        test = pt + self.ds*grad

        if updown == 0:
            if grad[2] > 0.0:
                updown = 1
            else:
                updown = -1
        while go:
            if not self.inbounds(pt):
                go = False
            elif mag < 1e-2:
                go = False
            else:
                grad, mag = self.interpolate_field(pts[-1])
                pt = pt + updown*self.ds*grad
                pts.append([pt[0],pt[1],pt[2]])

        return pts

    def inbounds(self, pt):
        #returns true if the point is still within the domain. Otherwise doesn't
        if pt[0] < self.x0 or pt[0] > self.x1 or pt[1] < self.y0 or pt[1] > self.y1 or pt[2] < 10.0 or pt[2] > self.z1:
            return False
        else:
            return True



