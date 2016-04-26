__author__ = 'Walrus'

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy
from PIL import Image

def buildUpB(rho, dt, dx, dy, u, v):
    b = numpy.zeros_like(u)
    b[1:-1,1:-1]=rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))-\
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2-\
                      2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))-\
                      ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)

    ####Periodic BC Pressure @ x = 2
    b[1:-1,-1]=rho*(1/dt*((u[1:-1,0]-u[1:-1,-2])/(2*dx)+(v[2:,-1]-v[0:-2,-1])/(2*dy))-\
                    ((u[1:-1,0]-u[1:-1,-2])/(2*dx))**2-\
                    2*((u[2:,-1]-u[0:-2,-1])/(2*dy)*(v[1:-1,0]-v[1:-1,-2])/(2*dx))-\
                    ((v[2:,-1]-v[0:-2,-1])/(2*dy))**2)

    ####Periodic BC Pressure @ x = 0
    b[1:-1,0]=rho*(1/dt*((u[1:-1,1]-u[1:-1,-1])/(2*dx)+(v[2:,0]-v[0:-2,0])/(2*dy))-\
                    ((u[1:-1,1]-u[1:-1,-1])/(2*dx))**2-\
                    2*((u[2:,0]-u[0:-2,0])/(2*dy)*(v[1:-1,1]-v[1:-1,-1])/(2*dx))-\
                    ((v[2:,0]-v[0:-2,0])/(2*dy))**2)

    return b

def presPoissPeriodic(p, dx, dy):
    pn = numpy.empty_like(p)

    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

        ####Periodic BC Pressure @ x = 2
        p[1:-1,-1] = ((pn[1:-1,0]+pn[1:-1,-2])*dy**2+(pn[2:,-1]+pn[0:-2,-1])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,-1]

        ####Periodic BC Pressure @ x = 0
        p[1:-1,0] = ((pn[1:-1,1]+pn[1:-1,-1])*dy**2+(pn[2:,0]+pn[0:-2,0])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,0]

        ####Wall boundary conditions, pressure
        p[-1,:] =p[-2,:]     ##dp/dy = 0 at y = 2
        p[0,:] = p[1,:]      ##dp/dy = 0 at y = 0

    return p

##variable declarations
nx = 200
ny = 200
nt = 30
nit = 500
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
X,Y = numpy.meshgrid(x,y)

im = Image.open("mywing.png") #Can be many different formats.
pix = im.load()
nulls = []
print pix[0,0] #Get the RGBA Value of the a pixel of an image
for i in range(ny):
    for j in range(ny):
        if pix[i,j] != (255,255,255,255):
            nulls.append((200-j,i))

##physical variables
rho = 1
nu = .0005
F = 1
dt = .0001

#initial conditions
u = numpy.zeros((ny,nx)) ##create a XxY vector of 0's
un = numpy.zeros((ny,nx)) ##create a XxY vector of 0's

v = numpy.zeros((ny,nx)) ##create a XxY vector of 0's
vn = numpy.zeros((ny,nx)) ##create a XxY vector of 0's

p = numpy.ones((ny,nx)) ##create a XxY vector of 0's
pn = numpy.ones((ny,nx)) ##create a XxY vector of 0's

b = numpy.zeros((ny,nx))

udiff = 1
stepcount = 0

R = 20
c_i = 50
c_j = 50

while udiff > .005:
    un = u.copy()
    vn = v.copy()

    b = buildUpB(rho, dt, dx, dy, u, v)
    p = presPoissPeriodic(p, dx, dy)

    u[1:-1,1:-1] = un[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
        vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
        dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])+\
        nu*(dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
        dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))+F*dt

    v[1:-1,1:-1] = vn[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
        dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+\
        nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+\
        (dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])))

    ####Periodic BC u @ x = 2
    u[1:-1,-1] = un[1:-1,-1]-\
        un[1:-1,-1]*dt/dx*(un[1:-1,-1]-un[1:-1,-2])-\
        vn[1:-1,-1]*dt/dy*(un[1:-1,-1]-un[0:-2,-1])-\
        dt/(2*rho*dx)*(p[1:-1,0]-p[1:-1,-2])+\
        nu*(dt/dx**2*(un[1:-1,0]-2*un[1:-1,-1]+un[1:-1,-2])+\
        dt/dy**2*(un[2:,-1]-2*un[1:-1,-1]+un[0:-2,-1]))+F*dt

    ####Periodic BC u @ x = 0
    u[1:-1,0] = un[1:-1,0]-\
        un[1:-1,0]*dt/dx*(un[1:-1,0]-un[1:-1,-1])-\
        vn[1:-1,0]*dt/dy*(un[1:-1,0]-un[0:-2,0])-\
        dt/(2*rho*dx)*(p[1:-1,1]-p[1:-1,-1])+\
        nu*(dt/dx**2*(un[1:-1,1]-2*un[1:-1,0]+un[1:-1,-1])+\
        dt/dy**2*(un[2:,0]-2*un[1:-1,0]+un[0:-2,0]))+F*dt

    ####Periodic BC v @ x = 2
    v[1:-1,-1] = vn[1:-1,-1]-\
        un[1:-1,-1]*dt/dx*(vn[1:-1,-1]-vn[1:-1,-2])-\
        vn[1:-1,-1]*dt/dy*(vn[1:-1,-1]-vn[0:-2,-1])-\
        dt/(2*rho*dy)*(p[2:,-1]-p[0:-2,-1])+\
        nu*(dt/dx**2*(vn[1:-1,0]-2*vn[1:-1,-1]+vn[1:-1,-2])+\
        (dt/dy**2*(vn[2:,-1]-2*vn[1:-1,-1]+vn[0:-2,-1])))

    ####Periodic BC v @ x = 0
    v[1:-1,0] = vn[1:-1,0]-\
        un[1:-1,0]*dt/dx*(vn[1:-1,0]-vn[1:-1,-1])-\
        vn[1:-1,0]*dt/dy*(vn[1:-1,0]-vn[0:-2,0])-\
        dt/(2*rho*dy)*(p[2:,0]-p[0:-2,0])+\
        nu*(dt/dx**2*(vn[1:-1,1]-2*vn[1:-1,0]+vn[1:-1,-1])+\
        (dt/dy**2*(vn[2:,0]-2*vn[1:-1,0]+vn[0:-2,0])))


    ####Wall BC: u,v = 0 @ y = 0,2
    u[0,:] = 0
    u[-1,:] = 0
    v[0,:] = 0
    v[-1,:]=0

    for i in nulls:
        u[i] = 0
        v[i] = 0
        p[i] = 0

    udiff = (numpy.sum(u)-numpy.sum(un))/numpy.sum(u)
    print udiff
    stepcount += 1

fig = pyplot.figure(figsize = (10,10), dpi=100)
#pyplot.quiver(X, Y, u, v);
pyplot.pcolor(X,Y,p)
pyplot.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
#pyplot.scatter(X, Y, p**3, alpha=0.5)
pyplot.show()