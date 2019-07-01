import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import pi, k

num = 1000 #Number of particles
name_of_particle = 'ten_micron_patricle'
r = 0.000016 #radi of particles
M = 561 * r**3  #Mass of particles
N = 651 #loops
p = 10
q = 65
dt = 1
T = 360 #Water Temp
X = 1.00 #x size
Y = 1.00 #y size
Z = 0.15 #z sixe
vis = 0.0003142
g = 9.8
gama = 6*pi*vis*r/M
gat = gama*dt
b = 2*gama*k*T/M #Betta^2

def meanxy(i,v0):
    if v0 == 0:
        return 0
    else:
        return (v0*(1-np.exp(-1*gat*i))/gama)

def vxt(v0):
    a = (b/2*gama)*(1-np.exp(-2*gama*dt))
    v = v0*np.exp(-1*gama*dt) + ((np.sqrt(a))*np.random.normal(0,1))
    return v

def meanz(i,v0):
    m1 = v0*(1-np.exp(-1*gat*i))/gama
    m2 = g*((i*gat)+np.exp(-1*i*gat)-1)/gama**2
    return m1-m2

def var(i):
    q = b/(gama**3)
    v = gat*i - 2*(1-np.exp(-1*gat*i))+0.5*(1-np.exp(-2*gat*i))
    return q*v

def xyt(x0,i,v0,ma):
    #x = x0 + np.random.normal(meanxy(i,v0),var(i))
    x = x0 + vxt(v0)
    #if x >= ma :
    #    return ma
    #elif x <= (-1*ma) :
    #    return (-1*ma)
    #else:
    return x

def zt(z0, i,v0,ma):
    z = z0 + np.random.normal(meanz(i,v0),var(i))
    if z >= 0 :
        return 0
    elif z <= (-1*ma) :
        return (-1*ma)
    else:
        return z


x = np.zeros((N,num))
y = np.zeros((N,num))
z = np.zeros((N,num))


x[0][:] = np.random.uniform(-1*X,X,num)
y[0][:] = np.random.uniform(-1*Y,Y,num)
z[0][:] = np.random.uniform(-1*Z,0,num)

vx = np.zeros((N,num))
vy = np.zeros((N,num))
vz = np.zeros((N,num))

for i in range(1,N):
    for j in range(0,num):
        vx[i][j] = vxt(vx[i-1][j])
        vy[i][j] = vxt(vy[i - 1][j])
        x[i][j] = xyt(x[i-1][j],i,vx[i-1][j],X)
        y[i][j] = xyt(y[i-1][j],i,vy[i-1][j],Y)
        z[i][j] = zt(z[0][j],i,vz[i-1][j],Z)

for i in range (0,q+1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.set_xlim3d(-10, 10)
    # ax.set_ylim3d(-10,10)
    ax.set_zlim3d(-0.15, 0)
    ax.scatter(x[p*i][:], y[p*i][:], z[p*i][:], c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    title = 'ten micron particles Brownian Motion t=(%i)s' %((p*i*dt))
    plt.title(title)
    plt.savefig(title+'.png')
    plt.show()


file_namex = 'X, N %d, dt %d, T %d, vis %d, %s, %d .csv' %(N,dt,T,vis,name_of_particle,num)
np.savetxt(file_namex, np.column_stack((x)), delimiter=",", fmt='%s')
file_namey = 'Y, N %d, dt %d, T %d, vis %d, %s, %d .csv' %(N,dt,T,vis,name_of_particle,num)
np.savetxt(file_namey, np.column_stack((y)), delimiter=",", fmt='%s')
file_namex = 'Z, N %d, dt %d, T %d, vis %d, %s, %d .csv' %(N,dt,T,vis,name_of_particle,num)
np.savetxt(file_namex, np.column_stack((z)), delimiter=",", fmt='%s')
