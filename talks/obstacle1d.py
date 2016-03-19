#!/usr/bin/env python

from pylab import *
from sys import exit

figdebug = False
def figsave(name):
    if figdebug:
        show()  # debug
    else:
        print "saving %s ..." % name
        savefig(name,bbox_inches='tight')

def psi(x):
    s = 5.0 * (x - 0.5)
    s2 = s * s
    ps = 0.6 + 0.1 * s + 0.5 * (0.5 * s2 - 0.25 * s2 * s2)
    delta = 0.15
    return maximum(ps, - delta * (1.0 - exp(-abs(ps)/delta)))

def getu(N):
    x = linspace(0.0,1.0,N+1)
    mpsi = psi(x).max()
    u = mpsi * ones(shape(x))
    u[0] = 0.0
    u[N] = 0.0
    gap = 0.01
    tol = 1.0e-4
    while True:
        uold = u.copy()
        for j in range(1,N):
            u[j] = max(psi(x[j])+gap, 0.5 * (u[j-1] + u[j+1]))
        if abs(u-uold).max() < 1.0e-6:
            break
    return x,u

xL = -0.1
xR = 1.1
yT = 1.0
yB = -0.2

N = 101
x,u = getu(N)

figure(figsize=(10,4))
hold(True)
plot([xL,xR],[0.0,0.0],'k',lw=2.0)
plot([0.0,0.0],[yB,yT],'k',lw=2.0)
plot(x,psi(x),'b',lw=2.0)
plot(x,u,'r',lw=2.0)
#plot(x,u,'r.-',ms=2.0)
#plot(x[1:N],u[1:N],'r.-',ms=15.0)
#text(0.1,0.4,r'$u_j$',color='red',fontsize=20.0)
text(0.1,0.5,r'$u(x)$',color='red',fontsize=20.0)
text(0.2,0.2,r'$\psi(x)$',color='blue',fontsize=20.0)
axis('off')
axis('tight')
axis([xL,xR,yB,yT])
figsave('obstacle1d.pdf')

N = 16
x,u = getu(N)
P = psi(x)
z = u - P

figure(figsize=(10,4))
hold(True)
plot([xL,xR],[0.0,0.0],'k',lw=2.0)
plot([0.0,0.0],[yB,yT],'k',lw=2.0)
xx = linspace(xL,xR,201)
plot(x,z,'r.-',ms=2.0)
plot(x[1:N],z[1:N],'r.-',ms=15.0)
text(0.05,0.5,r'$z_j$',color='red',fontsize=20.0)
F = zeros(shape(x))
dx = 1.0 / N
for j in range(1,N):
    F[j] = - (z[j+1]-2.0*z[j]+z[j-1])/dx**2.0 - (P[j+1]-2.0*P[j]+P[j-1])/dx**2.0
F = 0.025 * F  # scale for fit
plot(x,F,'.-',color='green',ms=2.0)
plot(x[1:N],F[1:N],'.-',color='green',ms=15.0)
text(0.2,0.6,r'$F_j$',color='green',fontsize=20.0)
axis('off')
axis('tight')
axis([xL,xR,-0.1,2.0])
figsave('ncp1d.pdf')

figure(figsize=(10,4))
hold(True)
plot([xL,xR],[0.0,0.0],'k',lw=2.0)
plot([0.0,0.0],[yB,yT],'k',lw=2.0)
xx = linspace(0.0,1.0,201)
plot(xx,psi(xx),'b')
plot(x[1:N],psi(x[1:N]),'b.',ms=15.0)
text(0.2,0.2,r'$\psi(x_j)$',color='blue',fontsize=20.0)
plot(x,u,'r.-',ms=2.0)
plot(x[1:N],u[1:N],'r.-',ms=15.0)
text(0.1,0.4,r'$u_j$',color='red',fontsize=20.0)
axis('off')
axis('tight')
axis([xL,xR,yB,yT])
figsave('vi1d.pdf')

