#!/usr/bin/env python

from pylab import *
from sys import exit

from scipy import interpolate

th = linspace(0.0,2.0*pi,1001)

pert = 0.6*exp(-(th-4.3)**2/0.3)
xOm = cos(th) + pert
yOm = sin(th) + pert

xs = xOm[600]
ys = yOm[600]
xe = xOm[800]
ye = yOm[800]
xmid = 0.0
ymid = 0.5
xRe = xOm[850]
yRe = yOm[850]
Gam0x = array([xs, -0.5, xmid, 0.4, xe])
Gam0y = array([ys, 0.2, ymid, 0.0, ye])
GamRx = array([xmid, 0.5, 0.8, 0.8, xRe])
GamRy = array([ymid, 0.35, 0.0, -0.3, yRe])

t = np.linspace(0.0, 4.0, 5)
tfine = np.linspace(0.0, 4.0, 1000)
spl0x = interpolate.UnivariateSpline(t,Gam0x,s=0)
spl0y = interpolate.UnivariateSpline(t,Gam0y,s=0)
splRx = interpolate.UnivariateSpline(t,GamRx,s=0)
splRy = interpolate.UnivariateSpline(t,GamRy,s=0)

figure(figsize=(6,6))
plot(xOm, yOm, '-k', lw=2.0)
plot(xOm[600:800], yOm[600:800], '-k', lw=4.5)
hold(True)
#plot(Gam0x, Gam0y, 'ok')  # debug
plot(spl0x(tfine), spl0y(tfine), '-k', lw=1.5)
#plot(GamRx, GamRy, 'dk')  # debug
plot(splRx(tfine), splRy(tfine), '--k', lw=1.5)

myfontsize=18.0
text(xOm[150]+0.12,yOm[150]+0.12,r'$\Omega$',fontsize=20.0)
text(xOm[250],yOm[250]-0.3,r'$\Omega_n^0$',fontsize=myfontsize)
text(GamRx[2]-0.2,GamRy[2]-0.2,r'$\Omega_n^r$',fontsize=myfontsize)
text(xs+0.3,ys+0.3,r'$\Omega_n$',fontsize=myfontsize)

arrow(xOm[650]+0.15,yOm[650]-0.25,-0.1,0.2)
text(xOm[650]+0.12,yOm[650]-0.4,r'$\Gamma_n^N$',fontsize=myfontsize)
arrow(Gam0x[1]-0.08,Gam0y[1]+0.08,0.05,-0.05)
text(Gam0x[1]-0.2,Gam0y[1]+0.15,r'$\Gamma_n^0$',fontsize=myfontsize)
hold(False)

axis('off')
savefig('domains-fig.pdf',bbox_inches='tight')
#show()  # debug

