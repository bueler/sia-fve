#!/usr/bin/env python

from pylab import *
from sys import exit

from scipy import interpolate

th = linspace(0.0,2.0*pi,1001)

pert = 0.3*exp(-(th-4.3)**2/0.3)
#pert = 0.6*exp(-(th-4.3)**2/0.3)
xOm = cos(th) + pert
yOm = sin(th) + pert

pertnx = 0.05*sin(6*th)
pertny = 0.05*sin(6*th+13.0)
xn = 0.7*cos(th) + pertnx
yn = 0.5*sin(th) + pertny

xr = xn[50:500]
yr = yn[50:500] + 0.3*exp(-(th[50:500]-1.5)**2/0.2)

figure(figsize=(6,6))
plot(xOm, yOm, '-k', lw=2.0)
hold(True)
plot(xn, yn, '-k', lw=1.5)
plot(xr, yr, '--k', lw=1.5)

myfontsize=18.0
text(0.7,0.9,r'$\Omega$',fontsize=20.0)
text(0.5,-0.6,r'$\Omega_n^0$',fontsize=myfontsize)
text(0.0,0.6,r'$\Omega_n^r$',fontsize=myfontsize)
text(-0.2,-0.1,r'$\Omega_n$',fontsize=myfontsize)

hold(False)

axis('off')
savefig('domains-fig.pdf',bbox_inches='tight')
#show()  # debug

