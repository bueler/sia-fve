#!/usr/bin/env python

from pylab import *
from sys import exit

th = linspace(0.0,2.0*pi,1001)

pert = 0.6*exp(-(th-4.3)**2/0.3)
xOm = cos(th) + pert
yOm = sin(th) + pert

figure(figsize=(6,6))
plot(xOm, yOm, '-k', lw=2.0)
#hold(True)

text(xOm[150]+0.12,yOm[150]+0.12,r'$\Omega$',fontsize=20.0)

#arrow(xOm[650]+0.15,yOm[650]-0.25,-0.1,0.2)

axis('off')
#savefig('domains-fig.pdf',bbox_inches='tight')
show()  # debug

