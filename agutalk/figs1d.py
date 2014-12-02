#!/usr/bin/env python

from pylab import *
from sys import exit

x = linspace(0.0,10.0,1001)
b = 0.07*(x-3.0)**2 + 0.2*sin(2.0*x)

h0 = 3.0
L = 3.0
h = maximum(0.0, h0*(-0.2 + sqrt(maximum(0.0,1.0 - (x-5)**2/L**2))) )
s = b + h

figure(figsize=(6,2))
plot(x, b, '-k', lw=2.0)
hold(True)
plot(x, s, 'k', lw=2.0)

#text(xOm[150]+0.12,yOm[150]+0.12,r'$\Omega$',fontsize=20.0)
#arrow(xOm[650]+0.15,yOm[650]-0.25,-0.1,0.2)

axis('off')
#savefig('domains-fig.pdf',bbox_inches='tight')
show()  # debug

