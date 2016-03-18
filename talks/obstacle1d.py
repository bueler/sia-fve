#!/usr/bin/env python

from pylab import *
from sys import exit

figdebug = False
def figsave(name):
    if figdebug:
        show()  # debug
    else:
        savefig(name,bbox_inches='tight')

def psi(x):
    s = 5.0 * (x - 0.5)
    s2 = s * s
    return 0.6 + 0.1 * s + 0.5 * (0.5 * s2 - 0.25 * s2 * s2)

N = 9
x = linspace(0.0,1.0,N+1)

mpsi = psi(x).max()
H = mpsi * ones(shape(x))
H[0] = 0.0
H[N] = 0.0

figure(figsize=(10,4))
hold(True)

gap = 0.01
tol = 1.0e-4
while True:
    Hold = H.copy()
    for j in range(1,N):
        H[j] = max(psi(x[j])+gap, 0.5 * (H[j-1] + H[j+1]))
    if abs(H-Hold).max() < 1.0e-4:
        break

xL = -0.1
xR = 1.1
yT = 1.0
yB = -0.2

plot([xL,xR],[0.0,0.0],'k',lw=2.0)
plot([0.0,0.0],[yB,yT],'k',lw=2.0)
xx = linspace(xL,xR,201)
plot(xx,psi(xx),'b',lw=2.0)
plot(x,H,'r.-',ms=15.0)
axis('off')
axis('tight')
axis([xL,xR,yB,yT])

figsave('obstacle1d.pdf')

