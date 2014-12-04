#!/usr/bin/env python

import matplotlib.font_manager as font_manager
from pylab import *
from sys import exit

x = linspace(0.0,10.0,1001)
b = 0.07*(x-3.0)**2 + 0.2*sin(2.0*x)

h0 = 3.0
L = 3.0
h = maximum(0.0, h0*(-0.2 + sqrt(maximum(0.0,1.0 - (x-5)**2/L**2))) )
s = b + h

figure(figsize=(10,4))

plot(x, b, '--k', lw=3.0)
hold(True)
plot(x, s, 'k', lw=3.0)

arrow(x[600],b[600],0.0,h[600],lw=1.0,head_width=0.1,color='k',length_includes_head=True)
arrow(x[600],s[600],0.0,-h[600],lw=1.0,head_width=0.1,color='k',length_includes_head=True)
t = text(x[600]-0.4,b[600]+0.4*h[600],'h')
font = font_manager.FontProperties(family='sans serif', style='italic', size=24)
t.set_font_properties(font)

arrow(x[450],b[600]+0.3*h[600],-1.0,-0.1,lw=2.0,head_width=0.2,color='k',length_includes_head=True)
t = text(x[400],b[600]+0.45*h[600],'q')
font = font_manager.FontProperties(family='sans serif', style='normal', weight='bold', size=24)
t.set_font_properties(font)

for j in range(10):
  xarr = x[50+100*j]
  if j>0:
    magarr = 0.6*sin(pi/2 + 0.6*xarr)
  else:
    magarr = 0.1
  arrow(xarr,s.max()+0.2,0.0,magarr,lw=1.5,head_width=0.1,color='k')
t = text(x[0]-0.5,s.max()+0.2,'f')
font = font_manager.FontProperties(family='sans serif', style='italic', size=24)
t.set_font_properties(font)

axis([0.0,10.0,-1.0,5.0])
axis('off')

savefig('layer-cartoon.pdf',bbox_inches='tight')
#show()  # debug

