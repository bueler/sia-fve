#!/usr/bin/env python

import matplotlib.font_manager as font_manager
from pylab import *
from sys import exit

def genbasicfig(showqh=True):
    x = linspace(0.0,10.0,1001)
    b = 0.07*(x-3.0)**2 + 0.2*sin(2.0*x)

    h0 = 3.0
    L = 3.0
    h = maximum(0.0, h0*(-0.2 + sqrt(maximum(0.0,1.0 - (x-5)**2/L**2))) )
    s = b + h

    plot(x, b, '--k', lw=3.0)
    hold(True)
    plot(x, s, 'k', lw=3.0)

    if showqh:
        arrow(x[600],b[600],0.0,h[600],lw=1.0,head_width=0.1,color='k',
              length_includes_head=True)
        arrow(x[600],s[600],0.0,-h[600],lw=1.0,head_width=0.1,color='k',
              length_includes_head=True)
        t = text(x[600]-0.4,b[600]+0.4*h[600],'h')
        font = font_manager.FontProperties(family='sans serif', style='italic',
                                           size=24)
        t.set_font_properties(font)
        arrow(x[450],b[600]+0.3*h[600],-1.0,-0.1,
              lw=2.0,head_width=0.2,color='k',length_includes_head=True)
        t = text(x[400],b[600]+0.45*h[600],'q')
        font = font_manager.FontProperties(family='sans serif', style='normal',
                                           weight='bold', size=24)
        t.set_font_properties(font)

    axis([0.0,10.0,-0.5,4.5])
    axis('off')
    return x, s, b

def drawclimate(x,s,mycolor):
    for j in range(10):
        xarr = x[50+100*j]
        if j>0:
            magarr = 0.6*sin(pi/2 + 0.6*xarr)
        else:
            magarr = 0.05
        arrow(xarr,s.max()+0.2,0.0,magarr,lw=1.5,head_width=0.1,color=mycolor)

debug = False
def mysave(name):
    if debug:
        show()  # debug
    else:
        savefig(name,bbox_inches='tight')

# basic figure
figure(figsize=(10,4))
genbasicfig()
hold(False)
mysave('cartoon-layer.pdf')

figure(figsize=(10,4))
x, s, _ = genbasicfig()
drawclimate(x,s,'b')
hold(False)
mysave('cartoon-wclimate.pdf')

# for circles
th = linspace(0.0,2.0*pi+pi/100.0,201)

figure(figsize=(10,4))
x, s, _ = genbasicfig()
drawclimate(x,s,'b')
plot(x[850] + 0.4*cos(th), s.max()+0.3 + 0.4*sin(th),'-r',lw=3.0)
hold(False)
mysave('cartoon-sensitive-one.pdf')

figure(figsize=(10,4))
x, s, _ = genbasicfig()
drawclimate(x,s,'b')
plot(x[50] + 0.4*cos(th), s.max()+0.25 + 0.4*sin(th),'-r',lw=3.0)
hold(False)
mysave('cartoon-sensitive-two.pdf')

figure(figsize=(10,4))
x, s, _ = genbasicfig()
drawclimate(x,s,'b')
plot(x[795] + 0.4*cos(th), s[795] + 0.65*sin(th),'-r',lw=3.0)
hold(False)
mysave('cartoon-sensitive-three.pdf')

figure(figsize=(10,4))
x, s, b = genbasicfig(showqh=False)
h0 = 2.8
L = 3.5
h = maximum(0.0, h0*(-0.2 + sqrt(maximum(0.0,1.0 - (x-4.2)**2/L**2))) )
plot(x[h>0], b[h>0] + h[h>0], ':k', lw=3.5)
arrow(x[200]-0.4,b[200] + h[200]+0.4,0.39,-0.39,
      lw=1.5,head_width=0.1,color='k',length_includes_head=True)
font = font_manager.FontProperties(family='sans serif', style='italic', size=24)
t1 = text(x[200]-0.7,b[200] + h[200]+0.55,r'$H_{n-1}$')
t1.set_font_properties(font)
arrow(x[700]+0.2,s[700]+0.6,-0.19,-0.58,
      lw=1.5,head_width=0.1,color='k',length_includes_head=True)
t2 = text(x[700]+0.0,s[700]+0.7,r'$H_n$')
t2.set_font_properties(font)
hold(False)
mysave('cartoon-sets.pdf')

