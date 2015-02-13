#!/usr/bin/env python

from sympy import *

init_printing(use_unicode=True)

s, t, p = symbols('s t p')

f = (1 - (t**(p-1)+t)*s + t**p) / (1 - 2*s*t + t**2)**(p/2)

dfds = diff(f,s)

y = dfds / (t * (1 - 2*s*t + t**2)**(-1-p/2))

z = powsimp( collect(simplify(y), s) )

print z

print collect(z.coeff(s), t)

print collect(simplify(z - s*z.coeff(s)), t)

