#!/usr/bin/env python3
"""Compare the k=2,3,4 IO equations printed in PSP supplement S9 (v2, p.23,
transcribed by hand from the PDF) with -P from the paper's eq. (P)."""
import sympy as sp
a, b, g, y, y1, y2 = sp.symbols('alpha beta gamma y y1 y2')
supp = {
 2: "-2*alpha^2*gamma*y^3*y1 - alpha^2*y^3*y2 - alpha^2*y^2*y1^2 - 2*alpha*gamma*y^2*y1 - 2*alpha*y^2*y2 - beta*gamma*y^3 - beta*y^2*y1 - y*y2 + y1^2",
 3: "-3*alpha^3*gamma*y^4*y1 - alpha^3*y^4*y2 - 2*alpha^3*y^3*y1^2 - 6*alpha^2*gamma*y^3*y1 - 3*alpha^2*y^3*y2 - 3*alpha^2*y^2*y1^2 - 3*alpha*gamma*y^2*y1 - 3*alpha*y^2*y2 - beta*gamma*y^3 - beta*y^2*y1 - y*y2 + y1^2",
 4: "-4*alpha^4*gamma*y^5*y1 - alpha^4*y^5*y2 - 3*alpha^4*y^4*y1^2 - 12*alpha^3*gamma*y^4*y1 - 4*alpha^3*y^4*y2 - 8*alpha^3*y^3*y1^2 - 12*alpha^2*gamma*y^3*y1 - 6*alpha^2*y^3*y2 - 6*alpha^2*y^2*y1^2 - 4*alpha*gamma*y^2*y1 - 4*alpha*y^2*y2 - beta*gamma*y^3 - beta*y^2*y1 - y*y2 + y1^2",
}
for k, s in supp.items():
    u = 1 + a*y
    P = y*u**k*y2 - (1-(k-1)*a*y)*u**(k-1)*y1**2 + (k*a*g*y**2*u**(k-1) + b*y**2)*y1 + b*g*y**3
    d = sp.expand(sp.sympify(s.replace('^', '**'), locals={'alpha': a, 'beta': b, 'gamma': g}) + P)
    print(f'k={k}: supplement + P =', d)
