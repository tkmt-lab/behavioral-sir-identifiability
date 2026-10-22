"""Checks of the written proof by a separate reader, written without reading exact/check_exact.py
or its results. Equation numbers in labels and comments refer to the proof notes (not distributed):
(3), (4) = article Eq. (4) before and after collecting terms; (6) = the values B(0,v), B(-1/alpha,v)
of Lemma 2 (named Q here); (7) = article Eq. (5); (9) = Eq. (7); (10) = Eq. (8); (12) = Eq. (9).
Run: python3 -B independent/check_proof_algebra.py   (needs sympy)"""
import sys, random
import sympy as sp

a, b, g, S, I, y, v, w, k = sp.symbols('alpha beta gamma S I y v w k')
ok = True
def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)

def P3(kk):  # eq (3)
    u = 1 + a*y
    return y*u**kk*w + g*y*v*u**kk - (v+g*y)*v*(1+(1-kk)*a*y)*u**(kk-1) + b*(v+g*y)*y**2
def P4(kk):  # eq (4)
    u = 1 + a*y
    return y*u**kk*w - (1-(kk-1)*a*y)*u**(kk-1)*v**2 + (kk*a*g*y**2*u**(kk-1) + b*y**2)*v + b*g*y**3

def model_jets(kk):
    u = 1 + a*I
    dS = -b*S*I/u**kk
    dI = b*S*I/u**kk - g*I
    D = lambda F: sp.diff(F, S)*dS + sp.diff(F, I)*dI
    return I, dI, D(dI)

# 1. P vanishes on the model, (3)==(4), for k=2,3,5; symbolic k via U=(1+a*I)^k as a separate symbol
for kk in [2, 3, 5]:
    Y, V, W = model_jets(kk)
    for lab, P in [('(3)', P3(kk)), ('(4)', P4(kk))]:
        val = sp.cancel(sp.together(P.subs({y: Y, v: V, w: W}, simultaneous=True)))
        check(f'P{lab} vanishes on model, k={kk}', val == 0)
    check(f'(3)-(4)==0, k={kk}', sp.expand(P3(kk) - P4(kk)) == 0)
U = sp.symbols('U')   # U stands for (1+a*I)^k, k symbolic; D U = k*a*U/(1+a*I) * D I
uI = 1 + a*I
dS = -b*S*I/U; dI = b*S*I/U - g*I
Dk = lambda F: sp.diff(F, S)*dS + sp.diff(F, I)*dI + sp.diff(F, U)*k*a*U/uI*dI
Y, V = I, dI; W = Dk(dI)
for lab, P in [('(3)', lambda: y*U*w + g*y*v*U - (v+g*y)*v*(1+(1-k)*a*y)*U/uI + b*(v+g*y)*y**2),
               ('(4)', lambda: y*U*w - (1-(k-1)*a*y)*U/uI*v**2 + (k*a*g*y**2*U/uI + b*y**2)*v + b*g*y**3)]:
    expr = P().subs({y: Y, v: V, w: W}, simultaneous=True)
    check(f'P{lab} vanishes on model, symbolic k', sp.cancel(sp.together(expr)) == 0)

# 2. coefficient table (7) for k=2..9; inverse map (9)
B = lambda n, r: sp.binomial(n, r) if 0 <= r <= n else 0
for kk in range(2, 10):
    Pp = sp.Poly(sp.expand(P4(kk)), y, v, w)
    pred = {}
    for j in range(0, kk+1):
        pred[(j+1, 0, 1)] = B(kk, j)*a**j
        pred[(j, 2, 0)] = sp.expand(-(B(kk-1, j) - (kk-1)*B(kk-1, j-1))*a**j)
    for j in range(1, kk+1):
        pred[(j+1, 1, 0)] = kk*B(kk-1, j-1)*g*a**j + (b if j == 1 else 0)
    pred[(3, 0, 0)] = b*g
    pred = {m: c for m, c in pred.items() if c != 0}
    got = {m: sp.expand(c) for m, c in Pp.terms()}
    check(f'(7) table exact and exhaustive, k={kk}', got == {m: sp.expand(c) for m, c in pred.items()})
    c21, d2, d3 = got[(2, 0, 1)], got[(2, 1, 0)], got[(3, 1, 0)]
    al = c21/kk; ga = d3/(kk*(kk-1)*al**2); be = d2 - kk*al*ga
    check(f'(9) inverse map, k={kk}', sp.simplify(al-a) == 0 and sp.simplify(ga-g) == 0 and sp.simplify(be-b) == 0)

# 3. content / irreducibility checks, eq (6)
for kk in [2, 3, 5]:
    P = sp.expand(P4(kk)); A = y*(1+a*y)**kk; Q = sp.expand(P - A*w)
    check(f'Q(0,v)=-v^2, k={kk}', sp.expand(Q.subs(y, 0)) == -v**2)
    check(f'Q(-1/a,v), k={kk}', sp.simplify(Q.subs(y, -1/a) - b/a**2*(v-g/a)) == 0)
    fl = sp.factor_list(P)
    check(f'P irreducible in Q[a,b,g,y,v,w], k={kk}', len(fl[1]) == 1 and fl[1][0][1] == 1)
# k=1 counterpart: (6) second formula fails / d3 vanishes
P1 = sp.Poly(sp.expand(P4(1)), y, v, w)
check('k=1: [y^3 v]P = 0 (d3 vanishes)', P1.coeff_monomial(y**3*v) == 0)

# 4. Lemma 5 of the article: ingredients (10), (12)
u = sp.symbols('u')
Ival = (u-1)/a
dS = -b*S*Ival/u**2  # placeholder overwritten below
for kk in [2, 3, 5, k]:
    dS = -b*S*Ival/u**kk; dI = b*S*Ival/u**kk - g*Ival
    du = a*dI
    scale = u**kk/Ival
    check(f'(10) ES=-bS, k={kk}', sp.simplify(scale*dS + b*S) == 0)
    check(f'(10) Eu, k={kk}', sp.simplify(scale*du - (a*b*S - a*g*u**kk)) == 0)
    m = sp.symbols('m'); f = sp.Function('f')(u)
    E = lambda F: -b*S*sp.diff(F, S) + (a*b*S - a*g*u**kk)*sp.diff(F, u)
    lead = sp.expand(E(f*S**m) / S**m)
    lead0 = sp.simplify(lead.subs(S, 0))
    check(f'(12) leading coefficient, k={kk}', sp.simplify(lead0 - (-b*m*f - a*g*u**kk*sp.diff(f, u))) == 0)

# 5. spot check: no nonconstant polynomial first integral of degree <=6 in (S,u), k=2,3, random numeric params
random.seed(1)
for kk in [2, 3]:
    pa, pb, pg = [sp.Rational(random.randint(1, 9), random.randint(1, 9)) for _ in range(3)]
    deg = 6
    mons = [S**i*u**j for i in range(deg+1) for j in range(deg+1-i) if i+j > 0]
    cs = sp.symbols(f'c0:{len(mons)}')
    F = sum(c*mm for c, mm in zip(cs, mons))
    EF = sp.expand(-pb*S*sp.diff(F, S) + (pa*pb*S - pa*pg*u**kk)*sp.diff(F, u))
    eqs = sp.Poly(EF, S, u).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    check(f'no poly first integral deg<=6, k={kk}', sol == [{c: 0 for c in cs}] or all(vv == 0 for vv in sol[0].values()) and len(sol[0]) == len(cs))

# 6. single-trajectory check (exact rationals): Taylor series of one trajectory; all polynomials
#    sum c_{ijl} y^i (y')^j (y'')^l over the box i<=k+1, j+2l<=2 vanishing on it form a 1-dim space.
#    Since any competitor (a',b',g') with the same k gives a P' in this box vanishing on the same y,
#    and [y y'']P'=1, this forces P'=P, hence (a',b',g')=(a,b,g) by (9), for this trajectory.
from fractions import Fraction as Fr
def mul(p, q, N): 
    r = [Fr(0)]*N
    for i, x in enumerate(p):
        if x:
            for j in range(N-i): r[i+j] += x*q[j]
    return r
def inv(p, N):
    r = [Fr(0)]*N; r[0] = 1/p[0]
    for n in range(1, N): r[n] = -sum(p[i]*r[n-i] for i in range(1, n+1))/p[0]
    return r
def deriv(p): return [i*p[i] for i in range(1, len(p))] + [Fr(0)]
def solve_series(al, be, ga, S0, I0, kk, N):
    Ss = [Fr(S0)]+[Fr(0)]*(N-1); Is = [Fr(I0)]+[Fr(0)]*(N-1)
    for _ in range(N):
        uu = [1+al*Is[0]]+[al*x for x in Is[1:]]
        uk = [Fr(1)]+[Fr(0)]*(N-1)
        for _ in range(kk): uk = mul(uk, uu, N)
        inc = mul(mul(Ss, Is, N), inv(uk, N), N); inc = [be*x for x in inc]
        dS = [-x for x in inc]; dI = [inc[i]-ga*Is[i] for i in range(N)]
        Ss = [Fr(S0)]+[dS[i]/(i+1) for i in range(N-1)]
        Is = [Fr(I0)]+[dI[i]/(i+1) for i in range(N-1)]
    return Is
for kk in [2, 3, 5]:
    al, be, ga, S0, I0 = Fr(1, 3), Fr(2, 5), Fr(1, 7), Fr(5, 2), Fr(1, 2)
    N = 40
    Ys = solve_series(al, be, ga, S0, I0, kk, N); Vs = deriv(Ys); Ws = deriv(Vs)
    M = 30  # trust first M series coefficients
    pw = lambda p, e: (lambda r: r)(__import__('functools').reduce(lambda x, _: mul(x, p, N), range(e), [Fr(1)]+[Fr(0)]*(N-1)))
    mons = [(i, j, l) for i in range(kk+2) for j in range(3) for l in range(2) if j+2*l <= 2]
    cols = [mul(mul(pw(Ys, i), pw(Vs, j), N), pw(Ws, l), N)[:M] for (i, j, l) in mons]
    Mat = sp.Matrix(M, len(mons), lambda r, c: sp.Rational(cols[c][r].numerator, cols[c][r].denominator))
    ns = Mat.nullspace()
    good = len(ns) == 1
    if good:
        vec = ns[0] / ns[0][mons.index((1, 0, 1))]
        coef = dict(zip(mons, vec))
        alr = coef[(2, 0, 1)]/kk; gar = coef[(3, 1, 0)]/(kk*(kk-1)*alr**2); ber = coef[(2, 1, 0)] - kk*alr*gar
        good = (alr, ber, gar) == (sp.Rational(1, 3), sp.Rational(2, 5), sp.Rational(1, 7))
    check(f'one trajectory: nullspace dim {len(ns)}, params recovered exactly, k={kk}', good)

print('ALL PASS' if ok else 'SOME FAIL')
sys.exit(0 if ok else 1)
