#!/usr/bin/env python3
"""Recomputation by a separate reader from the article text (draft version 1; same equation numbers as 1.1).

Written from the article text only (not from exact/check_exact.py or the proof notes).
Run: python3 -B independent/check_article_algebra.py   (needs sympy)
"""
import sympy as sp
from sympy import binomial, Rational as Q

a, b, g, S, I, y, v, w, k = sp.symbols('alpha beta gamma S I y v w k')


def P_of(kk):
    u = 1 + a * y
    return sp.expand(y * u**kk * w - (1 - (kk - 1) * a * y) * u**(kk - 1) * v**2
                     + (kk * a * g * y**2 * u**(kk - 1) + b * y**2) * v + b * g * y**3)


def model_D(kk):
    u = 1 + a * I
    dS = -b * S * I / u**kk
    dI = b * S * I / u**kk - g * I
    return lambda F: sp.diff(F, S) * dS + sp.diff(F, I) * dI


ok = True

# 1. P vanishes on the model, symbolic k (u^k treated as a genuine power with symbolic k)
D = model_D(k)
Y = I
V = D(Y)
W = D(V)
u = 1 + a * I
Psym = (I * u**k * W - (1 - (k - 1) * a * I) * u**(k - 1) * V**2
        + (k * a * g * I**2 * u**(k - 1) + b * I**2) * V + b * g * I**3)
# sympy does not normalise (1+aI)^(k-1)*(1+aI)^(-2k) for symbolic k; test symbolic-k
# expression at exact rational points with several k (incl. non-integer k, positive base),
# and exactly for integer k = 1..8.
import random
random.seed(1)
r_ok = True
for kk in [Q(5, 2), 7, Q(13, 3)]:
    for _ in range(3):
        pt0 = {a: Q(random.randint(1, 9), 7), b: Q(random.randint(1, 9), 5), g: Q(random.randint(1, 9), 11),
               S: Q(random.randint(1, 9), 3), I: Q(random.randint(1, 9), 4), k: kk}
        val = sp.nsimplify(Psym.subs(pt0))
        val = sp.N(val, 50)
        r_ok &= abs(val) < 1e-40
for kk in range(1, 9):
    Dk = model_D(kk)
    Vk = Dk(I); Wk = Dk(Vk); uk = 1 + a * I
    e = (I * uk**kk * Wk - (1 - (kk - 1) * a * I) * uk**(kk - 1) * Vk**2
         + (kk * a * g * I**2 * uk**(kk - 1) + b * I**2) * Vk + b * g * I**3)
    r_ok &= sp.cancel(e) == 0
print('[1] P vanishes on model (exact k=1..8; symbolic-k expr at rational points):', r_ok)
ok &= r_ok

# 1b. derive P independently: S from y,v; D(S formula) + (v+gamma y) = 0, times beta y^2
#     use generic k with symbol U = u^k, DU = k alpha U v / u
U = sp.Symbol('U')
uu = 1 + a * y
Sform = (v + g * y) * U / (b * y)
# total derivative along jets: Dy=v, Dv=w, DU = k*alpha*U*v/uu
DS = sp.diff(Sform, y) * v + sp.diff(Sform, v) * w + sp.diff(Sform, U) * k * a * U * v / uu
expr = sp.together((DS + (v + g * y)) * b * y**2 * uu / U)  # multiply by uu/U to get polynomial times u^{k-1}
Pk = (y * U * w - (1 - (k - 1) * a * y) * U / uu * v**2 + (k * a * g * y**2 * U / uu + b * y**2) * v + b * g * y**3)
diff = sp.simplify(sp.together((DS + (v + g * y)) * b * y**2 - Pk))
print('[1b] derived relation minus P (symbolic k):', diff)
ok &= (diff == 0)

# 2. Coefficient table, integer k = 1..12
for kk in range(1, 13):
    P = sp.Poly(P_of(kk), y, v, w)
    table = {}
    for j in range(0, kk + 1):
        table[(j + 1, 0, 1)] = binomial(kk, j) * a**j
    for j in range(1, kk + 1):
        table[(j + 1, 1, 0)] = table.get((j + 1, 1, 0), 0) + kk * binomial(kk - 1, j - 1) * g * a**j + (b if j == 1 else 0)
    for j in range(0, kk + 1):
        c = -(binomial(kk - 1, j) - (kk - 1) * (binomial(kk - 1, j - 1) if j >= 1 else 0)) * a**j
        table[(j, 2, 0)] = c
    table[(3, 0, 0)] = table.get((3, 0, 0), 0) + b * g
    got = dict(P.terms())
    keys = set(got) | set(table)
    bad = [m for m in keys if sp.expand(got.get(m, 0) - table.get(m, 0)) != 0]
    four = (got[(1, 0, 1)], got[(2, 0, 1)], got[(2, 1, 0)], got.get((3, 1, 0), 0))
    inv_ok = True
    if kk >= 2:
        al = four[1] / kk
        ga = four[3] / (kk * (kk - 1) * al**2)
        be = four[2] - kk * al * ga
        inv_ok = all(sp.simplify(x - t) == 0 for x, t in [(al, a), (ga, g), (be, b)])
    print(f'[2] k={kk}: table mismatches={bad}, four={four}, inverse ok={inv_ok}, [y v^2]={got.get((1,2,0),0)}')
    ok &= (not bad) and inv_ok and got.get((1, 2, 0), 0) == 0

# 2b. table for symbolic k via binomial expansion check at j=0..3 coefficients (generic formula)
# compare closed formula with series coefficient for symbolic k using sympy series in y
for kk in [k]:
    Pexpr = P_of(kk)
    # not polynomial for symbolic k; check low-order coefficients via taylor
    ser = sp.series(Pexpr.subs({w: sp.Symbol('W1')}), y, 0, 5).removeO()
print('[2b] (symbolic-k series skipped; integer table covers 1..12)')

# 3. Irreducibility over Q(alpha,beta,gamma) for k=2..6
for kk in range(2, 7):
    f = sp.factor_list(P_of(kk), y, v, w)
    nontriv = [(fac, e) for fac, e in f[1] if sp.Poly(fac, y, v, w).total_degree() > 0]
    print(f'[3] k={kk}: content={f[0]}, factors over Q[a,b,g][y,v,w]: {len(f[1])} -> nonconstant in (y,v,w): {len(nontriv)}')
    ok &= len(nontriv) == 1 and nontriv[0][1] == 1
    # B(0,v) and B(-1/alpha,v)
    Pp = P_of(kk)
    B = sp.expand(Pp - sp.Poly(Pp, w).coeff_monomial(w) * w)
    b0 = sp.simplify(B.subs(y, 0))
    bm = sp.factor(B.subs(y, -1 / a))
    print(f'     B(0,v)={b0}, B(-1/a,v)={bm}')
    ok &= sp.simplify(bm - b / a**2 * (v - g / a)) == 0 and sp.simplify(b0 + v**2) == 0

# 4. derivation E in (S,u)
uS = sp.Symbol('u')
for kk in [2, 3, 5, k]:
    Dk = model_D(kk)
    uI = 1 + a * I
    ES = sp.simplify(uI**kk / I * Dk(S))
    Eu = sp.simplify(uI**kk / I * Dk(uI))
    c1 = sp.simplify(ES + b * S)
    c2 = sp.simplify(Eu - (a * b * S - a * g * uI**kk))
    print(f'[4] k={kk}: ES+beta S={c1}, Eu-(..)={c2}')
    ok &= c1 == 0 and c2 == 0

# 5. k=1 swap invariance and k>=2 non-invariance
sw = {b: a * g, g: b / a}
for kk in range(1, 8):
    P = sp.Poly(P_of(kk), y, v, w)
    P2 = sp.Poly(P_of(kk).subs(sw, simultaneous=True), y, v, w)
    same = all(sp.simplify(c1 - P2.coeff_monomial(m)) == 0 for m, c1 in P.terms()) and \
        all(sp.simplify(c2 - P.coeff_monomial(m)) == 0 for m, c2 in P2.terms())
    print(f'[5] k={kk}: coefficients invariant under (beta,gamma)->(alpha gamma, beta/alpha): {same}')
    ok &= same == (kk == 1)

# 6. single trajectory, exact power series, k=2,3
t = sp.Symbol('t')
pt = {a: Q(1, 3), b: Q(2, 5), g: Q(1, 7)}
S0, I0 = Q(5, 2), Q(1, 2)
N = 26
for kk in [2, 3]:
    # Picard-type recursion on truncated series with exact rationals
    Sser = [S0] + [0] * N
    Iser = [I0] + [0] * N

    def mul(p, q):
        out = [0] * (N + 1)
        for i, pi in enumerate(p):
            if pi == 0:
                continue
            for j in range(0, N + 1 - i):
                out[i + j] += pi * q[j]
        return out

    def inv(p):
        out = [0] * (N + 1)
        out[0] = 1 / sp.Rational(p[0])
        for n in range(1, N + 1):
            out[n] = -sum(p[i] * out[n - i] for i in range(1, n + 1)) / p[0]
        return out

    for n in range(N):
        user = [1 + pt[a] * Iser[0]] + [pt[a] * c for c in Iser[1:]]
        upk = [1] + [0] * N
        for _ in range(kk):
            upk = mul(upk, user)
        rate = mul(mul(Sser, Iser), inv(upk))
        rate = [pt[b] * c for c in rate]
        Sser[n + 1] = -rate[n] / (n + 1)
        Iser[n + 1] = (rate[n] - pt[g] * Iser[n]) / (n + 1)
    ypoly = sum(c * t**i for i, c in enumerate(Iser))
    yp = sp.diff(ypoly, t)
    ypp = sp.diff(yp, t)
    mons = [(i, j, l) for i in range(kk + 2) for j in range(3) for l in range(2) if j + 2 * l <= 2]
    cs = sp.symbols(f'c0:{len(mons)}')
    rel = sum(c * ypoly**i * yp**j * ypp**l for c, (i, j, l) in zip(cs, mons))
    relp = sp.Poly(sp.expand(rel), t)
    eqs = [relp.coeff_monomial(t**n) for n in range(N - 4)]
    M = sp.Matrix([[sp.diff(e, c) for c in cs] for e in eqs])
    ns = M.nullspace()
    vec = ns[0] / ns[0][mons.index((1, 0, 1))]
    coef = dict(zip(mons, vec))
    al = coef[(2, 0, 1)] / kk
    ga = coef[(3, 1, 0)] / (kk * (kk - 1) * al**2)
    be = coef[(2, 1, 0)] - kk * al * ga
    Srec = (Iser[1] + ga * I0) * (1 + al * I0)**kk / (be * I0)
    print(f'[6] k={kk}: nullity={len(ns)}, recovered (a,b,g)=({al},{be},{ga}), S0={Srec}')
    ok &= len(ns) == 1 and (al, be, ga, Srec) == (pt[a], pt[b], pt[g], S0)

print('ALL OK' if ok else 'SOME CHECK FAILED')
