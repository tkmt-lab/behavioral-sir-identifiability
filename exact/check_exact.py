#!/usr/bin/env python3
"""Exact elimination from the vector field; imports no other project code.

Variables are (y,v=y',w=y'',alpha,beta,gamma,S). Sparse integer
polynomials are sufficient; Fraction is used only for exact specialization.
The source of the IO polynomial here is the resultant of F and D(F),
not a transcription of the proposed IO polynomial. Finite checks only.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

N = 7
ONE = {(0,) * N: 1}


def var(i):
    e = [0] * N
    e[i] = 1
    return {tuple(e): 1}


def plus(*ps):
    out = {}
    for p in ps:
        for m, c in p.items():
            out[m] = out.get(m, 0) + c
    return {m: c for m, c in out.items() if c}


def scale(p, n):
    return {m: c * n for m, c in p.items() if c * n}


def times(*ps):
    out = ONE
    for p in ps:
        new = {}
        for m, c in out.items():
            for n, d in p.items():
                e = tuple(a + b for a, b in zip(m, n))
                new[e] = new.get(e, 0) + c * d
        out = {m: c for m, c in new.items() if c}
    return out


def power(p, n):
    assert n >= 0
    out = ONE
    for _ in range(n):
        out = times(out, p)
    return out


def diff(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            e = list(m)
            e[i] -= 1
            out[tuple(e)] = c * m[i]
    return out


def divide_exact(p, q):
    """One-divisor multivariate division over Z with a unit lead coefficient."""
    out, rest = {}, dict(p)
    qm = max(q)
    qc = q[qm]
    assert abs(qc) == 1
    while rest:
        rm = max(rest)
        e = tuple(a - b for a, b in zip(rm, qm))
        assert min(e) >= 0, (rm, qm)
        t = {e: rest[rm] // qc}
        out = plus(out, t)
        rest = plus(rest, scale(times(t, q), -1))
    assert times(out, q) == p
    return out


y, v, w, a, b, g, S = [var(i) for i in range(N)]
u = plus(ONE, times(a, y))


def eliminated(k):
    # F = u^k(v+g*y)-b*S*y = F0+F1*S, and D(F)=0 on the model.
    uk = power(u, k)
    f0 = times(uk, plus(v, times(g, y)))
    f1 = scale(times(b, y), -1)
    # G=u^k D(F), with D(S)=-b*S*y/u^k, D(y)=v, D(v)=w.
    g0 = times(uk, plus(times(diff(f0, 0), v), times(diff(f0, 1), w)))
    g1 = plus(times(uk, diff(f1, 0), v), scale(times(f1, b, y), -1))
    resultant = plus(times(f1, g0), scale(times(g1, f0), -1))
    return divide_exact(resultant, scale(times(b, uk), -1))


def coeff(p, ey, ev, ew):
    return {m: c for m, c in p.items() if m[:3] == (ey, ev, ew)}


def parameter_coeff(p, ey, ev, ew):
    out = {}
    for m, c in coeff(p, ey, ev, ew).items():
        assert m[6] == 0
        n = (0, 0, 0) + m[3:]
        out[n] = c
    return out


def choose(n, j):
    return math.comb(n, j) if 0 <= j <= n else 0


def expected(k):
    # Separate closed formulas, used only on the comparison side.
    out = {}
    def put(ey, ev, ew, ea, eb, eg, c):
        nonlocal out
        out = plus(out, {(ey, ev, ew, ea, eb, eg, 0): c})
    for j in range(k + 1):
        put(j + 1, 0, 1, j, 0, 0, choose(k, j))
        put(j, 2, 0, j, 0, 0,
            -(choose(k - 1, j) - (k - 1) * choose(k - 1, j - 1)))
    for j in range(1, k + 1):
        put(j + 1, 1, 0, j, 0, 1, k * choose(k - 1, j - 1))
    put(2, 1, 0, 0, 1, 0, 1)
    put(3, 0, 0, 0, 1, 1, 1)
    return out


def substitute(p, i, replacement):
    out = {}
    powers = {e: power(replacement, e) for e in {m[i] for m in p}}
    for m, c in p.items():
        e = list(m)
        e[i] = 0
        out = plus(out, times({tuple(e): c}, powers[m[i]]))
    return out


def exact_vector_field_check(p, k):
    # Independently differentiate v=H/u^k along the original S,I vector field.
    H = plus(times(b, S, y), scale(times(g, y, power(u, k)), -1))
    J = plus(times(u, plus(scale(times(b, S, y, diff(H, 6)), -1),
                               times(H, diff(H, 0)))),
             scale(times(a, H, H), -k))
    # v=H/u^k, w=J/u^(2k+1). Multiply the whole substitution by u^(2k+1).
    out = {}
    for m, c in p.items():
        assert m[1] <= 2 and m[2] <= 1
        e = list(m)
        ev, ew = e[1], e[2]
        e[1] = e[2] = 0
        exponent = 2 * k + 1 - k * ev - (2 * k + 1) * ew
        assert exponent >= 0
        out = plus(out, times({tuple(e): c}, power(H, ev), power(J, ew),
                              power(u, exponent)))
    assert out == {}


def at_y_minus_inv_a(p):
    # Polynomial in the other variables, allowing negative powers of a.
    out = {}
    for m, c in p.items():
        e = list(m)
        ey = e[0]
        e[0] = 0
        e[3] -= ey
        out = plus(out, {tuple(e): c * (-1) ** ey})
    return out


def specialize_parameters(p, av, bv, gv):
    out = {}
    for m, c in p.items():
        e = list(m)
        ea, eb, eg = e[3:6]
        e[3:6] = [0, 0, 0]
        out = plus(out, {tuple(e): c * av ** ea * bv ** eb * gv ** eg})
    return out


def hash_poly(p):
    raw = json.dumps(sorted(p.items()), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def printed_s9():
    # New transcription from the fetched supplement, p23 (code blocks).
    # Each tuple is (integer coefficient, y,v,w,a,b,g exponents).
    rows = {
      2: [(-2,3,1,0,2,0,1),(-1,3,0,1,2,0,0),(-1,2,2,0,2,0,0),
          (-2,2,1,0,1,0,1),(-2,2,0,1,1,0,0),(-1,3,0,0,0,1,1),
          (-1,2,1,0,0,1,0),(-1,1,0,1,0,0,0),(1,0,2,0,0,0,0)],
      3: [(-3,4,1,0,3,0,1),(-1,4,0,1,3,0,0),(-2,3,2,0,3,0,0),
          (-6,3,1,0,2,0,1),(-3,3,0,1,2,0,0),(-3,2,2,0,2,0,0),
          (-3,2,1,0,1,0,1),(-3,2,0,1,1,0,0),(-1,3,0,0,0,1,1),
          (-1,2,1,0,0,1,0),(-1,1,0,1,0,0,0),(1,0,2,0,0,0,0)],
      4: [(-4,5,1,0,4,0,1),(-1,5,0,1,4,0,0),(-3,4,2,0,4,0,0),
          (-12,4,1,0,3,0,1),(-4,4,0,1,3,0,0),(-8,3,2,0,3,0,0),
          (-12,3,1,0,2,0,1),(-6,3,0,1,2,0,0),(-6,2,2,0,2,0,0),
          (-4,2,1,0,1,0,1),(-4,2,0,1,1,0,0),(-1,3,0,0,0,1,1),
          (-1,2,1,0,0,1,0),(-1,1,0,1,0,0,0),(1,0,2,0,0,0,0)]}
    return {k: {tuple(row[1:]) + (0,): row[0] for row in rs}
            for k, rs in rows.items()}


def main():
    source = printed_s9()
    rows = []
    for k in range(1, 13):
        p = eliminated(k)
        assert p == expected(k), ('closed coefficients', k)
        exact_vector_field_check(p, k)
        assert parameter_coeff(p, 1, 0, 1) == ONE
        assert parameter_coeff(p, 2, 0, 1) == scale(a, k)
        assert parameter_coeff(p, 2, 1, 0) == plus(b, scale(times(a, g), k))
        assert parameter_coeff(p, 3, 0, 0) == times(b, g)
        if k >= 2:
            assert parameter_coeff(p, 3, 1, 0) == scale(times(a, a, g), k * (k - 1))
        # Exact generic recovery is the above identities plus divisions by k,k-1,a^2.
        # Also exercise rational specialization, without floating point.
        av, bv, gv = Fraction(2), Fraction(3), Fraction(5)
        same_alt = (specialize_parameters(p, av, bv, gv) ==
                    specialize_parameters(p, av, av * gv, bv / av))
        assert same_alt == (k == 1)
        # Primitive in w: only possible common factors are y and u.
        q = {m: c for m, c in p.items() if m[2] == 0}
        q_at_zero = {m: c for m, c in q.items() if m[0] == 0}
        assert q_at_zero == scale(times(v, v), -1)
        q_at_pole = at_y_minus_inv_a(q)
        if k >= 2:
            assert q_at_pole == {(0,1,0,-2,1,0,0): 1,
                                 (0,0,0,-3,1,1,0): -1}
        else:
            assert q_at_pole  # k=1 is still primitive but has an extra -v^2.
        if k in source:
            assert plus(p, source[k]) == {}, ('source S9', k)
        # Verify scaled derivation: (u^k/y)D(S)=-bS;
        # (u^k/y)D(u)=a*b*S-a*g*u^k, as numerator polynomial identities.
        H = plus(times(b, S, y), scale(times(g, y, power(u, k)), -1))
        assert times(a, H) == times(y, plus(times(a, b, S),
                                                      scale(times(a, g, power(u, k)), -1)))
        row = {'k': k, 'expanded_terms': len(p), 'polynomial_sha256': hash_poly(p),
               'resultant_equals_closed': True, 'vector_field_identity': True,
               'recovery_coefficient_identities': True, 'primitive_checks': True,
               'same_alt_parameters': same_alt,
               'source_s9_equal_up_to_sign': True if k in source else None}
        rows.append(row)
        print(json.dumps(row, sort_keys=True))
    print(json.dumps({'status': 'PASS', 'k_range': [1, 12],
                      'arithmetic': 'Z multivariate polynomials and Fraction',
                      'own_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope': 'finite algebra checks; not a check of the general Laurent or OPT proofs'}))


if __name__ == '__main__':
    main()
