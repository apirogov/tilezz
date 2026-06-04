#!/usr/bin/env python3
"""
Independent CROSS-CHECK enumerator using the cyclotomic ring Z[zeta],
zeta = exp(i*pi/5), basis {1, zeta, zeta^2, zeta^3}, with reduction by
Phi_10(x) = x^4 - x^3 + x^2 - x + 1  (so zeta^4 = zeta^3 - zeta^2 + zeta - 1).

Vertices are exact elements of Z[zeta] (4 integers). This module's PURPOSE
is to dedup found polygons by EXACT geometric point-set congruence under the
full lattice symmetry group:
  - rotations: multiply every vertex by zeta^k, k=0..9 (zeta^5 = -1 etc.)
  - reflection: complex conjugation, conj(zeta) = zeta^9 = zeta^{-1}
  - translation: subtract a fixed vertex / canonicalize by sorted vertex set.
Two polygons are the same free polygon iff their vertex SETS coincide after
some lattice isometry. We compute a canonical key = lexicographically minimal
sorted tuple of exact vertex coordinates over all 20 isometries, after
translating so the set is anchored. This is fully independent of the
turn-sequence canonical used in zz10_fast.py, providing a second opinion.
"""

# element of Z[zeta] : tuple (c0,c1,c2,c3) = c0 + c1 z + c2 z^2 + c3 z^3
def zmul_by_z(a):
    """multiply by zeta once. z*(c0+c1 z+c2 z^2+c3 z^3)
       = c0 z + c1 z^2 + c2 z^3 + c3 z^4
       z^4 = z^3 - z^2 + z - 1
       => + c3*(z^3 - z^2 + z - 1)
       const: -c3 ; z: c0 + c3 ; z^2: c1 - c3 ; z^3: c2 + c3"""
    c0,c1,c2,c3 = a
    return (-c3, c0 + c3, c1 - c3, c2 + c3)

# precompute zeta^k as operators by repeated mult; represent zeta^k as the
# image of basis; easier: function rotate(point, k) applies zmul_by_z k times.
def zpow_mul(a, k):
    for _ in range(k % 10):
        a = zmul_by_z(a)
    return a

def zadd(a,b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3])
def zsub(a,b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3])
def zneg(a):
    return (-a[0],-a[1],-a[2],-a[3])

# complex conjugation: conj(z^k) = z^{-k} = z^{10-k}.
# conj(c0 + c1 z + c2 z^2 + c3 z^3) = c0 + c1 z^9 + c2 z^8 + c3 z^7.
# express z^7,z^8,z^9 in basis:
def _basis_pow(k):
    a=(1,0,0,0)
    for _ in range(k):
        a=zmul_by_z(a)
    return a
_Z = [_basis_pow(k) for k in range(10)]
def zconj(a):
    c0,c1,c2,c3 = a
    res=(c0,0,0,0)
    # + c1 * z^9 + c2 * z^8 + c3 * z^7
    r=list(res)
    for coef,k in ((c1,9),(c2,8),(c3,7)):
        b=_Z[k]
        r[0]+=coef*b[0]; r[1]+=coef*b[1]; r[2]+=coef*b[2]; r[3]+=coef*b[3]
    return tuple(r)

# step vector in facing direction d = zeta^d = _Z[d]
STEP = [_Z[d] for d in range(10)]
ORIG = (0,0,0,0)

def vertices_from_facings(facings):
    pts=[ORIG]; cur=ORIG
    for d in facings:
        cur=zadd(cur, STEP[d])
        pts.append(cur)
    return pts[:-1]   # drop the closing dup of origin

def canon_pointset(facings):
    pts = vertices_from_facings(facings)
    best=None
    # 20 isometries: rot k (0..9) optionally composed with conjugation.
    for conj in (False, True):
        base = [zconj(p) if conj else p for p in pts]
        for k in range(10):
            trans=[zpow_mul(p,k) for p in base]
            # translate so that the lexicographically minimal vertex -> origin
            mn=min(trans)
            shifted=tuple(sorted(zsub(p,mn) for p in trans))
            if best is None or shifted<best:
                best=shifted
    return best
