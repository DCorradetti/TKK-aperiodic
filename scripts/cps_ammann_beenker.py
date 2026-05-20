"""cps_ammann_beenker.py
Article_001 — Ammann-Beenker 2D model set CPS.

Lattice: Z^4, identified with Z[sqrt 2] + i Z[sqrt 2] in C.
For lambda = (m, n, p, q) in Z^4:
  physical projection:  x = (m + n*sqrt(2),  p + q*sqrt(2)) in R^2
  internal projection: x* = (m - n*sqrt(2),  p - q*sqrt(2)) in R^2
Star map: Galois on sqrt(2) (sqrt(2) -> -sqrt(2)).

Window: regular octagon of apothem `apothem` (default 1+1/sqrt(2)).
The octagon is the intersection of 8 half-planes of the form
  cos(k*pi/4) X + sin(k*pi/4) Y  <=  apothem,    k = 0,...,7.
The regular-octagon window is D_8-invariant.

Lattice automorphisms preserving the CPS (and acting on Sigma):
  rotation by pi/4 in physical space induces the lattice automorphism
      e_0 -> e_1 -> e_2 -> e_3 -> -e_0
  and simultaneously rotates the internal projection by -3*pi/4 (because the
  Galois conjugate flips the sign of sqrt(2), inverting the rotation angle).
  We rely on the simpler description: the *internal* window must be
  rotation-by-pi/4 invariant, i.e., the regular octagon, for the rotation to
  preserve Sigma.

Version: 1
"""
import math
from typing import List, Tuple

SQRT2 = math.sqrt(2.0)


def phys(m: int, n: int, p: int, q: int) -> Tuple[float, float]:
    return (m + n * SQRT2, p + q * SQRT2)


def star(m: int, n: int, p: int, q: int) -> Tuple[float, float]:
    return (m - n * SQRT2, p - q * SQRT2)


def octagon_membership(x: float, y: float, apothem: float, eps: float = 0.0) -> bool:
    """True if (x, y) lies in the closed regular octagon of given apothem,
    centered at origin, with one face perpendicular to the +x axis.
    eps allows a tiny shrinkage to enforce a half-open boundary uniformly."""
    a = apothem - eps
    inv_sqrt2 = 1.0 / SQRT2
    return (
        x <= a and -x <= a and
        y <= a and -y <= a and
        (x + y) * inv_sqrt2 <= a and
        (x - y) * inv_sqrt2 <= a and
        (-x - y) * inv_sqrt2 <= a and
        (-x + y) * inv_sqrt2 <= a
    )


def generate(R: float, apothem: float = None, eps_boundary: float = 1e-9):
    """Return list of (m, n, p, q, x, y, xs, ys) with physical |(x, y)| <= R
    and internal (xs, ys) in the regular octagonal window."""
    if apothem is None:
        # Standard A-B apothem giving density 1 in some normalizations:
        apothem = 1.0 + 1.0 / SQRT2  # ~ 1.7071
    pts = []
    # Bound on lattice indices: |m + n sqrt 2| <= R and |p + q sqrt 2| <= R
    # implies |m|, |n*sqrt 2| <= R + R = 2R approx.
    bnd = int(R / SQRT2) + 3
    bnd_int = int(R) + 3
    for m in range(-bnd_int, bnd_int + 1):
        for n in range(-bnd, bnd + 1):
            x = m + n * SQRT2
            if abs(x) > R:
                continue
            for p in range(-bnd_int, bnd_int + 1):
                for q in range(-bnd, bnd + 1):
                    y = p + q * SQRT2
                    if abs(y) > R:
                        continue
                    if x * x + y * y > R * R:
                        continue
                    xs = m - n * SQRT2
                    ys = p - q * SQRT2
                    if octagon_membership(xs, ys, apothem, eps=eps_boundary):
                        pts.append((m, n, p, q, x, y, xs, ys))
    pts.sort(key=lambda r: (r[4], r[5]))
    return pts


def physical_norm(t):
    """Physical norm of a lattice translation t = (m, n, p, q)."""
    m, n, p, q = t
    x = m + n * SQRT2
    y = p + q * SQRT2
    return math.sqrt(x * x + y * y)


if __name__ == "__main__":
    import sys
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
    pts = generate(R)
    print(f"Ammann-Beenker patch R={R}: {len(pts)} points")
    print("First 5 points:")
    for p in pts[:5]:
        print(f"  m,n,p,q=({p[0]:3d},{p[1]:3d},{p[2]:3d},{p[3]:3d})  "
              f"phys=({p[4]:8.4f},{p[5]:8.4f})  star=({p[6]:8.4f},{p[7]:8.4f})")
    print("Last 5 points:")
    for p in pts[-5:]:
        print(f"  m,n,p,q=({p[0]:3d},{p[1]:3d},{p[2]:3d},{p[3]:3d})  "
              f"phys=({p[4]:8.4f},{p[5]:8.4f})  star=({p[6]:8.4f},{p[7]:8.4f})")
