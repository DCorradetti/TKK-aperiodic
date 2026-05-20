"""cps_fibonacci.py
Article_001 — G2: Fibonacci cut-and-project model set.

We use the standard Fibonacci CPS:
    Lambda = Z[phi] (rank-2 lattice in R x R)
    physical projection: (a, b) |-> a + b * phi
    internal  projection: (a, b) |-> a + b * phi'   (algebraic conjugate)
    acceptance window Omega = [-1, 1/phi)            (standard half-open choice)

The model set Sigma is then
    Sigma = { a + b*phi : (a, b) in Z^2, a + b*phi' in Omega }.

This is the standard Fibonacci chain, with two gaps S = 1/phi and L = 1
between consecutive points and density 1/phi per unit length.

Version: 1
"""
from typing import List, Tuple
from fields import PHI, PHI_STAR, INV_PHI


def generate(R: float, window: Tuple[float, float] = (-1.0, INV_PHI)) -> List[Tuple[int, int, float, float]]:
    """Return sorted list of (a, b, x, x_star) for the Fibonacci patch with
    -R <= x <= R and window[0] <= x_star < window[1]."""
    pts = []
    a_max = int(R) + 5
    b_max = int(R / PHI) + 5
    for a in range(-a_max, a_max + 1):
        for b in range(-b_max, b_max + 1):
            x = a + b * PHI
            xs = a + b * PHI_STAR
            if -R <= x <= R and window[0] <= xs < window[1]:
                pts.append((a, b, x, xs))
    pts.sort(key=lambda p: p[2])
    return pts


def physical_gaps(points) -> List[float]:
    """Return successive gaps x_{k+1} - x_k along the patch."""
    return [points[i+1][2] - points[i][2] for i in range(len(points) - 1)]


if __name__ == "__main__":
    import sys
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    pts = generate(R)
    print(f"Fibonacci patch R={R}: {len(pts)} points")
    print(f"  first 5: {[(a, b, round(x, 6)) for (a, b, x, _) in pts[:5]]}")
    print(f"  last  5: {[(a, b, round(x, 6)) for (a, b, x, _) in pts[-5:]]}")
    gaps = physical_gaps(pts)
    unique_gaps = sorted(set(round(g, 8) for g in gaps))
    print(f"  unique gaps: {unique_gaps}")
