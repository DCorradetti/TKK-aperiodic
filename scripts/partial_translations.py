"""partial_translations.py
Article_001 — G4: build the partial translation matrices T_t on l^2(Sigma_R).

For a finite patch points = [(a_0, b_0, x_0, _), ..., (a_{n-1}, b_{n-1}, x_{n-1}, _)],
we index basis vectors {delta_x : x in patch} by 0..n-1, and define
    (T_t)_{j, i} = 1 iff x_j = x_i + t and x_j in patch,
i.e. T_t delta_{x_i} = delta_{x_i + t} if x_i + t in patch, else 0.

We work with t expressed in lattice coordinates t = (da, db) so that the
test x_j = x_i + t reduces to (a_j, b_j) = (a_i + da, b_i + db) in Z^2.
Floating-point arithmetic plays no role at this layer.

Propagation filter Delta_M(Sigma_R) = { t in Sigma_R - Sigma_R : |t|_phys <= M }.

Version: 1
"""
from typing import Dict, List, Tuple
from fractions import Fraction
from fields import phi_value
from linalg import Mat, ONE


def delta_M(points, M: float) -> List[Tuple[int, int]]:
    """Difference set filtered by physical norm M. Includes (0, 0)."""
    seen = set()
    out = []
    for i, p in enumerate(points):
        for j, q in enumerate(points):
            da, db = q[0] - p[0], q[1] - p[1]
            x = phi_value(da, db)
            if abs(x) <= M:
                key = (da, db)
                if key not in seen:
                    seen.add(key)
                    out.append(key)
    out.sort(key=lambda t: (phi_value(*t), t))
    return out


def T_matrix(points, t: Tuple[int, int]) -> Mat:
    """Sparse 0/1 partial translation matrix T_t indexed by patch positions."""
    coord_to_idx = {(p[0], p[1]): k for k, p in enumerate(points)}
    out: Mat = {}
    da, db = t
    for i, (a, b, _x, _xs) in enumerate(points):
        target = (a + da, b + db)
        if target in coord_to_idx:
            j = coord_to_idx[target]
            out[(j, i)] = ONE
    return out


def build_generators(points, M: float) -> Dict[Tuple[int, int], Mat]:
    """Return dict t -> T_t for all t in Delta_M(points)."""
    deltas = delta_M(points, M)
    return {t: T_matrix(points, t) for t in deltas}


if __name__ == "__main__":
    import sys
    from cps_fibonacci import generate
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
    M = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    pts = generate(R)
    deltas = delta_M(pts, M)
    print(f"Fibonacci R={R}, M={M}: {len(pts)} points, |Delta_M|={len(deltas)}")
    print(f"  Delta_M (lattice): {deltas}")
    print(f"  Delta_M (physical): {[round(phi_value(*t), 6) for t in deltas]}")
