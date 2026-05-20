"""fields.py
Article_001 — TKK Constructions for Aperiodic Jordan Structures over Model Sets
Plan reference: 0. Plans/Article_001/Plan_v1.md  (G1)

Purpose
-------
Real-arithmetic constants for Z[phi] = ring of integers of Q(sqrt(5)).
Floats with double precision are sufficient for window-membership tests
because:
  - patch coordinates (a, b) are bounded integers with |a|, |b| < 200;
  - the window boundary is rational or expressible in Z[phi];
  - the conjugate phi' is irrational, so a + b*phi' never lies exactly
    on a rational/Z[phi] boundary for nonzero (a, b); the separation
    from the boundary is bounded below by ~1/(max|b|) >> 1e-12.

For Jordan/TKK linear algebra all matrix entries are rational (the
partial translation matrices have 0/1 entries; the Jordan product
involves only the factor 1/2). The Z[phi] structure lives in the
geometry of the model set, not in the algebra entries.

Version: 1
"""
import math

SQRT5 = math.sqrt(5.0)
PHI = (1.0 + SQRT5) / 2.0          # golden ratio
PHI_STAR = (1.0 - SQRT5) / 2.0     # algebraic conjugate of phi
INV_PHI = 1.0 / PHI                # = phi - 1

def phi_value(a: int, b: int) -> float:
    """Real value of a + b*phi."""
    return a + b * PHI

def phi_star(a: int, b: int) -> float:
    """Star (algebraic conjugate) value: a + b*phi'."""
    return a + b * PHI_STAR

def norm(a: int, b: int) -> int:
    """Integer norm N(a + b*phi) = (a + b*phi)(a + b*phi') = a^2 + a*b - b^2."""
    return a * a + a * b - b * b
