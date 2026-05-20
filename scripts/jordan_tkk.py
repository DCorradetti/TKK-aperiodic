"""jordan_tkk.py
Article_001 — given a basis of the special Jordan algebra J, compute:
  - structure constants of the Jordan product;
  - the inner structure algebra Instr(J) = span{V_{a,b}} as a subspace of End(J);
  - the dimension of the TKK Lie algebra g(J) = J^- + Instr(J) + J^+;
  - identity verifications (Jordan and Jordan triple) by exhaustive or random
    sampling.

Mathematical note. For a SPECIAL Jordan algebra J = A_sa with A an associative
*-algebra, the Jordan identity, the Jordan triple identity, and the Jacobi
identity in the TKK Lie algebra g(J) all hold by classical theorems
(Jacobson 1968 §I.4; Loos 1975 §1; the TKK Jacobi identity reduces to the
Jordan triple identity, which itself follows from associativity of A).

We therefore TEST these identities for IMPLEMENTATION SANITY, not as
mathematical claims. A failure indicates a bug, not a deep new phenomenon.

Conventions:
  a ∘ b   = (a*b + b*a)/2
  {a,b,c} = (a∘b)∘c + a∘(b∘c) − b∘(a∘c)
  V_{a,b}(x) = {a, b, x}

Version: 1
"""
from typing import Dict, List, Tuple
from fractions import Fraction
import itertools, random

from linalg import (
    Mat, mat_to_svec, mat_anticomm_half, mat_jordan_triple, mat_add, mat_sub,
    IndependentSpan, ONE, ZERO,
)

Frac = Fraction

# ----- Jordan structure -------------------------------------------------------

def jordan_structure_constants(basis_J: List[Mat], span_J: IndependentSpan, n: int):
    """Return c[i][j] = tuple of Fractions with
       basis_J[i] ∘ basis_J[j] = sum_l c[i][j][l] basis_J[l].
    Symmetric: c[i][j] = c[j][i]."""
    k = len(basis_J)
    c = [[None] * k for _ in range(k)]
    for i in range(k):
        for j in range(i, k):
            prod = mat_anticomm_half(basis_J[i], basis_J[j])
            coefs = span_J.coords(mat_to_svec(prod, n))
            c[i][j] = coefs
            c[j][i] = coefs
    return c


def verify_jordan_identity(basis_J: List[Mat], n: int, sample=None):
    """Verify a^2 ∘ (a∘b) = a ∘ (a^2 ∘ b) for all (a, b) in basis_J^2, or on
    a random sample of given size. Returns (ok, n_checked, failures)."""
    k = len(basis_J)
    pairs = [(i, j) for i in range(k) for j in range(k)]
    if sample is not None and sample < len(pairs):
        random.seed(0xABCDEF)
        pairs = random.sample(pairs, sample)
    failures = []
    for (i, j) in pairs:
        a = basis_J[i]; b = basis_J[j]
        a2 = mat_anticomm_half(a, a)
        ab = mat_anticomm_half(a, b)
        lhs = mat_anticomm_half(a2, ab)
        a2b = mat_anticomm_half(a2, b)
        rhs = mat_anticomm_half(a, a2b)
        if lhs != rhs:
            failures.append((i, j))
    return (len(failures) == 0, len(pairs), failures)


def verify_jordan_triple_identity(basis_J: List[Mat], n: int, sample: int = 50):
    """Random sample check of
       {a,b,{x,y,z}} = {{a,b,x},y,z} − {x,{b,a,y},z} + {x,y,{a,b,z}}.
    Returns (ok, n_checked, failures). Default sample=50 quintuples."""
    k = len(basis_J)
    random.seed(0x123456)
    quins = [tuple(random.randrange(k) for _ in range(5)) for _ in range(sample)]
    failures = []
    for (ai, bi, xi, yi, zi) in quins:
        a = basis_J[ai]; b = basis_J[bi]; x = basis_J[xi]; y = basis_J[yi]; z = basis_J[zi]
        xyz = mat_jordan_triple(x, y, z)
        lhs = mat_jordan_triple(a, b, xyz)
        abx = mat_jordan_triple(a, b, x)
        bay = mat_jordan_triple(b, a, y)
        abz = mat_jordan_triple(a, b, z)
        rhs = mat_add(mat_sub(mat_jordan_triple(abx, y, z),
                              mat_jordan_triple(x, bay, z)),
                      mat_jordan_triple(x, y, abz))
        if lhs != rhs:
            failures.append((ai, bi, xi, yi, zi))
    return (len(failures) == 0, len(quins), failures)


# ----- Inner structure algebra -----------------------------------------------

def build_instr(basis_J: List[Mat], span_J: IndependentSpan, n: int):
    """Build a basis of Instr(J) = span{V_{a,b} : a, b in J} ⊆ End(J).
    Each candidate V_{i,j} is represented as a k×k rational matrix
    (the matrix of V_{a,b} in basis_J), then vectorized to length k*k.
    Returns (basis_pairs, span_instr) where basis_pairs is a list of (i, j)
    such that V_{basis_J[i], basis_J[j]} is linearly independent of earlier ones.
    """
    k = len(basis_J)
    span_instr = IndependentSpan()
    basis_pairs: List[Tuple[int, int]] = []
    for i in range(k):
        for j in range(k):
            a = basis_J[i]; b = basis_J[j]
            sv = {}
            for m in range(k):
                Vx = mat_jordan_triple(a, b, basis_J[m])
                coefs = span_J.coords(mat_to_svec(Vx, n))
                for l, c in enumerate(coefs):
                    if c != 0:
                        sv[l * k + m] = c
            if span_instr.append_if_new(sv):
                basis_pairs.append((i, j))
    return basis_pairs, span_instr


# ----- TKK assembly ----------------------------------------------------------

def tkk_dimension(dim_J: int, dim_instr: int) -> int:
    return 2 * dim_J + dim_instr


# ----- Symmetry verification --------------------------------------------------

def apply_permutation_to_mat(M: Mat, perm: List[int]) -> Mat:
    """Conjugate sparse matrix M by the permutation: out[(perm[i], perm[j])] = M[(i,j)]."""
    return {(perm[i], perm[j]): v for (i, j), v in M.items()}


def verify_basis_invariance(basis_J: List[Mat], perm: List[int], span_J: IndependentSpan, n: int):
    """Test whether the action M -> P M P^{-1} preserves span(basis_J).
    Returns (ok, action_matrix_on_J_basis_or_None).
    If ok, returns the k×k matrix representing the action on the J-basis."""
    k = len(basis_J)
    action = [[ZERO] * k for _ in range(k)]
    for j, M in enumerate(basis_J):
        Mp = apply_permutation_to_mat(M, perm)
        try:
            coefs = span_J.coords(mat_to_svec(Mp, n))
        except ValueError:
            return False, None
        for i, c in enumerate(coefs):
            action[i][j] = c
    return True, action
