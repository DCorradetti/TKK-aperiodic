"""linalg.py
Exact rational linear algebra utilities for Article_001.

Vectors are represented as DICTIONARIES {index -> Fraction}. Matrices used
in this project (partial-translation matrices and products thereof) are 0/1
with at most one nonzero per row, so their vectorized form is extremely
sparse and a dict representation is efficient.

This module provides:
  - sparse vector arithmetic (sv_add, sv_sub, sv_scalar, sv_eq, sv_norm0)
  - an IndependentSpan class that maintains a basis of a linear subspace
    and supports two operations:
      * append_if_new(v): return True if v was added (i.e. v not in current
        span), False if v is already in span; in the latter case stores
        the linear-combination coefficients in self.last_coords.
      * coords(v): returns a tuple of Fractions of length len(self.basis)
        with v = sum_j coords[j] * self.basis[j] (raises if v not in span).

Version: 1
"""
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

Frac = Fraction
SVec = Dict[int, Frac]   # sparse vector

ZERO = Frac(0)
ONE = Frac(1)
HALF = Frac(1, 2)

def sv_copy(u: SVec) -> SVec:
    return dict(u)

def sv_eq(u: SVec, v: SVec) -> bool:
    return u == v

def sv_norm0(u: SVec) -> int:
    return sum(1 for x in u.values() if x != 0)

def sv_add_inplace(u: SVec, v: SVec, scale: Frac = ONE) -> None:
    """u += scale * v."""
    for k, x in v.items():
        new = u.get(k, ZERO) + scale * x
        if new == 0:
            if k in u:
                del u[k]
        else:
            u[k] = new

def sv_scalar(c: Frac, u: SVec) -> SVec:
    if c == 0:
        return {}
    return {k: c * x for k, x in u.items() if x != 0}

def sv_first_nonzero(u: SVec) -> Optional[int]:
    if not u:
        return None
    keys = [k for k, x in u.items() if x != 0]
    if not keys:
        return None
    return min(keys)


class IndependentSpan:
    """Incrementally maintain a linearly independent set of sparse vectors.

    Tracks a parallel transformation so that we can express any new vector
    that lies in the span as an explicit combination of the stored vectors.

    Internal invariants:
      basis[j]       — the j-th original stored vector (sparse).
      rref[i]        — sparse vector, the i-th row of the reduced echelon form,
                       normalized so that rref[i][pivots[i]] == 1.
      pivots[i]      — pivot column of rref[i]; the pivots list is strictly
                       increasing with i.
      T[i]           — sparse vector of length len(basis); rref[i] equals
                       sum_j T[i][j] * basis[j].
    """

    def __init__(self, track_transformation: bool = True):
        self.basis: List[SVec] = []
        self.rref: List[SVec] = []
        self.pivots: List[int] = []
        self.T: List[SVec] = []
        self.track_T = track_transformation
        self.last_coords: Optional[Tuple[Frac, ...]] = None

    def __len__(self) -> int:
        return len(self.basis)

    def _reduce_with_combination(self, v: SVec):
        """Reduce v against rref, returning (u, comb) where:
          - u is the residual (sparse vector); u==0 iff v is in current span.
          - comb is a sparse vector of indices into basis; if u==0 then
            v = sum_j comb[j] * basis[j].
        comb is computed only when track_T is True; otherwise it's an empty dict.
        """
        u = sv_copy(v)
        comb: SVec = {}
        if self.track_T:
            for i, (r, p) in enumerate(zip(self.rref, self.pivots)):
                coef = u.get(p, ZERO)
                if coef != 0:
                    sv_add_inplace(u, r, -coef)
                    sv_add_inplace(comb, self.T[i], coef)
        else:
            for r, p in zip(self.rref, self.pivots):
                coef = u.get(p, ZERO)
                if coef != 0:
                    sv_add_inplace(u, r, -coef)
        return u, comb

    def append_if_new(self, v: SVec) -> bool:
        """Test v: if v is linearly independent of current basis, append it and
        return True. Otherwise return False and set self.last_coords to the
        tuple of coefficients expressing v in the existing basis.
        """
        u, comb = self._reduce_with_combination(v)
        if not u:
            # v already in span
            self.last_coords = tuple(comb.get(j, ZERO) for j in range(len(self.basis)))
            return False
        # find first nonzero
        p = sv_first_nonzero(u)
        if p is None:
            self.last_coords = tuple(comb.get(j, ZERO) for j in range(len(self.basis)))
            return False
        # normalize u
        upiv = u[p]
        if upiv != ONE:
            inv = ONE / upiv
            for k in list(u.keys()):
                u[k] *= inv
            for k in list(comb.keys()):
                comb[k] *= inv
        # new basis index
        new_idx = len(self.basis)
        # u = v_new + sum_j comb[j] * basis[j]  (in original basis coords)
        # but we want: u = sum_{j <= new_idx} T_new[j] * basis_after[j]
        # where basis_after = basis + [v].
        # So T_new[new_idx] = 1, T_new[j] = comb[j] for j < new_idx ... wait sign:
        # We reduced u from v_new and ADDED coef*T[i] to comb each step.
        # So u_final = v - sum_i coef_i * rref[i] = v - sum_j comb[j] * basis[j] (??)
        # Let me re-verify: at step i, we did u -= coef * rref[i] and comb += coef * T[i].
        # Total: u_final = v - sum_i coef_i * rref[i] = v - sum_i coef_i * sum_j T[i][j] basis[j]
        #               = v - sum_j (sum_i coef_i T[i][j]) basis[j]
        #               = v - sum_j comb[j] basis[j]    (since comb[j] = sum_i coef_i T[i][j])
        # So u_final = v - sum_j comb[j] basis[j]   →   u_final + sum_j comb[j] basis[j] = v
        # i.e. v = u_final + sum_j comb[j] basis[j].
        # In the new basis with v appended as basis_after[new_idx]:
        #   basis_after[new_idx] = v = u_final + sum_j comb[j] basis[j].
        #   So u_final = basis_after[new_idx] - sum_j comb[j] basis[j] = ... we want u_final expressed in basis_after.
        # Equivalently: new_rref_row = u_norm = (1/upiv) * u_final
        #   u_norm = (1/upiv)(basis_after[new_idx] - sum_j comb_pre_norm[j] basis_after[j])
        #          = basis_after[new_idx]/upiv - (1/upiv) sum_j comb_pre_norm[j] basis_after[j]
        # We already multiplied comb by 1/upiv above (loop `comb[k] *= inv`).
        # So u_norm = basis_after[new_idx]/upiv - sum_j comb[j] basis_after[j]
        # T_new[new_idx] = 1/upiv, T_new[j] = -comb[j] for j < new_idx.
        T_new: SVec = {}
        for j, c in comb.items():
            if c != 0:
                T_new[j] = -c
        T_new[new_idx] = ONE / upiv
        # back-substitute the new pivot into existing rref rows
        for i, (r, q) in enumerate(zip(self.rref, self.pivots)):
            coef = r.get(p, ZERO)
            if coef != 0:
                sv_add_inplace(r, u, -coef)
                if self.track_T:
                    sv_add_inplace(self.T[i], T_new, -coef)
        # insert in pivot-ascending order
        idx = 0
        while idx < len(self.pivots) and self.pivots[idx] < p:
            idx += 1
        self.basis.append(sv_copy(v))
        self.rref.insert(idx, u)
        self.pivots.insert(idx, p)
        if self.track_T:
            self.T.insert(idx, T_new)
        self.last_coords = None
        return True

    def coords(self, v: SVec) -> Tuple[Frac, ...]:
        """Express v in self.basis; raises ValueError if v not in span."""
        u, comb = self._reduce_with_combination(v)
        if u:
            # check if u is "all zero"
            if any(x != 0 for x in u.values()):
                raise ValueError("vector not in span")
        return tuple(comb.get(j, ZERO) for j in range(len(self.basis)))


# --- Matrix utilities on sparse 0/1 matrices for partial translations -------
# We store matrices as dict {(i,j): Fraction}.

Mat = Dict[Tuple[int, int], Frac]

def mat_to_svec(A: Mat, n: int) -> SVec:
    return {i * n + j: v for (i, j), v in A.items() if v != 0}

def svec_to_mat(v: SVec, n: int) -> Mat:
    out: Mat = {}
    for k, x in v.items():
        if x != 0:
            out[(k // n, k % n)] = x
    return out

def mat_zero() -> Mat:
    return {}

def mat_eq(A: Mat, B: Mat) -> bool:
    return A == B

def mat_add(A: Mat, B: Mat) -> Mat:
    out = dict(A)
    for k, v in B.items():
        new = out.get(k, ZERO) + v
        if new == 0:
            if k in out:
                del out[k]
        else:
            out[k] = new
    return out

def mat_sub(A: Mat, B: Mat) -> Mat:
    out = dict(A)
    for k, v in B.items():
        new = out.get(k, ZERO) - v
        if new == 0:
            if k in out:
                del out[k]
        else:
            out[k] = new
    return out

def mat_scalar(c: Frac, A: Mat) -> Mat:
    if c == 0:
        return {}
    return {k: c * v for k, v in A.items()}

def mat_mul(A: Mat, B: Mat) -> Mat:
    """Sparse matrix product. Indexes columns of A and rows of B via index maps."""
    # Build map from row->{(col,val)} for B
    B_rows: Dict[int, Dict[int, Frac]] = {}
    for (i, j), v in B.items():
        B_rows.setdefault(i, {})[j] = v
    out: Mat = {}
    for (i, k), aval in A.items():
        if k not in B_rows:
            continue
        for j, bval in B_rows[k].items():
            key = (i, j)
            new = out.get(key, ZERO) + aval * bval
            if new == 0:
                if key in out:
                    del out[key]
            else:
                out[key] = new
    return out

def mat_transpose(A: Mat) -> Mat:
    return {(j, i): v for (i, j), v in A.items()}

def mat_anticomm_half(A: Mat, B: Mat) -> Mat:
    """(A B + B A) / 2."""
    AB = mat_mul(A, B)
    BA = mat_mul(B, A)
    s = mat_add(AB, BA)
    return mat_scalar(HALF, s)

def mat_commutator(A: Mat, B: Mat) -> Mat:
    AB = mat_mul(A, B)
    BA = mat_mul(B, A)
    return mat_sub(AB, BA)

def mat_jordan_triple(A: Mat, B: Mat, C: Mat) -> Mat:
    """Jordan triple from Jordan product: (a∘b)∘c + a∘(b∘c) - b∘(a∘c)."""
    ab = mat_anticomm_half(A, B)
    bc = mat_anticomm_half(B, C)
    ac = mat_anticomm_half(A, C)
    t1 = mat_anticomm_half(ab, C)
    t2 = mat_anticomm_half(A, bc)
    t3 = mat_anticomm_half(B, ac)
    return mat_sub(mat_add(t1, t2), t3)
