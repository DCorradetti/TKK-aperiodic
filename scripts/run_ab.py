"""run_ab.py
Article_001 — Ammann–Beenker 2D experiment driver (iteration v3).

Pipeline mirrors run_fibonacci.py but adapted to the 4D lattice indexing
(m, n, p, q) ∈ Z^4 with the Z[sqrt 2] cut-and-project scheme of
cps_ammann_beenker.py.

Outputs a JSON manifest to 7. Results/Article_001/run_ab_R<R>_M<M>_<ts>/.
"""
import argparse, json, math, os, sys, time, hashlib, itertools, random
from datetime import datetime
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "7. Results", "Article_001"))
sys.path.insert(0, HERE)

from cps_ammann_beenker import generate, physical_norm, SQRT2
from linalg import (
    Mat, mat_to_svec, mat_anticomm_half, mat_jordan_triple, mat_add, mat_sub,
    mat_mul, IndependentSpan, ONE, ZERO,
)


def delta_M_4d(points, M):
    n_lat = 4
    seen = set()
    out = []
    for p in points:
        plat = (p[0], p[1], p[2], p[3])
        for q in points:
            qlat = (q[0], q[1], q[2], q[3])
            d = (qlat[0]-plat[0], qlat[1]-plat[1], qlat[2]-plat[2], qlat[3]-plat[3])
            if d in seen:
                continue
            if physical_norm(d) <= M + 1e-12:
                seen.add(d)
                out.append(d)
    out.sort(key=lambda t: (physical_norm(t), t))
    return out


def T_matrix_4d(points, t):
    coord_to_idx = {(p[0], p[1], p[2], p[3]): k for k, p in enumerate(points)}
    out: Mat = {}
    for i, p in enumerate(points):
        target = (p[0]+t[0], p[1]+t[1], p[2]+t[2], p[3]+t[3])
        if target in coord_to_idx:
            j = coord_to_idx[target]
            out[(j, i)] = ONE
    return out


def build_associative_star_algebra_inline(gens_dict, n):
    span = IndependentSpan(track_transformation=False)
    basis_mats = []
    for t, T in gens_dict.items():
        v = mat_to_svec(T, n)
        if span.append_if_new(v):
            basis_mats.append(dict(T))
    gens_list = list(gens_dict.values())
    frontier = list(basis_mats)
    while frontier:
        new_added = []
        for A in frontier:
            for g in gens_list:
                for P in (mat_mul(A, g), mat_mul(g, A)):
                    if not P:
                        continue
                    v = mat_to_svec(P, n)
                    if span.append_if_new(v):
                        basis_mats.append(P)
                        new_added.append(P)
        frontier = new_added
    return basis_mats, span


def build_jordan_inline(basis_A, n):
    from linalg import mat_transpose, mat_scalar, HALF
    span_J = IndependentSpan()
    basis_J = []
    for M in basis_A:
        Mt = mat_transpose(M)
        S = mat_scalar(HALF, mat_add(M, Mt))
        if not S:
            continue
        v = mat_to_svec(S, n)
        if span_J.append_if_new(v):
            basis_J.append(S)
    return basis_J, span_J


def verify_jordan_identity_sample(basis_J, n, sample=64):
    k = len(basis_J)
    pairs = [(i, j) for i in range(k) for j in range(k)]
    if sample is not None and sample < len(pairs):
        random.seed(0xABCDEF)
        pairs = random.sample(pairs, sample)
    failures = 0
    for (i, j) in pairs:
        a = basis_J[i]; b = basis_J[j]
        a2 = mat_anticomm_half(a, a)
        ab = mat_anticomm_half(a, b)
        lhs = mat_anticomm_half(a2, ab)
        a2b = mat_anticomm_half(a2, b)
        rhs = mat_anticomm_half(a, a2b)
        if lhs != rhs:
            failures += 1
    return (failures == 0, len(pairs), failures)


def apply_d8_rotation(points, k):
    """Apply rotation by 2*pi*k/8 to the patch: lattice automorphism induced
    by physical-space rotation by k*pi/4. In the standard A-B CPS the rotation
    by pi/4 cycles the basis vectors e_0 -> e_1 -> e_2 -> e_3 -> -e_0, where
    e_j = (cos(j*pi/4), sin(j*pi/4))^T. The lattice action is therefore
    (m, n, p, q) -> (-q, m, n, p) for k=1 (one quarter turn of the A-B basis).
    For k=4 (i.e., rotation by pi) it is (m, n, p, q) -> (-m, -n, -p, -q),
    which corresponds to inversion through origin (CPT-like).

    For simplicity we test the simplest D_8 element: inversion (k=4),
    which acts on the lattice by (m, n, p, q) -> (-m, -n, -p, -q) and on the
    patch by x -> -x. The octagonal window is invariant under this map.
    """
    if k != 4:
        raise NotImplementedError("only inversion (k=4) implemented in v3")
    coord_to_idx = {(p[0], p[1], p[2], p[3]): i for i, p in enumerate(points)}
    perm = []
    for p in points:
        target = (-p[0], -p[1], -p[2], -p[3])
        perm.append(coord_to_idx.get(target, -1))
    return perm


def verify_basis_invariance(basis_J, perm, span_J, n):
    from linalg import ZERO
    k = len(basis_J)
    for j, M in enumerate(basis_J):
        Mp = {(perm[i], perm[ii]): v for (i, ii), v in M.items()}
        try:
            span_J.coords(mat_to_svec(Mp, n))
        except ValueError:
            return False
    return True


def hash_obj(obj) -> str:
    h = hashlib.sha256()
    h.update(json.dumps(obj, sort_keys=True, default=str).encode("utf-8"))
    return h.hexdigest()[:16]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--R", type=float, default=4.0)
    p.add_argument("--M", type=float, default=1.0)
    p.add_argument("--apothem", type=float, default=1.0 + 1.0 / SQRT2)
    p.add_argument("--jordan_sample", type=int, default=64)
    p.add_argument("--no_instr", action="store_true")
    args = p.parse_args()

    t0 = time.time()
    print(f"=== Ammann-Beenker experiment: R={args.R}, M={args.M}, apothem={args.apothem:.4f} ===")
    pts = generate(args.R, apothem=args.apothem)
    n = len(pts)
    print(f"  patch points (|Sigma_R|)        : {n}")

    deltas = delta_M_4d(pts, args.M)
    print(f"  generators count (|Delta_M|)    : {len(deltas)}")
    if len(deltas) <= 20:
        for d in deltas:
            print(f"     lattice {d}  phys_norm={physical_norm(d):.4f}")
    gens = {t: T_matrix_4d(pts, t) for t in deltas}

    t1 = time.time()
    basis_A, span_A = build_associative_star_algebra_inline(gens, n)
    t2 = time.time()
    dim_A = len(basis_A)
    print(f"  dim A_R^(M)                     : {dim_A}     (built in {t2-t1:.2f}s)")

    basis_J, span_J = build_jordan_inline(basis_A, n)
    t3 = time.time()
    dim_J = len(basis_J)
    print(f"  dim J_R^(M) = (A_R^(M))_sa      : {dim_J}     (built in {t3-t2:.2f}s)")

    ok_jord, n_jord, fails_jord = verify_jordan_identity_sample(basis_J, n, args.jordan_sample)
    t4 = time.time()
    print(f"  Jordan identity (sample)        : {'PASS' if ok_jord else 'FAIL'} (checked {n_jord} pairs in {t4-t3:.2f}s)")

    dim_instr = None
    if not args.no_instr and dim_J <= 25:
        from jordan_tkk import build_instr
        basis_pairs, span_instr = build_instr(basis_J, span_J, n)
        dim_instr = len(basis_pairs)
        print(f"  dim Instr(J)                    : {dim_instr}")
    else:
        print(f"  dim Instr(J)                    : skipped (J too large or --no_instr)")
    dim_TKK = (2 * dim_J + dim_instr) if dim_instr is not None else None

    # Inversion symmetry (k=4 only; perfect octagon invariance is sufficient).
    perm = apply_d8_rotation(pts, 4)
    perm_total = all(p >= 0 for p in perm)
    sym = {"k": 4, "perm_total": perm_total}
    if perm_total:
        ok_sym = verify_basis_invariance(basis_J, perm, span_J, n)
        sym["jordan_invariance"] = bool(ok_sym)
    print(f"  inversion (k=4) symmetry        : {sym}")

    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    out_dir = os.path.join(RESULTS_ROOT, f"run_ab_R{int(args.R)}_M{int(args.M)}_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    manifest = {
        "script": "run_ab.py",
        "version": 1,
        "timestamp_utc": ts,
        "command": f"python run_ab.py --R {args.R} --M {args.M}",
        "python_version": sys.version.split()[0],
        "model_set": "Ammann-Beenker",
        "lattice": "Z^4 via Z[sqrt(2)] + i Z[sqrt(2)]",
        "window": f"regular octagon, apothem {args.apothem:.6f}",
        "field": "Z[sqrt(2)]/Q",
        "patch_radius_R": args.R,
        "propagation_max_M": args.M,
        "n_points": n,
        "delta_M_size": len(deltas),
        "delta_M_samples": deltas[:20],
        "dim_A": dim_A,
        "dim_J": dim_J,
        "dim_instr": dim_instr,
        "dim_TKK": dim_TKK,
        "jordan_identity_sample": {
            "verified": bool(ok_jord),
            "n_checked": n_jord,
            "n_failures": fails_jord,
        },
        "inversion_symmetry": sym,
        "runtime_seconds": round(time.time() - t0, 3),
    }
    manifest["self_hash"] = hash_obj({k: v for k, v in manifest.items() if k != "self_hash"})
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
    print(f"  certificate written to: {out_dir}")
    return manifest


if __name__ == "__main__":
    main()
