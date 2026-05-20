"""run_fibonacci.py
Article_001 — main orchestration for Fibonacci experiments.

Pipeline (Plan_v1 §Strategy):
  G2  generate Fibonacci patch Sigma_R
  G3  identify interior (collared) Sigma_R^int
  G4  build partial-translation operators T_t for t in Delta_M
  G5  close under multiplication -> *-algebra A_R^(M); take symmetric part J
  G6  verify Jordan identity (exhaustive or sampled)
  G7  build Instr(J) and compute dim TKK
  Symmetry: verify Z_2 inversion x -> -x acts by automorphism

Outputs:
  - JSON certificate written to
        7. Results/Article_001/run_fibonacci_<timestamp>/manifest.json
  - human-readable summary to stdout

Usage:
  python run_fibonacci.py --R 20 --M 1
  python run_fibonacci.py --R 40 --M 1
  python run_fibonacci.py --R 20 --M 2
"""
import argparse, json, os, sys, time, hashlib
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "7. Results", "Article_001"))

sys.path.insert(0, HERE)

from fields import PHI, INV_PHI, phi_value
from cps_fibonacci import generate, physical_gaps
from partial_translations import delta_M, build_generators
from star_algebra import build_associative_star_algebra, build_jordan_subalgebra
from jordan_tkk import (
    jordan_structure_constants, verify_jordan_identity,
    verify_jordan_triple_identity, build_instr, tkk_dimension,
    verify_basis_invariance,
)


def hash_obj(obj) -> str:
    h = hashlib.sha256()
    h.update(json.dumps(obj, sort_keys=True, default=str).encode("utf-8"))
    return h.hexdigest()[:16]


def fibonacci_inversion_permutation(points):
    """Inversion x -> -x is the natural Z_2 acting on Fibonacci (with
    symmetric window). For our standard window [-1, 1/phi) it does NOT preserve
    the patch (window asymmetric). Try the symmetric window [-1/phi, 1/phi)
    if needed; for the default window we still report whether the permutation
    is total or partial.
    Returns perm: list of length len(points) where perm[i] is the patch index
    of -points[i], or -1 if not in patch.
    """
    coord_to_idx = {(p[0], p[1]): k for k, p in enumerate(points)}
    perm = []
    for (a, b, _x, _xs) in points:
        target = (-a, -b)
        perm.append(coord_to_idx.get(target, -1))
    return perm


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--R", type=float, default=10.0)
    p.add_argument("--M", type=float, default=1.0)
    p.add_argument("--window", choices=("standard", "symmetric"), default="symmetric",
                   help="standard = [-1, 1/phi)  (Z_2 broken); symmetric = [-phi/2, phi/2)")
    p.add_argument("--jordan_sample", type=int, default=None,
                   help="if set, sample this many basis pairs for Jordan identity")
    p.add_argument("--jts_sample", type=int, default=30)
    p.add_argument("--no_instr", action="store_true",
                   help="skip building Instr(J) (use for very large J)")
    args = p.parse_args()

    if args.window == "standard":
        window = (-1.0, INV_PHI)
        window_label = "[-1, 1/phi)"
    else:
        window = (-PHI / 2.0, PHI / 2.0)
        window_label = "[-phi/2, phi/2)"

    t0 = time.time()
    print(f"=== Fibonacci experiment: R={args.R}, M={args.M}, window={window_label} ===")
    pts = generate(args.R, window=window)
    print(f"  patch points (|Sigma_R|)        : {len(pts)}")
    print(f"  unique physical gaps           : {sorted(set(round(g, 6) for g in physical_gaps(pts)))}")

    deltas = delta_M(pts, args.M)
    print(f"  generators count (|Delta_M|)    : {len(deltas)}")
    gens = build_generators(pts, args.M)

    t1 = time.time()
    basis_A, span_A = build_associative_star_algebra(gens, len(pts))
    t2 = time.time()
    dim_A = len(basis_A)
    print(f"  dim A_R^(M)                     : {dim_A}     (built in {t2-t1:.2f}s)")

    basis_J, span_J = build_jordan_subalgebra(basis_A, len(pts))
    t3 = time.time()
    dim_J = len(basis_J)
    print(f"  dim J_R^(M) = (A_R^(M))_sa      : {dim_J}     (built in {t3-t2:.2f}s)")

    # Identity verifications
    ok_jord, n_jord, fails_jord = verify_jordan_identity(basis_J, len(pts), args.jordan_sample)
    t4 = time.time()
    print(f"  Jordan identity                 : {'PASS' if ok_jord else 'FAIL'} (checked {n_jord} pairs in {t4-t3:.2f}s)")
    if not ok_jord:
        print(f"  Jordan failures (first 5)       : {fails_jord[:5]}")

    ok_jts, n_jts, fails_jts = verify_jordan_triple_identity(basis_J, len(pts), args.jts_sample)
    t5 = time.time()
    print(f"  Jordan triple identity (sample) : {'PASS' if ok_jts else 'FAIL'} (sample {n_jts} in {t5-t4:.2f}s)")
    if not ok_jts:
        print(f"  JTS failures (first 5)          : {fails_jts[:5]}")

    instr_basis_pairs = []
    dim_instr = 0
    if not args.no_instr:
        instr_basis_pairs, span_instr = build_instr(basis_J, span_J, len(pts))
        dim_instr = len(instr_basis_pairs)
        t6 = time.time()
        print(f"  dim Instr(J)                    : {dim_instr}    (built in {t6-t5:.2f}s)")
    dim_TKK = tkk_dimension(dim_J, dim_instr)
    print(f"  dim TKK(J)                      : {dim_TKK}    (= 2*dim J + dim Instr)")

    # Z_2 inversion symmetry
    perm = fibonacci_inversion_permutation(pts)
    perm_total = all(p >= 0 for p in perm)
    sym_certificate = {"perm_total": perm_total}
    if perm_total:
        ok_sym, action = verify_basis_invariance(basis_J, perm, span_J, len(pts))
        sym_certificate["jordan_invariance"] = bool(ok_sym)
    else:
        # window is asymmetric; report the deficit
        n_missing = sum(1 for p in perm if p < 0)
        sym_certificate["jordan_invariance"] = None
        sym_certificate["n_missing"] = n_missing
    print(f"  Z_2 inversion symmetry          : {sym_certificate}")

    # Write certificate
    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    out_dir = os.path.join(RESULTS_ROOT, f"run_fibonacci_R{int(args.R)}_M{int(args.M)}_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    manifest = {
        "script": "run_fibonacci.py",
        "version": 1,
        "timestamp_utc": ts,
        "command": f"python run_fibonacci.py --R {args.R} --M {args.M}",
        "python_version": sys.version.split()[0],
        "model_set": "Fibonacci",
        "window": window_label,
        "field": "Z[phi]/Q",
        "patch_radius_R": args.R,
        "propagation_max_M": args.M,
        "n_points": len(pts),
        "delta_M_size": len(deltas),
        "dim_A": dim_A,
        "dim_J": dim_J,
        "dim_instr": dim_instr,
        "dim_TKK": dim_TKK,
        "jordan_identity": {
            "verified": bool(ok_jord),
            "n_checked": n_jord,
            "n_failures": len(fails_jord),
        },
        "jordan_triple_identity_sample": {
            "verified": bool(ok_jts),
            "n_checked": n_jts,
            "n_failures": len(fails_jts),
        },
        "z2_inversion_symmetry": sym_certificate,
        "runtime_seconds": round(time.time() - t0, 3),
        "patch_first_5_points": [(p[0], p[1], round(p[2], 6)) for p in pts[:5]],
        "patch_last_5_points":  [(p[0], p[1], round(p[2], 6)) for p in pts[-5:]],
    }
    manifest["self_hash"] = hash_obj({k: v for k, v in manifest.items() if k != "self_hash"})
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
    with open(os.path.join(out_dir, "log.txt"), "w", encoding="utf-8") as f:
        f.write(f"Fibonacci R={args.R}, M={args.M} — see manifest.json\n")
    print(f"  certificate written to: {out_dir}")
    return manifest


if __name__ == "__main__":
    main()
