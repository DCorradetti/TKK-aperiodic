"""verify_quasiaddition_jordan_failure.py
Article_001 — exhibit a concrete failure of the Jordan identity for the
quasiaddition algebra on a Fibonacci patch.

The quasiaddition algebra has basis {e_x : x in Sigma_R} and the truncated
commutative product
        e_x ∘ e_y = e_{x+y}  if x+y in Sigma_R, else 0,
extended bilinearly to span_Q{e_x}. This is exactly the "naive" combinatorial
Jordan-algebra candidate suggested by the Patera-Twarock tradition for Lie
algebras of quasicrystals (transposed to the symmetric Jordan setting).

We search for a triple (x, y, z) in Sigma_R such that the Jordan identity
        a^2 ∘ (a ∘ b) − a ∘ (a^2 ∘ b) ≠ 0
fails with a = e_x + e_y and b = e_z. By the polarized formal calculation
both sides equal a sum of formal monomials in the algebra generators, but
the truncation conditions (each monomial requires intermediate sums to lie
in Sigma_R) introduce SIGN-FREE differences in the coefficients that can
make LHS ≠ RHS.

If a failure is found we emit an explicit JSON certificate (negative C2)
to 7. Results/Article_001/run_quasi_failure_<timestamp>/manifest.json.

Version: 1
"""
import argparse, json, os, sys, time, itertools
from datetime import datetime
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "7. Results", "Article_001"))
sys.path.insert(0, HERE)

from fields import PHI, PHI_STAR, INV_PHI
from cps_fibonacci import generate


def in_window(p, window):
    """Window membership in infinite-Sigma model set."""
    star = p[0] + p[1] * PHI_STAR
    return window[0] <= star < window[1]


def mul_quasi(d_a, d_b, window):
    """Truncated multiplication: e_x * e_y = e_{x+y} if x+y in (INFINITE)
    model set, else 0. Truncation is governed by WINDOW membership, not
    patch finiteness, so any failure detected here is genuine (not a patch
    artifact)."""
    out = {}
    for p, ap in d_a.items():
        if ap == 0:
            continue
        for q, bp in d_b.items():
            if bp == 0:
                continue
            s = (p[0] + q[0], p[1] + q[1])
            if in_window(s, window):
                out[s] = out.get(s, Fraction(0)) + ap * bp
    return {k: v for k, v in out.items() if v != 0}


def add_dict(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fraction(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def scale(c, a):
    out = {}
    for k, v in a.items():
        nv = c * v
        if nv != 0:
            out[k] = nv
    return out


def jord(a, b, window):
    ab = mul_quasi(a, b, window)
    ba = mul_quasi(b, a, window)
    return scale(Fraction(1, 2), add_dict(ab, ba))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--R", type=float, default=20.0)
    p.add_argument("--window", choices=("standard", "symmetric"), default="symmetric")
    p.add_argument("--max_search", type=int, default=200000)
    args = p.parse_args()

    if args.window == "standard":
        window = (-1.0, INV_PHI)
        window_label = "[-1, 1/phi)"
    else:
        window = (-PHI / 2.0, PHI / 2.0)
        window_label = "[-phi/2, phi/2)"

    t0 = time.time()
    pts = generate(args.R, window=window)
    coord_to_idx = {(p_[0], p_[1]): k for k, p_ in enumerate(pts)}
    n = len(pts)
    print(f"Fibonacci patch R={args.R}, window={window_label}: {n} points")

    # Search for a (x, y, z) with x, y, z in Sigma such that
    #   a = e_x + e_y, b = e_z
    # gives a Jordan identity failure.
    # We allow x = y (then a = 2 e_x, but coefficient differences may still
    # appear at higher monomial level — keep for completeness).
    failures = []
    n_tested = 0
    for ix in range(n):
        for iy in range(ix, n):
            for iz in range(n):
                n_tested += 1
                if n_tested >= args.max_search:
                    break
                p_x = (pts[ix][0], pts[ix][1])
                p_y = (pts[iy][0], pts[iy][1])
                p_z = (pts[iz][0], pts[iz][1])
                a = {p_x: Fraction(1)}
                a = add_dict(a, {p_y: Fraction(1)})  # handles x == y by summing
                b = {p_z: Fraction(1)}
                a2 = jord(a, a, window)
                lhs = jord(a2, jord(a, b, window), window)
                rhs = jord(a, jord(a2, b, window), window)
                if lhs != rhs:
                    failures.append({
                        "x_latt": p_x, "y_latt": p_y, "z_latt": p_z,
                        "x_phys": round(pts[ix][2], 8),
                        "y_phys": round(pts[iy][2], 8),
                        "z_phys": round(pts[iz][2], 8),
                        "a2": {str(k): str(v) for k, v in a2.items()},
                        "lhs": {str(k): str(v) for k, v in lhs.items()},
                        "rhs": {str(k): str(v) for k, v in rhs.items()},
                    })
                    if len(failures) >= 5:
                        break
            if len(failures) >= 5 or n_tested >= args.max_search:
                break
        if len(failures) >= 5 or n_tested >= args.max_search:
            break

    print(f"  tested {n_tested} (x, y, z) triples in {time.time()-t0:.2f}s")
    print(f"  failures found: {len(failures)}")
    for i, f in enumerate(failures[:3]):
        print(f"  -- failure {i+1}: x={f['x_latt']}({f['x_phys']}) "
              f"y={f['y_latt']}({f['y_phys']}) z={f['z_latt']}({f['z_phys']})")
        print(f"     a^2 = {f['a2']}")
        print(f"     LHS = {f['lhs']}")
        print(f"     RHS = {f['rhs']}")

    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    out_dir = os.path.join(RESULTS_ROOT, f"run_quasi_failure_R{int(args.R)}_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    manifest = {
        "script": "verify_quasiaddition_jordan_failure.py",
        "version": 1,
        "timestamp_utc": ts,
        "command": f"python verify_quasiaddition_jordan_failure.py --R {args.R} --window {args.window}",
        "model_set": "Fibonacci",
        "window": window_label,
        "patch_radius_R": args.R,
        "n_points": n,
        "max_search": args.max_search,
        "n_tested": n_tested,
        "n_failures_found": len(failures),
        "failures": failures,
        "conclusion": ("Jordan identity FAILS for the quasiaddition algebra at the listed (x, y, z)"
                       if failures else "no failure detected within tested range"),
        "runtime_seconds": round(time.time() - t0, 3),
    }
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
    print(f"  certificate written to: {out_dir}")
    return manifest


if __name__ == "__main__":
    main()
