"""probe.py — quick dimension probe without TKK Instr build."""
import argparse, time, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fields import PHI, INV_PHI
from cps_fibonacci import generate
from partial_translations import build_generators
from star_algebra import build_associative_star_algebra, build_jordan_subalgebra

p = argparse.ArgumentParser()
p.add_argument("--R", type=float, default=20.0)
p.add_argument("--M", type=float, default=2.0)
args = p.parse_args()

window = (-PHI/2.0, PHI/2.0)
t0 = time.time()
pts = generate(args.R, window=window)
print(f"|Sigma_R={args.R}| = {len(pts)}    (gen in {time.time()-t0:.2f}s)")

t1 = time.time()
gens = build_generators(pts, args.M)
print(f"|Delta_M={args.M}| = {len(gens)}    (gen in {time.time()-t1:.2f}s)")
print(f"  generators (lattice coords): {sorted(gens.keys())[:20]}{'...' if len(gens) > 20 else ''}")

t2 = time.time()
basis_A, span_A = build_associative_star_algebra(gens, len(pts), track_transformation=False)
print(f"dim A = {len(basis_A)}    (built in {time.time()-t2:.2f}s)")

t3 = time.time()
basis_J, span_J = build_jordan_subalgebra(basis_A, len(pts))
print(f"dim J = {len(basis_J)}    (built in {time.time()-t3:.2f}s)")

print(f"TOTAL: {time.time()-t0:.2f}s")
