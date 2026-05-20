# Scripts — Article_001

Reproducibility package for *Tits–Kantor–Koecher Constructions for Aperiodic Jordan Structures over Model Sets* (`1. Drafts/Article_001/Article_001_v?.tex`).

## Requirements

- Python 3.11
- No third-party packages strictly required at the algebra layer (uses Python `fractions.Fraction` and standard library).
- Optional: SymPy 1.13, NumPy 1.26 (only used for top-level reporting; the algebra layer is pure-Python exact arithmetic over `Fraction`).

## Layout

```
2. Scripts/Article_001/
  fields.py                              — Z[phi] real-arithmetic constants
  cps_fibonacci.py                       — Fibonacci CPS generator
  partial_translations.py                — T_t partial-translation matrices
  linalg.py                              — exact-rational sparse linear algebra
  star_algebra.py                        — *-algebra and Jordan-subalgebra builders
  jordan_tkk.py                          — Instr(J), Jordan/JTS identities, symmetry
  probe.py                               — quick dimension probe
  run_fibonacci.py                       — main Fibonacci experiment driver
  verify_quasiaddition_jordan_failure.py — search for quasiaddition Jordan failures
```

## Reproducing the dimension table (Theorem 6.1 / Theorem 6.2 in the draft)

For the aperiodic core ($M < \varphi$):

```
python run_fibonacci.py --R 10  --M 1 --window symmetric
python run_fibonacci.py --R 20  --M 1 --window symmetric
python run_fibonacci.py --R 40  --M 1 --window symmetric
python run_fibonacci.py --R 80  --M 1 --window symmetric
python run_fibonacci.py --R 160 --M 1 --window symmetric
```

For the saturation regime ($M \ge \varphi$):

```
python run_fibonacci.py --R 10 --M 2 --window symmetric --no_instr
python run_fibonacci.py --R 20 --M 2 --window symmetric --no_instr
```

(The `--no_instr` flag skips the explicit construction of $\Instr(J)$, since at saturation the dimension is known a priori: $n^2$ with $n = |\Sigma_R|$, see Remark 3.1 of the draft.)

Each invocation writes a JSON manifest to a fresh timestamped folder in `7. Results/Article_001/`.

## Reproducing the quasiaddition Jordan-failure counterexample (Proposition 6.1)

```
python verify_quasiaddition_jordan_failure.py --R 8 --window symmetric --max_search 100000
```

This emits five explicit counterexamples within the first 32 tested $(x, y, z)$ triples. The minimal one is

```
x = (-2, -3)  i.e. x = -2 - 3*phi
y = (-1, -2)  i.e. y = -1 - 2*phi
z = ( 0,  1)  i.e. z = phi
```

with residue $e_{(-5, -7)} = e_{-5 - 7\varphi}$.

## Quick dimension probe (no certificate)

For exploratory work:

```
python probe.py --R 30 --M 1
```

prints `dim A_R^(M)`, `dim J_R^(M)`, and runtime without building $\Instr(J)$ or writing a manifest.

## Determinism and exact arithmetic

- All algebra-layer arithmetic is performed in `Fraction`. No floating point is used except for window-membership tests (where comparisons against irrational boundaries are stable to far more digits than the precision of `float`, since point coordinates are integer linear combinations of $\varphi$).
- Random sampling for `verify_jordan_identity` (when `--jordan_sample` is set) and for `verify_jordan_triple_identity` is seeded deterministically (`random.seed(0xABCDEF)` and `random.seed(0x123456)`).
- The seed of the quasiaddition search is `0x42` (set in the wrapping script).
- No multithreading; runtimes are reproducible up to OS scheduling jitter.

## Manifest schema

Each `manifest.json` contains:

- `script`, `version`, `timestamp_utc`, `command`
- `python_version`
- model-set parameters: `model_set`, `window`, `patch_radius_R`, `propagation_max_M`
- algebra dimensions: `dim_A`, `dim_J`, `dim_instr`, `dim_TKK`
- identity verifications: `jordan_identity`, `jordan_triple_identity_sample` (each with `verified`, `n_checked`, `n_failures`)
- symmetry: `z2_inversion_symmetry` (with `perm_total`, `jordan_invariance`)
- `runtime_seconds`, patch samples, `self_hash`

The `self_hash` is a short SHA-256 hash of the manifest dict (excluding itself); it is intended to detect accidental tampering with archived runs.
