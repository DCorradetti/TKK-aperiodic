# TKK Aperiodic
[![License: MIT](https://img.shields.io/github/license/DCorradetti/TKK-aperiodic?style=flat-square&color=blue)](LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/DCorradetti/TKK-aperiodic?style=flat-square)](https://github.com/DCorradetti/TKK-aperiodic/commits/main)

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![SymPy 1.13](https://img.shields.io/badge/SymPy-1.13-3B5526?style=flat-square)](https://www.sympy.org/)
[![NumPy 1.26](https://img.shields.io/badge/NumPy-1.26-013243?style=flat-square&logo=numpy&logoColor=white)](https://numpy.org/)
[![Exact arithmetic](https://img.shields.io/badge/arithmetic-exact%20%E2%84%9A-success?style=flat-square)](#determinism)
[![Reproducible](https://img.shields.io/badge/runs-reproducible-success?style=flat-square)](#determinism)

[![Companion paper (arXiv)](https://img.shields.io/badge/arXiv-2303.12219-b31b1b?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2303.12219)
[![Companion paper (DOI)](https://img.shields.io/badge/DOI-10.1016%2Fj.geomphys.2025.105645-blue?style=flat-square)](https://doi.org/10.1016/j.geomphys.2025.105645)
[![Topic](https://img.shields.io/badge/topic-Jordan%20algebras%20%E2%80%A2%20TKK%20%E2%80%A2%20quasicrystals-7B3FE4?style=flat-square)](#what-is-in-here)

Reproducibility package for

> D. Corradetti, *Tits–Kantor–Koecher Constructions for Aperiodic Jordan Structures over Model Sets*, 2026 (manuscript in preparation).

This repository contains the Python scripts, computational certificates and run manifests that support the manuscript. It is a companion to the prior publication

> D. Corradetti, D. Chester, R. Aschheim, K. Irwin, *Jordan algebras over icosahedral cut-and-project quasicrystals*, J. Geom. Phys. (2025) 105645; arXiv:[2303.12219](https://arxiv.org/abs/2303.12219). DOI: [10.1016/j.geomphys.2025.105645](https://doi.org/10.1016/j.geomphys.2025.105645).

with which it shares the cut-and-project / Jordan-algebra context but uses a fundamentally different (operator-algebraic, propagation-filtered) construction.

## What is in here

The manuscript constructs, for every regular cut-and-project model set $\Sigma$ and every propagation bound $M > 0$:

- a propagation-filtered $\ast$-algebra $\mathcal A_\Sigma^{(M)} \subset B(\ell^2(\Sigma))$ of partial translations,
- a special Jordan algebra $J_\Sigma^{(M)} = (\mathcal A_\Sigma^{(M)})_{\mathrm{sa}}$ with the symmetrized product $\tfrac12(AB + BA)$,
- the Tits–Kantor–Koecher Lie algebra $\mathfrak g_\Sigma^{(M)} = \TKK(J_\Sigma^{(M)})$, $3$-graded.

For the Fibonacci chain with symmetric window the paper proves:

- **Stability (Theorem 6.1).** For every $M < \varphi$ and patch radius $R \ge 6$:
  $$\dim_{\mathbb Q} \mathcal A_{\Sigma_R}^{(M)} = 5,\quad \dim_{\mathbb Q} J_{\Sigma_R}^{(M)} = 4,\quad \dim_{\mathbb Q} \mathfrak g_{\Sigma_R}^{(M)} = 13.$$
- **Local-core identification (Theorem 6.3).** $J_\Sigma^{(M < \varphi)} \cong \mathbb Q \oplus H_2(\mathbb Q)$, hence
  $$\mathfrak g_\Sigma^{(M < \varphi)} \cong \mathfrak{sl}(2, \mathbb Q) \oplus \mathfrak{sp}(4, \mathbb Q).$$
- **Saturation (Theorem 6.4).** For $M \ge \varphi$, $J_{\Sigma_R}^{(M)} = \Sym_n(\mathbb Q)$ and $\mathfrak g_{\Sigma_R}^{(M)} \cong \mathfrak{sp}(2n, \mathbb Q)$ where $n = |\Sigma_R|$.
- **Quasiaddition Jordan failure (Proposition 8.1).** The naive truncation $e_x \circ e_y = e_{x+y}$ when $x+y \in \Sigma$, else $0$, does *not* satisfy the Jordan identity. Explicit counterexample with residue $e_{-5 - 7\varphi}$.

This is a *complement* to [CCAI 2025], whose direct Jordan algebra on the basis $\{L_x\}_{x \in \Sigma}$ uses the non-trivial quasi-product $x \vdash y = \varphi^2 x - \varphi y$ and works for any convex acceptance window.

## Repository layout

```
scripts/        Python 3.11 sources — algebra layer + experiment drivers + README
certificates/   Mathematical certificates in Markdown
results/        Time-stamped JSON manifests, one folder per run
LICENSE         MIT
README.md       This file
```

## Quickstart

Run the Fibonacci dimension table (Table 7.1 of the manuscript):

```sh
python scripts/run_fibonacci.py --R 10  --M 1 --window symmetric
python scripts/run_fibonacci.py --R 20  --M 1 --window symmetric
python scripts/run_fibonacci.py --R 40  --M 1 --window symmetric
python scripts/run_fibonacci.py --R 80  --M 1 --window symmetric
python scripts/run_fibonacci.py --R 160 --M 1 --window symmetric
python scripts/run_fibonacci.py --R 10 --M 2 --window symmetric --no_instr
python scripts/run_fibonacci.py --R 20 --M 2 --window symmetric --no_instr
```

Reproduce the quasiaddition Jordan-failure counterexample (Proposition 8.1):

```sh
python scripts/verify_quasiaddition_jordan_failure.py --R 8 --window symmetric --max_search 100000
```

For the (companion-paper) Ammann–Beenker experiments:

```sh
python scripts/run_ab.py --R 4 --M 0.5
python scripts/run_ab.py --R 8 --M 0.5
python scripts/run_ab.py --R 4 --M 1 --no_instr
```

Each run writes a JSON manifest to a fresh time-stamped folder under `results/`. The exact runtime environment is **Python 3.11**, **SymPy 1.13**, **NumPy 1.26**, but the algebra layer is pure-Python `fractions.Fraction` and requires no third-party packages.

## Determinism

All randomized sampling uses fixed integer seeds (`0xABCDEF`, `0x123456`, `0xFACADE` in `jordan_tkk.py`); no floating-point arithmetic is used in algebraic-identity verifications; window-membership tests use floats with margins many orders of magnitude larger than double-precision noise. Every `manifest.json` ends with a `self_hash` SHA-256 truncated to 16 hex digits, computed over the manifest dictionary excluding the hash itself, so accidental tampering with archived runs is detectable.

## Citing

If you use this code or these data, please cite the manuscript and this repository:

```bibtex
@misc{TKK-aperiodic-2026,
  author       = {Daniele Corradetti},
  title        = {{TKK-aperiodic}: reproducibility package},
  year         = {2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/DCorradetti/TKK-aperiodic}
}
```

and the prior publication:

```bibtex
@article{CCAI2025,
  author  = {Corradetti, D. and Chester, D. and Aschheim, R. and Irwin, K.},
  title   = {Jordan algebras over icosahedral cut-and-project quasicrystals},
  journal = {J. Geom. Phys.},
  year    = {2025},
  doi     = {10.1016/j.geomphys.2025.105645}
}
```

## License

MIT — see [LICENSE](LICENSE).
