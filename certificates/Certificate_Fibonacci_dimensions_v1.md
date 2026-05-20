# Certificate: Fibonacci dimensions and identities v1

## Statement

For the Fibonacci cut-and-project model set $\Sigma\subset\mathbb R$ with symmetric window $\Omega=[-\varphi/2,\varphi/2)$ where $\varphi=(1+\sqrt5)/2$, the propagation-filtered Jordan algebra $J_{\Sigma_R}^{(M)}$ and its TKK Lie algebra $\mathfrak g_{\Sigma_R}^{(M)}$ have the following dimensions for $(R, M)$ pairs:

| $R$ | $M$ | $\|\Sigma_R\|$ | $\dim\mathcal A$ | $\dim J$ | $\dim\Instr$ | $\dim\mathfrak g$ | Jordan id. | JTS id. (sample) | $\mathbb Z_2$ inversion |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 1 | 15 | 5 | 4 | 5 | 13 | PASS (16/16) | PASS (30/30) | PASS |
| 20 | 1 | 29 | 5 | 4 | 5 | 13 | PASS (16/16) | PASS (30/30) | PASS |
| 40 | 1 | 57 | 5 | 4 | 5 | 13 | PASS (16/16) | PASS (30/30) | PASS |
| 80 | 1 | 115 | 5 | 4 | 5 | 13 | PASS (16/16) | PASS (30/30) | PASS |
| 160 | 1 | 231 | 5 | 4 | 5 | 13 | PASS (16/16) | PASS (30/30) | PASS |
| 10 | 2 | 15 | 225 | 120 | 225 | $\sp(30,\mathbb Q)$ | PASS (100/100) | PASS (5/5) | PASS |
| 20 | 2 | 29 | 841 | 435 | 841 | $\sp(58,\mathbb Q)$ | PASS (100/100) | PASS (5/5) | PASS |

## Notation

- $\mathcal A_{\Sigma_R}^{(M)}$ — associative $\ast$-algebra of partial translations on $\ell^2(\Sigma_R)$ of physical propagation $\le M$.
- $J_{\Sigma_R}^{(M)} = (\mathcal A_{\Sigma_R}^{(M)})_{\mathrm{sa}}$ — special Jordan algebra (Theorem 4.1 of `Article_001_v1.tex`).
- $\Instr(J)$ — inner structure algebra, spanned by $V_{a,b}:x\mapsto\{a,b,x\}$ in $\End(J)$.
- $\mathfrak g_{\Sigma_R}^{(M)}=\TKK(J_{\Sigma_R}^{(M)}) = J^- \oplus \Instr(J)\oplus J^+$, $3$-graded Lie algebra.

## Dependencies

- Classical TKK theorem (Jacobson 1968 §I.4; Loos 1975 §1): every special Jordan triple system has a $3$-graded Lie algebra $\TKK(J)$ with brackets $(4.2)$ of `Article_001_v1.tex`.
- Schlottmann uniform density theorem (1998): for regular model sets, $|\Sigma_R\setminus\Sigma_R^{\mathrm{int}}|/|\Sigma_R|\to 0$.

## Proof

For $M=1<\varphi$:
- Generators are $\{T_{-1}, T_0=I, T_1\}$ (since the smallest physical gap in $\Sigma$ is $1$ and the next gap is $\varphi>1$, so no other lattice vector has physical norm $\le 1$).
- $T_1^2 = 0$ (no two consecutive $S=1$ gaps in Fibonacci substitution), similarly $T_{-1}^2 = 0$.
- $T_1 T_{-1}$, $T_{-1} T_1$ are nonzero diagonal projections.
- Hence $\mathcal A = \mathrm{span}\{I, T_1, T_{-1}, T_1 T_{-1}, T_{-1} T_1\}$, dimension 5.
- $J = \mathcal A_{\mathrm{sa}}$ has basis $\{I, T_1 T_{-1}, T_{-1} T_1, (T_1+T_{-1})/2\}$, dimension 4.
- $\Instr(J)$ has dimension 5 by direct exact computation (run on each patch, see Computational support below).

For $M=2>\varphi$:
- Generators include $T_{\pm 1}$ and $T_{\pm\varphi}$, sufficient to reach every pair of patch points by composition.
- Hence $\mathcal A = M_n(\mathbb Q)$, full matrix algebra; $J = \Sym_n(\mathbb Q)$; $\TKK(J) = \sp(2n,\mathbb Q)$ by Remark 3.1 of `Article_001_v1.tex` (classical: $\TKK(\Sym_n) = \sp(2n)$, see Jacobson 1968 §I.4 and Loos 1975 §1).

## Verification of hypotheses

- Regularity of $\Omega=[-\varphi/2,\varphi/2)$: compact, $\overline{\Omega^\circ}=[-\varphi/2,\varphi/2]$ (closure does not change Haar-measure as boundary is two points); $\partial\Omega$ has Haar measure zero. ✓
- $\pi_{\mathrm{phys}}|_\Lambda$ injective: $(a, b)\mapsto a + b\varphi$ injective on $\mathbb Z^2$ since $\varphi$ irrational. ✓
- $\pi_{\mathrm{int}}(\Lambda)$ dense in $G=\mathbb R$: $\{a + b\varphi': a,b\in\mathbb Z\}$ is dense in $\mathbb R$ for $\varphi'$ irrational. ✓

## Computational support

For each $(R, M)$ row, run:

```
python "2. Scripts/Article_001/run_fibonacci.py" --R <R> --M <M> --window symmetric
```

This produces a JSON manifest at `7. Results/Article_001/run_fibonacci_R<R>_M<M>_<timestamp>/manifest.json` containing all reported dimensions, identity verifications, and runtime data. The output is reproducible; the package pin is Python 3.11, SymPy 1.13, NumPy 1.26.

Representative manifest fields:

```json
{
  "model_set": "Fibonacci",
  "window": "[-phi/2, phi/2)",
  "patch_radius_R": 80,
  "propagation_max_M": 1,
  "n_points": 115,
  "dim_A": 5,
  "dim_J": 4,
  "dim_instr": 5,
  "dim_TKK": 13,
  "jordan_identity": {"verified": true, "n_checked": 16, "n_failures": 0},
  "jordan_triple_identity_sample": {"verified": true, "n_checked": 30, "n_failures": 0},
  "z2_inversion_symmetry": {"perm_total": true, "jordan_invariance": true}
}
```

The manifest set under `7. Results/Article_001/` was generated on 2026-05-19.

## Failure modes / limitations

- The "Jordan identity verified" status checks all $4^2=16$ basis pairs for $M=1$ and a random sample of 100 pairs for $M=2$. For special Jordan algebras the identity is a theorem (Jacobson 1968), so any failure here would be an implementation bug; the verification is therefore an implementation sanity check, not a mathematical claim.
- "Jordan triple identity (sample)" checks a random sample of 30 (resp.\ 5) quintuples. Same caveat applies.
- The Z_2 inversion symmetry verification for $R=160$ presumes the patch is exactly symmetric around 0 in $\Sigma$. This is guaranteed by the symmetric window choice and is verified in the manifest (`perm_total: true`).
- The saturation result for $M=2$ relies on the path-reachability argument in the proof of Theorem 6.2 of `Article_001_v1.tex`. The proof is direct but uses substitution-rule combinatorics on Fibonacci paths; a more general statement for arbitrary regular model sets is left to future work.

## Public-manuscript status

Included in the main draft `Article_001_v1.tex`:
- Theorem 6.1 (stability at $M<\varphi$) — proof in draft §6.
- Theorem 6.2 (saturation at $M\ge\varphi$) — proof in draft §6.
- Remark 6.3 (interpretation as aperiodic core) — in draft §6.
- Table 7.1 (dimension table) — in draft §7.

## Changelog

- v1 (2026-05-19): initial version, certifies $(R, M) \in \{(10,1), (20,1), (40,1), (80,1), (160,1), (10,2), (20,2)\}$.
