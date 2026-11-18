# Sup-norm stability radii for a nodal face of the $B_3$ spherical discriminant

**Oleksiy Babanskyy · 2026 · Repository package 1.0.0**

[Manuscript PDF](paper/Nodal-Stability-Radii.pdf) · [LaTeX source](paper/main.tex) · [Reading room](https://www.oleksiybabanskyy.com/cymatics/paper) · [Cymatics Studio](https://www.oleksiybabanskyy.com/cymatics/)

How large can a real perturbation be before a truncated nodal line of the
$B_3$ spherical discriminant ceases to be a single interior $C^1$ graph?
The perturbations span all 100 real spherical polynomial modes of degree at
most nine. Their size is the supremum norm on the entire closed spherical
tube.

1. **Variational radius.** The graph-failure radius equals the minimum cost of
   simultaneously cancelling the field and its normal derivative at some point of the closed tube.
   The paper bounds this radius throughout the admissible parameter domain.
2. **Exact regimes.** An explicit spherical mixed-jet certificate gives a
   sharp radius in the uniform and adapted support regimes. The unique
   critical perturbation is a scaled degree-nine Chebyshev polynomial.
   Every perturbation in the closed radius ball preserves the graph, including
   the critical field with a singular zero level. Any larger amplification of
   that field gives an actual failure.
3. **Limits and conditioning.** Example M extends the exact formula beyond
   the uniform support condition. Example F proves an upper bound strictly
   below the candidate formula in a narrow sector. The exact radius there
   remains open. In either exact regime, a family strictly inside the ball
   has uniformly converging graphs and terminal slopes diverging as
   $\varepsilon^{-1/3}$ at minimizing terminals.

The results concern a signed real field within the specified band, tube and
normalization. Global nodal topology, plate response, particle transport and
material deposition lie outside this result. The [companion](companion/README.md)
replays selected symbolic and rational checks for examples M and F. The
[proof and source map](companion/PROOF_SOURCE_MAP.md) links these checks to the
analytic arguments for admissibility, continuous norm bounds and graph failure.

## Reproduce

From the repository root, use Python 3.14 with the pinned dependencies to
verify the distributed bytes and run all four replays.

```sh
python -m pip install -r companion/requirements.txt
python -B scripts/verify.py --replay
```

[Reproduction instructions](REPRODUCE.md) give the manuscript build, the
qualified tool versions and the expected results.

## Cite and reuse

The manuscript records the use of AI assistance in preparing the paper and
its companion. Citation metadata is in [CITATION.cff](CITATION.cff). The
original manuscript is CC BY 4.0. The original companion, scripts and
supporting repository material are MIT. [License scope](LICENSE_SCOPE.md)

The [optional GitHub Actions workflow](ci/verify.yml) runs the same four checks on Python 3.14. To enable it in a fork, copy it to the workflow directory at `.github/workflows/verify.yml`.
