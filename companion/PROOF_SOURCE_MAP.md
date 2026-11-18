# Proof and source map

This companion accompanies *Sup-norm stability radii for a nodal face of
the B3 spherical discriminant*. The table uses the article's LaTeX labels
and, where helpful, theorem numbers.

| Claim | Article locators | Finite support |
| --- | --- | --- |
| The real band has dimension 100 and the fixed closed tube defines a norm. | `eq:target`, `eq:domain`, `eq:chart`, `prop:basis`, `app:harmonics`, `lem:norm` | F checks the dimension used for fixed-point duality and counts the dimension of an auxiliary odd subspace. |
| Failure radius equals the minimum exact-zero-jet cost for this target. | `thm:main`, `lem:graph`, `sec:reduction`, `prop:exactjet`, `prop:jetonly` | F checks a finite jet correction identity. Its cubic negative control illustrates a degenerate jet on a smooth zero graph. |
| Fixed-point duality and the fibre inequality bound the whole-tube radius. | `prop:dual`, `eq:dualmax`, `prop:barrier`, `prop:strongfloor`, `eq:universalradius` | F checks sign and dimension formulas. M20–M25 check selected scalar and polynomial factors. The interval inequality applies in the open normal interval. The boundary proof uses the value constraint. |
| The moving atoms represent the scalar functional $J=(a/c_t)E_r+\lambda(t)E$ throughout the full band. The contact signs attain the norm. | Theorem 5.1, `eq:supportcondition`, `eq:roots`, `eq:atoms`, `eq:fourierbasis`, `eq:rationalinterp`, `thm:certificate` | M checks 55 formal moment identities and scalar factors. The proof places physical contacts and handles the strip tangency. |
| At supported central jets the value-and-derivative cancellation has a unique minimizing field. The closed ball preserves the graph under the exact-radius hypotheses. | `lem:rigidity`, `cor:pointlower`, `eq:Lcondition`, `eq:Mcondition`, `eq:Mright`, `thm:radius`, `prop:criticalgraph`, `sec:sharp` | M checks rational envelopes for the adapted example. The proof gives preservation at equality and actual failure along a sequence of amplifications converging to that threshold. |
| The adapted radius can hold beyond the uniform prescription. The narrow sector can have a strict upper bound below the terminal formula. | `ex:M`, `eq:Mexample`, `ex:F`, `eq:Fexample`, `eq:Fglobal`, `eq:Fterminalbound`, `eq:Ftubebound`, `eq:Fupper` | M checks the prescribed uniform support failure and adapted arithmetic. F checks narrow-sector constants and witness formulas. The radius proof uses the whole-tube norm and a continuous three-zero argument. |
| A strict norm gap below threshold controls the normal derivative and the graph. Slopes can diverge as a family approaches a limit with a smooth zero graph. | `prop:Dlip`, `prop:margin`, `cor:budget`, `thm:degeneration`, `eq:slopedegeneration`, `sec:limits` | These conclusions use analytic estimates in the article. |

## Literature

The [bibliography](references.bib) gives the article references. The
[source records](data/sources.json) select three locators used to
orient the companion.

- Erdélyi, Theorem 2.1, supplies the trigonometric value and derivative
  inequality applied to the normal fibre. [Publication](https://doi.org/10.1007/s00365-013-9208-0)
- Dressler et al., §§4.2–4.3, give classical multivariate extremal-polynomial
  context, including representation and sign complementarity. [Publication](https://doi.org/10.1137/24M1717853)
- Davydov and Schaback, Theorem 9, give fixed-node duality for minimal exact
  numerical differentiation formulas. [Publication](https://doi.org/10.1007/s00211-018-0973-3)
