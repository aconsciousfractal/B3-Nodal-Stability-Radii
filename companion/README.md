# Reproducibility companion for nodal-face stability

*Reader companion, version 1* · *Sup-norm stability radii for a nodal face
of the B3 spherical discriminant*

The article studies the infimum of the norms of perturbations that destroy
a single internal $C^1$ zero graph on a fixed spherical tube. The perturbation
space is the full real band through degree nine. The norm is the continuous
supremum on the closed tube. In the uniform and adapted exact-radius regimes,
every perturbation in the closed ball of the critical radius preserves the
graph. The failure infimum in those regimes is approached from above.

The [proof and source map](PROOF_SOURCE_MAP.md) locates the arguments in the
article and describes how the finite checks support them.

## Run the exact checks

Extract the ZIP into a fresh directory and open a terminal at its root.
The qualified versions are Python 3.14.6, SymPy 1.14.0 and mpmath 1.3.0.
Use a virtual environment if needed, then install the pins and run:

```sh
python -m pip install -r requirements.txt
mkdir output
python -B scripts/check_example_m.py --output output/m-normal.json
python -O -B scripts/check_example_m.py --output output/m-optimized.json
python -B scripts/check_example_f.py --output output/f-normal.json
python -O -B scripts/check_example_f.py --output output/f-optimized.json
```

The commands work in PowerShell and POSIX shells. The checks run offline
after installation. Run them serially. Each output must be a fresh file.
F resolves relative output paths from the root of the extracted bundle
and requires that its output remain there. Both scripts require an
existing parent directory.

| Check | Expected exit and row counts | Reference output |
| --- | --- | --- |
| M, adapted support | Exit 0. 97 exact and four negative rows. 55 formal complex-monomial identities represent 100 real components. | [M reference](evidence/m-normal.json) |
| F, narrow sector | Exit 0. 23 exact algebra or arithmetic rows, five negative controls, one bibliographic source-boundary row and one input schema/identity row. | [F reference](evidence/f-normal.json) |

Both outputs for each check should match its reference file byte for byte.
For example, run `sha256sum evidence/*.json output/*.json` or, in PowerShell,
`Get-ChildItem evidence/*.json,output/*.json | Get-FileHash -Algorithm SHA256`.
The [manifest](SHA256SUMS) records the digest of each supplied file other
than the manifest itself. Where available, run `sha256sum -c SHA256SUMS`.
The separate ZIP checksum can be compared with the SHA-256 of the download
before extraction.

## Example inputs and scope

Let $L=\omega-\eta$, $\delta=\arctan(\sin\eta\cos\eta/4)$,
$a=\sin\delta$, $B(\theta)=\sin^3(2\theta)\cos(2\theta)/64$ and
$A_\eta=aB(\eta)/9$. Both examples have $A_*=A_\eta$.

- **Adapted support (Example M).** Take $\eta=\arctan(1/100)$,
  $L=21\delta/5$ and $\theta_c=\eta+3\delta/2$. The exact inputs include
  $\tan\delta=25/10001$, $a^2=625/100020626$ and
  $B(\eta)=1249875000/10001^4$. The article proves the radius is $A_\eta$
  even though the prescribed uniform support condition fails. The finite
  check also carries rational estimates for a stronger far-zone margin.
- **Narrow sector (Example F).** Take $\tan\eta=499999999/10^9$ and
  $\tan\omega=1/2$. The input files supply [terminal coefficients and
  parameters](data/terminal.json), [duality constants](data/duality.json)
  and [radius constants](data/radius.json). The article bounds the norm
  over the whole tube and shows an actual failure with cost below
  $A_\eta$, yielding $\rho_\eta\le U<A_\eta$. The exact radius and
  attainment of the failure infimum in this sector remain open.

The bibliographic check validates the [three selected citation
records](data/sources.json) against the supplied [bibliography](references.bib),
including their declared uses and locators. Its coverage is the distributed
citation metadata. The algebra and arithmetic checks support specific steps.
Admissibility, continuous norm bounds, rigidity and the passage from jets
to graph failure follow from the article's proofs. The domain is the signed
spherical field on the specified tube.

## Generate the illustrative figures

For a fresh copy, install the [plotting pins](requirements-figures.txt)
and run from the root of the extracted bundle:

```sh
python -m pip install -r requirements-figures.txt
python -B illustrations/generate_figures.py
```

The script writes `illustrations/tube_support.pdf` and
`illustrations/fibre_threshold.pdf`, replacing files with those names.
The first shows the chart and support at $\theta=\eta$, $t=1/10$ and
chord $3\delta$. The second plots the terminal fibre at amplification
factors $0,1,2$ along $t=\sin r/a$, scaled by $1/(aB(\eta))$.
Both use $\eta=\arctan(1/100)$ and $\omega=\arctan(1/2)$.
The plots illustrate the analytic results. Plotting samples supply
neither a continuous norm bound nor an exhaustive root count.
