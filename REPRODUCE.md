# Reproduction

Run the Python commands from the repository root. The qualified versions
are Python 3.14.6, SymPy 1.14.0 and mpmath 1.3.0. The commands work in
PowerShell and POSIX shells.

## Verify and replay

```sh
python -B scripts/verify.py
python -m pip install -r companion/requirements.txt
python -B scripts/verify.py --replay
```

The first command uses the standard library to verify the files listed in
[MANIFEST_SHA256.txt](MANIFEST_SHA256.txt). The replay also verifies those hashes,
runs the four existing checker invocations serially and compares each complete
JSON result with its reference bytes. The replay creates a fresh temporary
output directory inside `companion/` and removes it on exit. It runs offline
with the Python packages installed.

| Checker | Qualified result in each Python mode | Reference |
| --- | --- | --- |
| M, adapted support | 97 exact checks and four negative controls. | [Expected M output](companion/evidence/m-normal.json) |
| F, narrow sector | 23 exact algebra or arithmetic rows, five negative controls, one bibliographic row and one input schema and identity row. | [Expected F output](companion/evidence/f-normal.json) |

Expect one manifest summary and, with `--replay`, one PASS line per checker
in ordinary and optimized (`-O`) Python. A failure returns a nonzero exit
status. The [companion guide](companion/README.md) gives individual commands
for retaining output. The finite checks evaluate prescribed identities and
inputs. The paper provides the analytic proofs over the continuous tube.

## Build the manuscript

Use pdfLaTeX and BibTeX from a TeX distribution with the standard `article`
class, T1 encoding support, Latin Modern, `geometry`, `graphicx`, the AMS
math packages and `hyperref`. LaTeX reads the supplied PDFs in
[`paper/figures/`](paper/figures/), so the Python plotting packages are
optional. Run from a fresh checkout.

```sh
cd paper
mkdir build
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build main.tex
bibtex build/main
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build main.tex
```

This creates `paper/build/main.pdf`, a 33-page A4 document. The qualified
build used MiKTeX pdfTeX 1.40.28 and BibTeX 0.99e with
`SOURCE_DATE_EPOCH=1794956400` and MiKTeX's `--disable-installer` option on
each pdfLaTeX pass. All 33 rendered pages matched the distributed PDF.
The build differed in its PDF trailer identifier.
Use the manifest to verify the exact distributed
[manuscript PDF](paper/Nodal-Stability-Radii.pdf).

## Companion provenance

The [reading room](https://www.oleksiybabanskyy.com/cymatics/paper) supplies
[reader companion v1](https://www.oleksiybabanskyy.com/cymatics/nodal-stability-reader-companion-v1.zip).
Its ZIP SHA-256 is
`b6724fa4719bb310d60f1b1d8d193f7d09c3b572440e945768f407ade85cab97`.
All 16 archive members appear byte for byte under
[`companion/`](companion/), including the companion's own
[manifest](companion/SHA256SUMS). The companion guide also describes optional
figure regeneration. The repository keeps the publication figures used for
the TeX build under [`paper/figures/`](paper/figures/).
