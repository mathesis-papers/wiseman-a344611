# Reverse-alternating partition sums

Haobo Ma — ChronoAI Pte Ltd, Singapore  
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore

Revision 2, November 5, 2026. License: CC BY 4.0 (`cc-by-4.0`).

## Statement and sources

The conjecture in OEIS A344611, revision 27, FORMULA, is
`A344611(n) <= A160786(n)`. The note proves it for every `n >= 0`,
with equality exactly for `n = 0, 1, 2, 3`. The general injection is
`C(m) <= O(m+1)`, where `C` counts partitions with nonnegative
reverse-alternating sum and `O` counts partitions with odd length.

The note also establishes both unrestricted monotonicity questions and their
exact plateaus, quantitative deficits for two comparisons, and a signed-endpoint
convolution extending Manyama's (2026) formula. Ordered-composition formulas
recover those recorded by Barry (2007, 2009) and Israel (2018), with exact index
shifts, and imply Mathar's (2014) recurrence. Source descriptions of distinct
negative-sum partitions and a path interpretation receive explicit corrections.
The unspecified right-hand argument in A344743 is retained as unspecified;
its common-index reading fails at `n=1`, while a shifted comparison is proved.

The exact public source locators and bibliography are in `note.tex` and
`references.json`. The strict odd-length product formula is due to Jovovic
(2002) and also appears in Andrews–El Bachraoui (2026), Theorem 5.

## Build

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, use a LaTeX installation providing `geometry`, `amsmath`,
`amssymb`, `amsthm`, `booktabs`, `array`, `fontenc`, `lmodern`, and `hyperref`:

```sh
pdflatex -no-shell-escape note.tex
pdflatex -no-shell-escape note.tex
```

The bibliography is included in the TeX source; BibTeX is unnecessary.

## Exact finite checks

Python 3.9 or later, standard library only:

```sh
python3 verify.py --max-n 20 --output reproduced.json
```

The expected status is `ALL_CHECKS_OK`. The program enumerates all 313,065
partitions of weights 0–42 and checks both injections through source weight 41.
Ordinary and distinct-part length-parity dynamic programs provide separate
count calculations. It checks 1,849 signed-endpoint convolution coefficients,
both exact plateau sets, the quantitative deficits, and the stated source
counterexamples. It also enumerates all 65,536 compositions of weights 0–16,
checks their full alternating-sum distributions, and checks the ordered
recurrence at indices 3–100. `verification.json` contains the expected counts.
These finite checks diagnose the formulas; the unbounded statements have
complete written proofs in the note.

## Lean declaration scope

At commit `d133fcb1d5e4878bdf5232fd0c07ef2abb3da1d9`, the declaration
`D5.S1.Words.Compositions.WisemanA344611OddLengthInjection.result` has conclusion
`forall n, C(2*n) <= O(2*n+1)`. Both this file and its recurrence dependency
are linked directly in the PDF. `trimmedSums` realizes the reverse-alternating
sum by `s_0=0`, `s_i=lambda_i-s_(i-1)`. The other declarations described there
concern distinct recurrence endpoints, not the sign comparison. General
weights, strictness, plateaus, signed convolutions, ordered compositions, and
corrections are written proofs, with the finite checks specified above.

## Files

`note.tex` is the sole main source. `metadata.json` gives the scientific
description and byline. `references.json` supplies exact source identities.
`verify.py` and `verification.json` supply finite checks. `LICENSE.txt` gives
the license. `MANIFEST.json` binds the source archive's member bytes.
