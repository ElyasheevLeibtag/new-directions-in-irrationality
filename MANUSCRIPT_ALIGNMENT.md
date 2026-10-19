# Alignment with the supplied manuscript

The manuscript snapshot is unmodified. The companion uses its trajectories,
normalizations and stable labels. The following repairs should be applied to
the publication source. `manuscript/RECOMMENDED_CHANGES.patch` supplies the
concrete wording/notation/rounding repairs described below; it does not alter
the manuscript's theorem statements or certify the supplied Lean developments.

## Numerical and notation repairs

1. In the proof of `lem:z2-growth`, the displayed coarse decimals sum to
   49.098639333..., not a number below 49.098639. Use the appendix's stronger
   carry bound γ > 12.1284279747910. The resulting sum is
   49.0986383585423... < 49.098639. The notebooks recompute the certificate.
2. In the ζ(3) part of `thm:z23-measures`, the displayed expression
   53.287584 + 30 log(3/2) + 57 − 30.458547 equals 91.992990243... .
   Use the appendix's certified combined estimate instead:
   65.451536565866 + 57 − 30.458547243805 = 91.992989322061 < 91.992990.
3. In the ζ(2) integer representatives, replace both `G_nH_n` and
   `D_n^car D_n^sh` by the already defined `mathcal D_n`.
4. The statement δ > 0.015 alone does not imply μ < 66.3. Use the actual
   bounds g > 112.8075 and h < 111.105, giving
   μ ≤ g/(g−h) < 112.8075/(112.8075−111.105) < 66.3.
5. Where h is only an upper bound on the growth of an integer representative,
   write δ ≥ g/h−1, rather than equating δ with g/h−1. A chosen common divisor
   need not remove the entire gcd. State the bounded-gap/nonvanishing inputs
   when invoking the irrationality-measure criterion.
6. The ζ(2) appendix proves g > 62.636575 and then restates only g > 62.63657.
   Retain the stronger bound in that sentence to match the main proof.

These are consistency repairs to the displayed argument. The computational
bounds in this package already use the stronger certified inputs. This is
not a full referee report or a certification of the entire manuscript.

## Supplement references

Use “the supplementary material” together with a description of the
calculation. Exact filenames belong in `COVERAGE.md` and the notebook guide.
The patch replaces the remaining explicit analytic script names and the
generic “supplementary Git repository” wording. It also repairs “scalculated”,
“supraliminal material”, and the duplicated incomplete pure-ζ(3) sentence.

The paper's statement that the appendix contains supplementary code should
distinguish its written proofs from the external computational package.
References to software names such as RamanujanTools and internal figure paths
are retained.

## Items to settle with the final release

- Lean: three exact source snapshots are now included. The source-level
  statement map and PNT assumptions are in `lean/README.md`. The Meijer
  exports give approximation-rate bounds, the ζ(2) export states a measure
  bound, and the independence development has separate main theorems.
  Avoid describing all current numerical measure bounds as formally
  verified: no independent build or full notebook-to-Lean bridge has been
  checked here. State explicitly that the supplied formalizations assume PNT.
- Plots: the source promises plotting data and a notebook in a supplementary
  repository. The archives supply figure exports, not those data or scripts.
  Add them or narrow the availability sentence. The patch deliberately does
  not invent a replacement claim of full plot reproducibility.
- Broad search data: the `ramanujandata` repository and exact revision are
  not included in this package.
- DOI and repository: insert the actual versioned Zenodo DOI and public
  repository URL after creating the records. `CITATION.cff` contains no
  invented identifier.
- The manuscript uses superscript affiliation markers without defining the
  corresponding affiliations or email addresses in the supplied main source.
  Complete those author details in the publication version.

The source and its local references were inspected; the full manuscript was
not compiled and the bibliography was not independently reverified in this
packaging pass.
