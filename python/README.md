# Section 6 Python notebooks

Nine notebooks reproduce the computational ingredients of Section 6 and
Appendices B.1–B.4 of the supplied manuscript. Each runs independently in a
fresh kernel using only the local helpers and data. Read the package-level
[scope](../COVERAGE.md) and [reproduction guide](../REPRODUCING.md) first.

| Notebook | Stable paper labels | Contents |
|---|---|---|
| [00_guide](00_guide.ipynb) | `thm:intro-zeta34`, `thm:intro-zeta23` | Proof overview, exact final comparisons, statement index |
| [joint/01](joint/01_matrices_sums_and_minors.ipynb) | `lem:z34-limits`, `lem:minors` | Trajectories, pairing, initial states, sum and minor identities |
| [joint/02](joint/02_recurrence_and_divisor.ipynb) | `prop:z34-divisor` | Telescoping, prime counts, exact cells, divisor density |
| [joint/03](joint/03_rates_and_conclusion.ipynb) | `lem:z34-selected-error`, `lem:z34-nonvanishing`, `thm:z34` | Spectrum, readout, nonvanishing, growth, μ < 66.3 |
| [zeta2/01](zeta2/01_zeta2_cmf_and_closed_forms.ipynb) | `lem:z2-limits`, `app:z2-blocks` | Four-block paths, closed forms, symbolic pairing |
| [zeta2/02](zeta2/02_zeta2_prime_arithmetic.ipynb) | `prop:zeta2-int`, `lem:z2-growth` | Carry profile, shell calculations, integral certificate |
| [zeta2/03](zeta2/03_zeta2_rates_and_measure.ipynb) | `lem:z2-selected-error`, `thm:z23-measures` | Rigorous rate enclosures, μ(ζ(2)) < 4.62675 |
| [pure_zeta3/01](pure_zeta3/01_pure_zeta3_approximants_and_rates.ipynb) | `app:zeta3-measure`, `lem:z34-selected-error` | Kernel, partial fractions, sum identification, spectrum |
| [pure_zeta3/02](pure_zeta3/02_pure_zeta3_divisor_and_measure.ipynb) | `app:zeta3-measure`, `thm:z23-measures` | Carry and shell certificates, density, μ(ζ(3)) < 5.4197 |

[PAPER_REFERENCES.md](PAPER_REFERENCES.md) gives the full label map. Printed
numbers are derived from the manuscript and supplied preamble; a compiled
`.aux` for this exact source has not been supplied. Labels are the stable
identifiers. `python paper_references.py --aux path/to/main.aux` can use
compiled reference numbers when available.

Run all notebooks from the package root with
`python python/run_notebooks.py --reverse`, or select one:

```sh
python python/run_notebooks.py joint/02_recurrence_and_divisor.ipynb
```

The runner uses the current Python interpreter and temporary kernel
specifications; no global kernel registration is required. Each notebook
starts in its own directory. Keep notebooks beside their helpers and data.
Executed copies and a report are written to `verification/` within this
Python directory. Source notebooks contain the outputs from the release run.

The earlier single `zeta3_zeta4_one_proof.ipynb` is replaced by the three
`joint/` notebooks. The historical `05_rates_and_linear_independence.ipynb`
maps to `joint/03_rates_and_conclusion.ipynb`. The two analytic scripts named
in Appendix B are included separately in `../additional/analytic/`.

The selected common divisor need not equal the full gcd. The notebooks use
inequalities for the approximation exponent, retain the bounded-gap
nonvanishing statement, and state the analytical and all-index arithmetic
inputs. See [NOTATION_AND_SCOPE.md](NOTATION_AND_SCOPE.md). Older source
provenance is retained as provenance, not as a new Lean verification.
