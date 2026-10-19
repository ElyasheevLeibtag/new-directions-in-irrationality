# New Directions in Irrationality — Supplementary Material

Elyasheev Leibtag, Hila Barkan, Rotem Kalisch, Shachar Weinbaum,
Michael Shalyt, and Ido Kaminer

Release candidate **1.0.0-rc2**, prepared 19 October 2026.

This package makes the computational ingredients of Section 6 and
Appendices B.1–B.4 reproducible. It also includes selected exact checks for
the π integral construction, a ζ(5) initial-state readout, and the
harmonic-weighted Gregory series. The manuscript source and its referenced
figures are included for identification and comparison.

The Section 6 notebooks recover the numerical conclusions

- μ(Aζ(3)+Bζ(4)) < 66.3 for each fixed nonzero rational pair (A,B);
- μ(ζ(2)) < 4.62675;
- μ(ζ(3)) < 5.4197.

These conclusions use the analytical and uniform arithmetic arguments stated
in the paper. Successful Python execution verifies the documented
computational ingredients; it is not a complete formal verification of the
paper. Exact identities, exhaustive finite certificates and finite examples
are distinguished in the notebooks.

## Lean source revisions and recorded dates

The Lean files are unchanged copies of the following original revisions.
The full commit identifiers and dependency pins are recorded in
[the upstream revision record](provenance/lean_upstreams.json); individual
file hashes are in [the source map](provenance/lean_file_origins.json).

| Original project / public source copy | Original commit | Recorded commit date (UTC) |
|---|---|---|
| [meijerg_zeta3_zeta4](lean/meijerg_zeta3_zeta4/) | `639a790897ec` | 13 October 2026, 19:40:07 |
| [zeta2_four_trajectory_matrices](lean/zeta2_four_trajectory_matrices/) | `a0e8ac417cdd` | 13 October 2026, 08:59:32 |
| [zetas_irr](lean/zetas_irr/) | `d9d52233b153` | 13 October 2026, 08:51:30 |

These dates are the timestamps recorded in the original Git commits. They
identify the selected source revisions; they are not independently certified
proof-completion dates or evidence of public availability on those dates.
Git allows commit dates to be specified by their creator; see the
[Git documentation](https://git-scm.com/docs/git-commit#_commit_information).

This paper repository was first published on **19 October 2026**. The three
original RamanujanMachine repositories and their histories remain private.
This public repository contains attributed source snapshots, not GitHub
forks or imported copies of the full upstream histories. All Lean source
files, Lake configuration files and pinned dependency manifests needed for
the documented build procedures are included. See the [Lean guide](lean/README.md)
for the exact theorem statements, assumptions and build status.

## Reading the supplement

Start with [the notebook guide](python/00_guide.ipynb), then follow the three
subjects in [the notebook map](python/README.md). The notebooks contain saved
outputs and can be read without running them.

| Topic | Entry point |
|---|---|
| Simultaneous ζ(3), ζ(4) approximation | [Matrices, sums and minors](python/joint/01_matrices_sums_and_minors.ipynb) |
| ζ(2) | [CMF and closed forms](python/zeta2/01_zeta2_cmf_and_closed_forms.ipynb) |
| Pure ζ(3) | [Approximants and rates](python/pure_zeta3/01_pure_zeta3_approximants_and_rates.ipynb) |
| Additional exact calculations | [Additional calculations](additional/README.md) |
| Lean sources and proof assumptions | [Lean guide](lean/README.md) |
| Match a paper claim to code | [Coverage and statement map](COVERAGE.md) |
| Conventions and proof inputs | [Notation and scope](python/NOTATION_AND_SCOPE.md) |
| What was actually rerun | [Verification](VERIFICATION.md) |
| Corrections needed in the supplied manuscript | [Manuscript alignment](MANUSCRIPT_ALIGNMENT.md) |

## Reproduce the calculations

Use Python 3.12. From this extracted directory:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_release.py
python run_all.py
```

On Windows, activate with `.venv\Scripts\activate`.
Internet access is needed to install dependencies; the calculations use
local inputs. Do not use `-O`, `-OO` or `PYTHONOPTIMIZE`, which disable
assertions. See [REPRODUCING.md](REPRODUCING.md) for individual commands,
outputs, environments and verification limits.


