# Verification of release candidate 1.0.0-rc2

Verified on 19 October 2026 with Python 3.12.14.

**All 9 notebooks (76 code cells) and all 8 additional script checks passed.**
The package was copied to a separate temporary directory with spaces in its
name, outside the project, and executed from an unrelated working directory.
Each notebook used a fresh kernel and the notebook order was reversed.
The complete computational run took 121.87 seconds on the test machine;
runtime depends on the machine. The already installed pinned environment
was reused; a new network installation was not tested.

| Notebook | Code cells | Seconds | Result |
|---|---:|---:|---|
| `zeta2/03_zeta2_rates_and_measure.ipynb` | 7 | 1.50 | Passed |
| `zeta2/02_zeta2_prime_arithmetic.ipynb` | 6 | 14.12 | Passed |
| `zeta2/01_zeta2_cmf_and_closed_forms.ipynb` | 8 | 1.72 | Passed |
| `pure_zeta3/02_pure_zeta3_divisor_and_measure.ipynb` | 7 | 9.15 | Passed |
| `pure_zeta3/01_pure_zeta3_approximants_and_rates.ipynb` | 7 | 9.36 | Passed |
| `joint/03_rates_and_conclusion.ipynb` | 9 | 6.18 | Passed |
| `joint/02_recurrence_and_divisor.ipynb` | 15 | 5.72 | Passed |
| `joint/01_matrices_sums_and_minors.ipynb` | 15 | 30.17 | Passed |
| `00_guide.ipynb` | 2 | 0.67 | Passed |

| Additional check | Seconds | Result |
|---|---:|---|
| `verify_4f3_meijerg_coboundary` | 3.81 | Passed |
| `verify_5f4_meijerg_recessive` | 5.19 | Passed |
| `pi_I_axis_certificate` | 29.51 | Passed |
| `pi_I_initial_check` | 3.01 | Passed |
| `pi_I_uniqueness_check` | 0.30 | Passed |
| `pi_K_seed_check` | 0.35 | Passed |
| `recovery_zeta5_readout_check` | 0.40 | Passed |
| `verify_gregory` | 0.46 | Passed |

The main machine-readable reports are `verification/full_run.json` and
`python/verification/execution_report.json`; per-script logs are in
`verification/`. Source notebooks contain the outputs from this clean-folder
run. Direct dependency versions are in both reports, with the full installed
package inventory in `verification/environment.txt`.

Additional packaging checks covered notebook schemas, Python syntax, local
documentation links, manuscript labels, unchanged manuscript assets, checksum
verification and ZIP integrity. The suggested manuscript patch was tested
with a dry run. Both execution entry points reject disabled assertions.

The Gregory checker certified the positive threshold margin at M=7165 as
8.9121102332795… × 10^−5, with the rigorous enclosure stored in
`additional/gregory/gregory_certificate.json`.

**Limits:** Python verification is not a complete formal proof. Analytical
and uniform arithmetic inputs remain as stated in the paper and notebooks.
Three Lean source snapshots are included and matched to upstream Git blob IDs.
No independent Lean build was performed. The full manuscript was not compiled.
Raw plot data and plotting notebooks were not supplied and are not regenerated.
See `COVERAGE.md` and `MANUSCRIPT_ALIGNMENT.md` for the exact boundaries.

## Combined archive checks

The mathematical Python routines, notebook code cells and saved computation outputs are unchanged
from the tested 1.0.0-rc1 package. The original execution reports are retained.
The rc2 changes add Lean source, build instructions, provenance and updated
reader documentation. The Python runner's report label now explains that Lean is packaged separately. All 256 upstream Lean files and nine configuration
files match their recorded Git blob identifiers. The additional Lean build
driver was exercised in plan mode only. See `verification/combined_package_checks.json` for archive and source-identity checks.
