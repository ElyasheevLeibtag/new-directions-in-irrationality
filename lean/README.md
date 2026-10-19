# Lean formalizations

This directory contains public source snapshots of three RamanujanMachine
projects selected by the authors. All **256 Lean source files**, the original
Lake configurations and dependency manifests are included without changes.
These are source copies at exact revisions, not GitHub forks. The upstream
repositories were private when copied; no access to them is needed to read
these snapshots. The dependencies named in the manifests are public.

**Verification status:** source identity and the statements below were checked.
The projects have not been independently compiled during preparation of this
release. A successful build would still leave the explicitly assumed prime
number theorem described below. Python results do not certify the Lean proofs.

## Which project proves which statement?

| Project | Source entry point | Statement in the source |
|---|---|---|
| ζ(3), ζ(4) approximation rates (70 files) | [Blueprint_5255.lean](meijerg_zeta3_zeta4/lean/Meijer5255/Blueprint_5255.lean) | `Meijer5255.main_joint`: liminf approximation exponent > 0.01634832; `Meijer5255.main_pure3`: > 0.22626234. These exports are approximation-rate statements. They do not themselves state irrationality or irrationality-measure bounds. |
| ζ(2) irrationality measure (45 files) | [MainTheorem.lean](zeta2_four_trajectory_matrices/Zeta2FourTrajectoryMatrices/MainTheorem.lean) | `zeta2_irrationality_measure`: μ(ζ(2)) ≤ 62636575/13537936 < 4.62675. |
| Linear independence (141 files) | [blueprint.lean](zetas_irr/blueprint.lean), [Blueprint_pfq54.lean](zetas_irr/pfq54_formalization_handoff/Blueprint_pfq54.lean) | `ZetasIrr.theorem_1_6`: linear independence of 1, ζ(3), ζ(4) over ℚ. `Pfq54.main_theorem` also states independence and an eventual approximation-exponent bound ≥ 1/375. |

The first project's definition assigns an infinite exponent to an exact
approximant. Its positive exponent conclusion alone therefore does not prove
irrationality. The independence project has a separate irrationality argument.
The older 1/375 rate is distinct from the current Section 6 notebook rates.
No formal bridge equating every notebook normalization and divisor with these
independent Lean developments has been checked in this packaging work.

## Explicit prime-number-theorem assumptions

- [Meijer5255/PNT.lean](meijerg_zeta3_zeta4/lean/Meijer5255/PNT.lean), line 18:
  one `sorry` in `Meijer5255.PNT.prime_number_theorem`.
- [arith_rate.lean](zeta2_four_trajectory_matrices/Zeta2FourTrajectoryMatrices/ArithmeticRate/arith_rate.lean), line 78:
  `axiom Zeta2ArithRate.prime_number_theorem`; no `sorry` was found in this project's source.
- [log_Delta_le.lean](zetas_irr/blueprints/log_Delta_le.lean), line 133:
  one `sorry` in the prime-number-theorem declaration.

All three assume the asymptotic formula θ(x)/x → 1. These are explicitly
unproved inputs in the supplied formalizations. The static scan in
[lean_source_audit.json](../verification/lean_source_audit.json) describes
source tokens; it is not independently reproduced `#print axioms` output.

## Build and inspect the actual theorem entry points

Install Lean via elan, use Python 3.12, and retain each `lake-manifest.json`.
All projects specify Lean `v4.35.0-rc4`, but use different pinned mathlib
revisions. Do not combine their dependency directories or run `lake update`
to replace the supplied pins.

From the supplement root, inspect the commands without running a build:

```sh
python lean/check_projects.py --plan
```

To fetch dependencies, build all declared libraries and check the substantive
entry points, including their printed axioms:

```sh
python lean/check_projects.py
```

This can take substantial time and disk space. Network access is required
for the pinned public dependencies and cached mathlib artifacts. Logs and a
machine-readable result are written to `verification/lean_build/`.

For `zetas_irr`, the default `lake build` only selects the small `ZetasIrr`
stub. The driver explicitly builds all 137 declared libraries and checks
both substantive blueprint files. This avoids mistaking the default build
for a check of the main theorems. The driver itself has only been checked in
plan mode here; no successful Lean build is claimed.

## Provenance and scope of the copies

Exact commits, toolchains and mathlib pins are in
[lean_upstreams.json](../provenance/lean_upstreams.json).
[lean_file_origins.json](../provenance/lean_file_origins.json) maps each file
to its original path, Git blob identifier and SHA-256 hash.

Only Lean source and build configuration files are copied. Historical
handoff notes, earlier manuscript versions, unrelated computational work,
logs and prebuilt binaries are omitted. Some preserved source comments refer
to those upstream documents or earlier development stages; the current
entry points and verification status are described above. Upstream authors
and attribution in the source are retained. No license file was found in
the selected upstream trees; see [licensing status](../LICENSE_STATUS.md).
