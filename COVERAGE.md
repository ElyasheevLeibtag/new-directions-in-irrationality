# Coverage and manuscript-to-code map

The reference manuscript is `manuscript/main.tex`, copied unchanged from the
received focused archive. Its exact LaTeX labels are the primary identifiers.
The package covers the computations below. It does not claim to reproduce
every experiment or formally verify every theorem in the paper.

| Paper material | Delivered implementation | Verification scope |
|---|---|---|
| `thm:intro-zeta34`, `subsec:z34`, `app:zeta34-proof` | `python/joint/01_…`, `02_…`, `03_…`; `additional/analytic/verify_5f4_meijerg_recessive.py` | Symbolic matrix/readout identities, telescoping, finite exact certificates, rate enclosures and final bound, with the stated proof inputs |
| `thm:intro-zeta23`, `app:z2-blocks`, `lem:z2-growth` | `python/zeta2/01_…`, `02_…`, `03_…`; `additional/analytic/verify_4f3_meijerg_coboundary.py` | Four-block formulas, pairing, carries/shells, exact rate enclosures and final ζ(2) bound |
| `app:zeta3-measure` | `python/pure_zeta3/01_…`, `02_…` | Identification with joint sums, finite coefficient diagnostics, spectral and divisor certificates, final ζ(3) bound |
| Integral CMF for π, `app:pi-build-recovery` | `additional/pi_integrals/` | Nine symbolic shift identities, initial states and generic uniqueness at a nondegenerate seed |
| ζ(5) initial-state recovery | `additional/recoveries/recovery_zeta5_readout_check.py` | One exact initial-state/readout identity; no ζ(5) irrationality theorem |
| `app:gregory-series`, `eq:hh-T`, `eq:hh-detT`, `eq:hh-Tinf`, `eq:hh-cyclic`, `lem:hh-roots`, `prop:hh-gap` | `additional/gregory/verify_gregory.py` | Printed matrix identities, cubic/root-bound algebra, rigorous positive threshold margin at M=7165 |
| Equal-accuracy plots and Gregory experimental exponent plot | Rendered figures in `manuscript/figures/` | Exported figures only; raw sampled data and plotting notebooks were not supplied |
| Broad CMF trajectory search / `ramanujandata` | Not included | Repository and exact search-data revision still needed for these experiments |
| Lean formalizations | Three pinned source copies; [Lean guide](lean/README.md) | Source identity checked; theorem scope and PNT assumptions documented; no independent Lean build or complete notebook-to-Lean bridge claimed |

The additional analytic scripts print numerical root approximations as well
as exact identities. Their decimal printouts are not interval certificates;
the rigorous logarithm/root enclosures for the claimed Section 6 bounds are
in the corresponding Python notebooks.

The Gregory checker was added during preparation of this release by
transcribing the displayed formulas and verifying their algebraic
relationships. It does not independently prove the trajectory construction,
recessiveness, nonvanishing for every index, integrality, or other literature
criteria. These arguments remain in the manuscript. It does not produce
the experimental fit or its sample points.

Several certificates embedded in the manuscript ZIP concern older
trajectories or numerical bounds. They are excluded from the current
computational package; see `provenance/README.md`. The complete original
ZIP remains identified by its checksum rather than being presented as a
second conflicting set of current certificates.
