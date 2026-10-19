# ζ(2) Python companion

These three notebooks follow the notation and numerical thresholds of the user-supplied paper. Each runs in a fresh kernel without running another notebook. Keep this directory intact: `companion.py` and `certificates/` are local dependencies, not outputs of a previous notebook.

| Notebook | Paper statements | What it recomputes |
|---|---|---|
| `01_zeta2_cmf_and_closed_forms.ipynb` | Theorem 6.2; Lemmas 6.17–6.18; Appendix `app:z2-blocks` | Paper-order CMF moves, two four-block paths, RamanujanTools comparison, closed-form readouts, finite generator identities, full symbolic dual pairing and initial value |
| `02_zeta2_prime_arithmetic.ipynb` | Lemmas 6.20–6.21; Appendix `app:z2-blocks` | Carry valuations, all eight shell subregions, exhaustive 390-cell arrangement certificate, supplied congruence diagnostics, integer pairs, exact carry integral |
| `03_zeta2_rates_and_measure.ipynb` | Lemmas 6.19 and 6.21; Theorems 6.2 and 6.16 | Forward limiting spectrum, exact logarithm enclosures, arithmetic growth, δ>0.2757 and μ(ζ(2))<4.62675 |

Printed numbers refer to the supplied manuscript; the stable LaTeX labels are included in the notebook Markdown as well.

## Run

Install the bundle's root requirements and use that environment's Python kernel. Open any notebook and run all cells from top to bottom. The working directory may be either this `zeta2` directory or the bundle root. The notebooks need no network access after installation and never import a notebook or a Lean file.

All default certificate calculations are enabled. The supplied Python certificates regenerate their JSON files under `certificates/`; their full console logs go into `generated/`. A failed assertion raises an error. Precomputed JSON is not accepted in place of running a certificate.

## Notation and interpretation

- The public helper uses the paper's upper order `(x0,x1,x2,x3)=(0,1,1,1)`. Original certificate modules retain their source order `(x1,x2,x3,x0)=(1,1,1,0)`. Upper-parameter symmetry identifies these trajectories; notebook 01 makes the conversion explicit.
- Row states multiply matrices on the right. Helper readouts return `(p_n,q_n)`; `READOUT_ROWS` puts the numerator row first.
- `C_n=∏ binomial(un,8n)`, with `u=9,10,11`, agrees with both the main text and appendix.
- The carry and shell divisors are `mathcal D_n^car` and `mathcal D_n^sh`, with product `mathcal D_n`. Reserve `G_n` for the Meijer function.
- The integer representatives need not be coprime. Their denominator is an upper bound for the reduced denominator.
- The limiting forward matrix is called `A_infinity`, not `Q_infinity`. Its spectral modulus ratio gives the same error rate under the paper's dual pairing.

The reproduced thresholds are `g>62.636575`, `h<49.098639`, `g-h>13.537936`, `δ>0.2757`, and `μ(ζ(2))<4.62675`. All threshold assertions use exact rational bounds. Decimal displays are informational.

## Evidence scope

Exact finite matrix and readout comparisons are diagnostics, not proofs of all-depth statements. The floor-arrangement and scalar-rate certificates are exhaustive finite computations with scope explained in notebook Markdown. The prime number theorem, Barnes-integral estimates, recurrence asymptotics, all-depth divisibility, and the irrationality-measure criterion are the paper's mathematical inputs.

The full symbolic pairing of Lemma 6.17 is reconstructed from the paper parameters in `dual_pairing.py` and verified exactly by notebooks 01 and 03. The original cited script is also supplied in `../../additional/analytic/`, with its own execution log. The source certificates' legacy `NUMERICAL_CANDIDATE` status fields are preserved without treating them as proof: each notebook reruns the source and asserts the needed exact inequalities.

`SOURCE_MANIFEST.json` records original paths and SHA-256 hashes for the unchanged supplied Python certificate modules. `companion.py` is a new documented Python interface in the paper's notation; `dual_pairing.py` derives the missing symbolic connection certificate from the paper's defining equations. These Python notebooks do not create or execute Lean. The combined supplement includes unchanged Lean source snapshots; see `../../lean/README.md`.

## Manuscript rounding correction

The proof of Lemma 6.21 displays `38.168733 + 25 − (12.128427 + 233/120) < 49.098639`. Its left side is actually `49.098639333…`. The regenerated exact lower bound `gamma > 12.1284279747910` fixes this last rounding step and proves the stated `h < 49.098639`; the final measure bound is unchanged. Notebook 03 records this explicitly.
