# Notation and verification scope

This companion follows the supplied manuscript section, not the older notation or
theorem numbering in the archived supporting articles. It contains Python
only. Historical provenance records mention earlier Lean sources. The combined package includes pinned Lean source copies in `../lean/`, with assumptions and build status documented in `../lean/README.md`. These notebooks do not execute Lean.

## Common notation

| Paper notation | Meaning in the notebooks |
|---|---|
| `M(n)` | Forward trajectory matrix, with row states multiplied on the right. |
| `Q(n)` | Special-function trajectory, or its relevant reduced block. |
| `Phi_M(n)` | `M(0) ... M(n-1)`, with the empty product equal to the identity. |
| `Phi_Q^vee(n)` | Inverse transpose of the entire forward product. |
| `K(n), c(n)` | Dual pairing: `M(n) K(n+1) = c(n) K(n) Q(n)^(-T)`. |
| `q_n, p_{s,n}` | Rational readouts of the forward walk. |
| **q**`_n`, **p**`_{s,n}` | Joint-proof readouts after multiplication by `C_n/12`. |
| `d_m` | `lcm(1,...,m)`, with `d_0=1`. |
| `D_n` | The joint normalization `6 d_(19n)^4`. |
| `mathcal D_n` | A proved common divisor, not necessarily the full gcd. |
| `g` | Exponential rate of decay of the rational approximation error. |
| `h` | Upper bound for exponential growth of an integer denominator representative. |
| `delta` | Approximation exponent of the reduced rational sequence. |
| `mu` | Irrationality measure of the target number. |

## Conventions that need care

1. **Rates give an inequality.** From an error upper bound and a denominator
   growth upper bound, the conclusion is `delta >= g/h - 1` on the relevant
   nonzero-error sequence, rather than an asserted equality. A chosen common
   divisor need not be the full gcd. The irrationality-measure conclusion also
   uses the nondegeneracy / bounded-gap hypotheses of the paper's criterion.
   Numerical examples alone do not establish that criterion.

2. **Joint nonvanishing.** Lemma 6.15 (`lem:z34-nonvanishing`) says that three
   consecutive errors cannot all be zero. It does not say that each error is
   nonzero. The notebooks retain this precise formulation.

3. **Joint final rounding.** The weaker displayed statement `delta > 0.015`
   by itself does not imply `mu < 66.3`. Use the actual rate bounds:
   `g > 112.8075`, `h < 111.105`, giving
   `mu <= 112.8075/(112.8075-111.105) < 66.3`.
   The guide verifies this as an exact rational inequality.

4. **The zeta(2) divisors.** Use `mathcal D_n^car` for the carry divisor,
   `mathcal D_n^sh` for the shell divisor, and `mathcal D_n` for their product.
   Reserve `G_n` for the Meijer function. The unmodified manuscript snapshot
   retains its earlier aliases; the notebooks use
   the notation of the preceding manuscript lemma.

5. **Coordinates and indices.** Some source code orders the zeta(2) upper
   parameters as `(x_1,x_2,x_3,x_0)`. Permuting those upper parameters leaves
   the hypergeometric function unchanged. The pure zeta(3) archive uses a
   native Meijer-G walk depth `N` and kernel index `n=N-1`. Its notebooks
   explicitly distinguish the native walk from the kernel formula in
   Appendix B.4.

6. **Printed reference numbers.** The notebooks use the supplied manuscript source
   with the shared, section-local theorem counter from the existing project
   preamble. A compiled `.aux` for that exact draft was not supplied. Both
   numbers and exact LaTeX labels are included; the labels remain unambiguous
   if the paper is renumbered. See `PAPER_REFERENCES.md`.

7. **Pure zeta(3) intermediate rounding.** The latest manuscript's chain
   `53.287584 + 30 log(3/2) + 57 - 30.458547` is approximately
   `91.9929902432449`, slightly above `91.992990`. The earlier truncation
   `12.163953` is also not an upper bound for `30 log(3/2)`. The notebooks
   certify the combined rate `log(Q_n C'_n)/n < 65.451536565866` and use
   the divisor density `>30.45854724380549...`. Consequently
   `h < 65.451536565866 + 57 - 30.45854724380549... < 91.992990`.
   This preserves the final claimed bound while repairing the intermediate
   rounding direction. Neither the paper snapshots nor the Lean files have
   been changed.

8. **Zeta(2) intermediate rounding.** The three displayed decimals alone give
   `38.168733 + 25 - (12.128427 + 233/120) = 49.098639333...`, which is
   slightly larger than `49.098639`. The exact carry-integral certificate
   supplies the stronger estimate `25 - gamma - 233/120 < 10.929906`.
   Combining that with `38.168733` correctly gives `h < 49.098639`.
   The notebooks use this certified route, preserving `mu(zeta(2)) < 4.62675`.

## What a successful notebook run means

Each notebook identifies whether a cell regenerates a symbolic identity,
recomputes an exact finite certificate, tests finitely many examples, or
combines previously established mathematical inputs. Arbitrary precision
approximations are not described as rigorous enclosures unless interval
arithmetic or explicit rational upper/lower bounds are used. Archived output
files are not substitutes for recomputation.

The asymptotic arguments using recurrence theory, Mellin–Barnes integrals,
uniform prime divisibility and the prime number theorem remain mathematical
arguments in the paper / supplied formalization. The Python notebooks make
their computational ingredients and numerical consequences reproducible.

## Identification of the pure and joint sums

The pure notebooks now include the all-index gamma-reflection argument:
`Cprime_n (Q_n,P3_n) = (-1)^n (bold_q_n,bold_p3_n)`.
There is no index shift between these two sum presentations. The native
Meijer walk convention remains separate. An independent binomial calculation
checks every Laurent coefficient and both sums for n=0 through 6, including
the simple-pole tails; it supplements the proof rather than replacing it.
