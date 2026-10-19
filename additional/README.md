# Additional calculations

Run all eight checks with `python run_all.py --additional-only` from the
package root. They use the same dependencies as the notebooks.

| Directory | Files and purpose |
|---|---|
| `analytic/` | The two scripts explicitly named in Appendix B, recovered from `RotemCodes/` in the manuscript archive. Filenames are normalized and invocation instructions updated. Their mathematical algorithms are retained. |
| `pi_integrals/` | Original shift matrices and four regenerating checks: nine polynomial identities, the I initial state, generic uniqueness at the seed, and the K initial state. Failed identity flags now also cause a nonzero exit. |
| `recoveries/` | The supplied exact ζ(5) initial-state readout check. |
| `gregory/` | A new explicit checker of the matrix, cubic, root-bound and threshold calculations printed in Appendix C. |

The script `pi_I_axis_original.py` is data/source loaded by the checks; it
is not a separate verification job. Output JSON files record what was
checked. The complete console logs are in the package-level `verification/`.

The analytic root printouts, finite examples and symbolic identities have
different scopes. See `../COVERAGE.md` before interpreting them as proof
evidence. No additional irrationality claim is inferred from these checks.
