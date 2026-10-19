Exact pi integral certificates

All scripts require only Python 3 and SymPy. No network access, RT, PSLQ or external files.

From the package root, run:
 python run_all.py --additional-only
Or enter additional/pi_integrals/ and run the individual scripts.

The axis file was extracted from bestDeltas.ipynb (historical source notebook; not a runtime dependency), zero-based cell22 (identical cell62), then changed only to import Matrix from SymPy instead of RT.
The original notebook uses the readout in cells25,66. Its cells64–65 use the plotting seed(1,2,1) and ray(2,3,2).

The other scripts were written independently for the 2026-10-07 audit. They check nine generic rational integration-by-parts identities, the exact three initial integrals, generic uniqueness of the coefficient specification, and the exact Bai-origin seed readout in the I basis. Each script writes its accompanying JSON.

Scope: these certificates establish contiguity and initial coefficient identities. They do not prove path admissibility, convergence or denominator-growth estimates for arbitrary rays. In particular the plotting seed is a b+10 contiguous variant of the odd-subsequence of the published Zeilberger–Zudilin integral, not the literal published sequence.
