# Source provenance and curation

The computational source packages are the two ZIPs supplied by the
authors for this task. The combined rc2 archive also includes the three
author-selected Lean repositories at the commits in `lean_upstreams.json`. Their names, sizes and SHA-256 hashes are recorded in
`input_archives.json`. `file_origins.json` maps unchanged and adapted files
to archive members; `source_changes.json` records changes to the Python
companion. Existing subject-level provenance describes earlier ancestry
and is not a claim that those historical sources were rerun here.

## Included

- The Section 6 Python companion, with fresh execution, updated documentation
  and the current manuscript snapshot.
- The active manuscript source, bibliography, preamble and referenced figure
  exports, unchanged.
- Both `RotemCodes/` analytic scripts; names and run instructions normalized.
- The π integral source and regenerating scripts, and the ζ(5) readout check.
- New package infrastructure, documentation, and the Gregory algebra checker.

## Excluded from the current computational route

- The manuscript ZIP's older `certificates/zeta2`, `zeta3`, `zeta34` and
  `revision_checks` sets, which document different or earlier bounds
  (including 4.58989684, 4.64628563 and 5.44019). The supplied Section 6
  companion is the selected implementation for the current stated bounds.
- Higher-zeta, Catalan and Kampé de Fériet certificates tied to other draft
  appendices, and historical recovery inventory tables. They are not needed
  by the selected notebooks and are not represented as current proof evidence.
- The nested Maple archive, which is not required by the selected Python
  reproduction path and has not been executed in this release.
- Unused manuscript drafts, backup files, cache/runtime directories, stale
  verification reports and previous editorial replacement fragments.

Nothing in the original ZIPs or synced project sources was modified.
The computational meaning and the limits of each included check are stated
in `../COVERAGE.md`.

Preparation note: package documentation, orchestration and the additional Gregory
checker were prepared with AI assistance. Their execution is recorded in the
verification reports; this does not change the provenance of the mathematical
proofs in the supplied manuscript.

## Lean copies added in rc2

`lean_file_origins.json` records each unchanged source/configuration copy.
All upstream `.lean` files are included. `file_origins.json` and
`source_changes.json` retain the rc1 computational provenance; rc2 packaging
documentation is described in `../RELEASE_NOTES.md`. Historical source
comments are preserved, while the current guide states what was actually
checked. The original private repositories were not modified.

The GitHub working copy adds the verified personal repository URL to the
README, citation metadata and release notes; computational content is unchanged
from the combined rc2 archive. Its manifests are regenerated for those edits.
