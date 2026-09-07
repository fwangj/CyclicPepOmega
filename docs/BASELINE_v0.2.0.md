# CyclicPepOmega v0.2.0 baseline

This release freezes an early, public-data-only architecture for connecting cyclic
peptide chemical modification, observed conformation, molecular recognition, and
future target-selectivity studies.

## Implemented foundation

- PDB and mmCIF ingestion with deposited identifiers and link provenance.
- Candidate detection based on molecular connectivity and peptide atom topology.
- CCD atom/bond information, polymer links, deposited covalent connections, and
  conservative coordinate fallbacks represented as distinct evidence sources.
- Versioned atom-level peptide chemical identities that collapse equivalent models
  and copies while retaining experimental instances.
- Head-to-tail, head-to-side-chain, lactam, disulfide, thioether, and
  multiple-crosslink classifications, plus explicit unresolved/linear outcomes.
- Canonical/noncanonical, atom-level D/L, and backbone N-methyl annotations.
- Backbone phi/psi/omega, radius of gyration, backbone coordinates, and atom-mapped
  RMSD with cyclic-permutation search.
- Same-identity comparison and candidate atom-mapped analogue change reports.
- A 20-case public validation manifest spanning positive, multi-model/multi-copy,
  and difficult negative-control structures.

## Known limitations

Incomplete CCD data, missing coordinates/connections, symmetry and biological
assembly interpretation, depsipeptides, cross-component chemistry, and
protonation/tautomer ambiguity require further work. Public validation covers
selected assertions and is not a PDB-wide accuracy estimate. Candidate analogue
mappings require expert chemical curation.

## Explicitly deferred

IMHB, SASA/polar SASA, ring-shape research, production interface fingerprints,
recognition-rule and selectivity analysis, full-PDB ingestion, databases/APIs/web,
energy models, permeability prediction, ML, docking, and molecular dynamics.
