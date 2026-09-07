# Public data model

The central rule is that one peptide is not one structure:

```text
PeptideChemicalIdentity
  -> many ExperimentalPeptideInstances
  -> many CoordinateObservations / Conformations
  -> zero or many InterfaceObservations with Targets
```

## Core entities

- **PeptideChemicalIdentity**: versioned ordered components, atom connectivity,
  cyclization/crosslink topology, stereochemistry, and backbone modifications.
- **ExperimentalPeptideInstance**: a deposited occurrence with PDB, entity, label,
  author, chain, and residue identifiers. Equivalent copies can share an identity.
- **CoordinateObservation**: coordinates for one model/copy, including missing-atom
  status and provenance.
- **Conformation**: descriptors computed from a coordinate observation, with method
  versions and atom selections.
- **ModificationEvent**: an atom-mapped, chemically interpretable change between
  identities; sequence similarity alone is insufficient.
- **MatchedAnaloguePair**: two identities, their atom mapping, curated modification
  events, and comparable conformational observations.
- **Target**: a normalized binding partner, optionally related to a target family.
- **InterfaceObservation**: a conformation/target/complex-specific interface.
- **InteractionObservation**: an atom/residue-level contact with type, geometry,
  method, and provenance.

Identifiers are deterministic only within an explicitly versioned normalization
scheme. Deposited label/auth identifiers are retained rather than replaced. NMR
models and biological-assembly copies must not generate new chemical identities
unless their chemistry actually differs.

Future analogue comparisons store deltas (torsions, RMSD, Rg, IMHB, SASA, polar
exposure) as results linked to source observations, not as properties of the
chemical identity itself.
