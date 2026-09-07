# Architecture

The initial pipeline has four deliberately separate layers:

1. `io/structures.py` converts PDB or mmCIF records into a format-neutral structure.
2. `core/schema.py` preserves atoms, author residue IDs, entity IDs, and covalent links.
3. `detect.py` reconstructs a residue connectivity graph and classifies cycles.
4. `cli.py` handles local files, RCSB retrieval, batches, and JSON serialization.

Deposited `struct_conn`, `LINK`, and `SSBOND` records have priority. Conservative
coordinate inference fills missing peptide C–N and disulfide S–S bonds. A residue
is peptide-like when its observed atoms contain an amino-acid backbone; residue
names and chain length are not the primary detector.

The scientific dependency chain is:

```text
chemical modification -> conformational response -> molecular recognition -> target selectivity
```

The package boundaries reflect that chain without implying that all stages are
already implemented:

- `core`: shared schemas, identifiers, and atom/link provenance.
- `io`: PDB/mmCIF ingestion and deposited metadata.
- `chemistry`: chemical identity, residue annotation, entity scope, and cycle topology.
- `conformation`: torsions, geometry, alignment, and backbone comparison.
- `comparison`: same-identity comparison and candidate modification mapping.
- `recognition`: interface observation schemas; detailed analysis is deferred.
- `atlas`: namespace reserved for a later public knowledgebase.

One chemical identity may have many experimental instances and conformational
observations. Coordinate models and symmetry-related copies are observations,
not automatically new chemical identities. Every chemical assertion should be
traceable to deposited, CCD, polymer, or coordinate-fallback evidence.

## Current scientific gaps

- expand/select biological assemblies and type target polymers
- reconstruct components when CCD or coordinate atoms are incomplete
- implement IMHB, SASA/polar SASA, ring-shape, and mature interface fingerprints
- curate chemically reviewed matched analogue pairs and target-comparison series
