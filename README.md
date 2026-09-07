# CyclicPepOmega

CyclicPepOmega is a public-data-driven computational platform for cyclic-peptide
chemical normalization, conformational analysis, structural comparison, and
future molecular-recognition and selectivity studies. It accepts a PDB ID or a
local PDB/mmCIF file and emits a machine-readable report of detected cyclic
peptides plus ambiguous or noncyclic peptide-like candidates.

## Implemented in v0.2.0

- PDB coordinates plus `LINK` and `SSBOND` records
- mmCIF coordinates plus `_struct_conn` through the optional `gemmi` backend
- residue recognition by backbone atom topology (`N`, `CA`, `C`), not a fixed
  amino-acid name list
- deposited covalent bonds, CCD component graphs, polymer links, and conservative
  coordinate fallbacks with explicit provenance
- atom-level, versioned chemical identities that collapse equivalent coordinate
  observations
- head-to-tail, head-to-side-chain, lactam, disulfide, thioether, and
  multiple-crosslink classification with unresolved cases retained
- atom-level D/L and N-methylation annotation and noncanonical-residue handling
- backbone phi/psi/omega, radius of gyration, backbone representation, and
  atom-mapped RMSD with cyclic-permutation search
- same-identity comparison and candidate atom-mapped analogue differences
- preservation of original chain/residue identifiers in JSON

## Install and run

```bash
python -m pip install -e '.[gemmi,test]'
cyclicpepomega analyze 8ABC
cpo analyze structure.cif -o report.json
cpo batch pdb_ids.txt -o reports.json
cpo validate validation/public_set.json -o validation-results.json
cpo validate validation/public_set_v2.json --structure-dir /adata/test/cyclicpepomega/mmcif
cpo compare 1CWA 2RMC
cpo map-modification 1CWA 1CWB
```

Output files are never overwritten implicitly.

Downloaded structures and generated reports are deliberately excluded from this
repository. The versioned public validation manifests contain only PDB IDs and
expected annotations.

## Development test

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Important current limitations

This is an early research-software baseline, not a complete atlas. Biological
assembly expansion, exhaustive CCD reconstruction, IMHB and SASA analysis,
ring-shape descriptors, production interface fingerprints, curated matched-pair
mining, databases, APIs, and the website are planned rather than implemented.
Pure disulfide topology cannot by itself distinguish a cyclic peptide from a
disulfide-rich protein; such cases require versioned entity-scope evidence or
remain unresolved. Coordinate-inferred bonds never silently override deposited
chemistry.

## Scientific direction

CyclicPepOmega is designed around **chemical modification → conformation →
molecular recognition → selectivity**. The current release establishes public-PDB
parsing and connectivity-based candidate detection; it does not yet calculate
conformational energies, predict permeability, infer binding affinity, or claim a
complete PDB census.

See [the architecture](docs/architecture.md), [scientific scope](docs/scientific_scope.md),
[data model](docs/data_model.md), [roadmap](docs/roadmap.md), and
[v0.2.0 baseline manifest](docs/BASELINE_v0.2.0.md).
