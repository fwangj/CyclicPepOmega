# Milestone 1 architecture

The initial pipeline has four deliberately separate layers:

1. `io/structures.py` converts PDB or mmCIF records into a format-neutral structure.
2. `core/schema.py` preserves atoms, author residue IDs, entity IDs, and covalent links.
3. `detect.py` reconstructs a residue connectivity graph and classifies cycles.
4. `cli.py` handles local files, RCSB retrieval, batches, and JSON serialization.

Deposited `struct_conn`, `LINK`, and `SSBOND` records have priority. Conservative
coordinate inference fills missing peptide C–N and disulfide S–S bonds. A residue
is peptide-like when its observed atoms contain an amino-acid backbone; residue
names and chain length are not the primary detector.

## Current scientific gaps

- consume `chem_comp_bond` graphs for non-standard atom naming and ligand entities
- classify D configuration and N-methylation chemically rather than by component ID
- expand/select biological assemblies and type target polymers
- assign explicit confidence and ambiguity evidence to every inferred bond
- curate and regress a multi-architecture positive/negative validation panel
