# CyclicPepOmega

CyclicPepOmega is a connectivity-first pipeline for detecting cyclic peptides in
experimental PDB structures. The first milestone accepts a PDB ID or a local
PDB/mmCIF file and emits a JSON report containing every detected cyclic peptide
plus ambiguous/noncyclic peptide-like candidates.

## Current scope

- PDB coordinates plus `LINK` and `SSBOND` records
- mmCIF coordinates plus `_struct_conn` through the optional `gemmi` backend
- residue recognition by backbone atom topology (`N`, `CA`, `C`), not a fixed
  amino-acid name list
- deposited covalent bonds augmented by conservative C–N and S–S geometry
- head-to-tail, head-to-side-chain, side-chain-to-side-chain, and
  disulfide-constrained classification
- preservation of original chain/residue identifiers in JSON

## Install and run

```bash
python -m pip install -e '.[gemmi,test]'
cyclicpepomega analyze 8ABC
cpo analyze structure.cif -o report.json
cpo batch pdb_ids.txt -o reports.json
cpo validate validation/public_set.json -o validation-results.json
```

Output files are never overwritten implicitly.

Development and validation are run on `192.168.10.201` with the dedicated
environment at `/adata/envs/cyclicpepomega`. Large structures and generated
reports belong under `/adata/test/cyclicpepomega`, not in this repository.

## Development test

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Important current limitations

This is the first vertical slice, not yet the full atlas. Biological assembly
expansion, complete chemical-component graph reconstruction, modification/D
stereochemistry/N-methyl annotation, target/interface analysis, descriptors,
database persistence, and conformational comparison remain later milestones.
Geometry-inferred bonds are deliberately conservative and retain linear or
unresolved components in the ambiguity section rather than silently dropping
them.

## Scientific direction

CyclicPepOmega is designed around **chemical modification → conformation →
molecular recognition → selectivity**. The current release establishes public-PDB
parsing and connectivity-based candidate detection; it does not yet calculate
conformational energies, predict permeability, infer binding affinity, or claim a
complete PDB census.
