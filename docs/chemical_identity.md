# Chemical identity and evidence (v0.2)

`cpo-chemical-identity-1` separates chemical identity from PDB occurrence and
coordinate model. Its SHA-256-derived ID is constructed from ordered component
graphs, CCD stereochemistry/N-methyl features, inter-residue bonds, cyclization,
crosslinks, and covalent connector components. PDB accession, model number,
chain/auth numbering, coordinates, and evidence-source labels are excluded from
the chemical hash.

For cyclic backbones, all sequence rotations and both traversal directions are
canonicalized before hashing. Traversal reversal is used for identity
canonicalization of an undirected atom graph; it is not interpreted as reversed
peptide chemistry during conformational alignment.

Evidence sources remain explicit:

- `embedded_pdbx_chem_comp`: deposited component atom/bond/stereo definition;
- `deposited_struct_conn`: deposited inter-component covalent connection;
- `polymer_connectivity`: adjacency from `entity_poly_seq` within one label asym;
- `coordinate_inferred`: conservative C–N or S–S fallback.

Deposited evidence is never overwritten by coordinate inference. Missing CCD
definitions, ambiguous L/D types, missing atoms and inferred closure are warnings,
not silently canonicalized facts.

The current identity is suitable for grouping the same deposited chemistry across
models and tested PDB entries. It is not yet a universal tautomer/protonation-aware
small-molecule canonicalization and does not resolve chemically equivalent CCD
aliases.

