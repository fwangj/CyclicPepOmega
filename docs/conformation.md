# Experimental conformation module (v0.2)

Implemented descriptors are backbone phi/psi/omega (degrees), geometric heavy-atom
radius of gyration (Å), an ordered N/CA/C coordinate representation, and Kabsch
backbone RMSD. Cyclic comparisons test every sequence rotation compatible with the
component mapping. Missing backbone atoms are skipped and the mapped atom/torsion
counts are reported.

`cyclicpepomega compare` requires the same versioned chemical identity. With one
NMR structure it compares the first two models of a shared identity; with two
sources it compares the first matching instance in each. `map-modification`
proposes a residue/CCD-atom mapping for equal-length cyclic analogues and reports
the modification plus backbone/torsion difference. That mapping is explicitly a
candidate requiring chemical curation.

Not implemented: reverse-sequence equivalence in alignment, side-chain symmetry,
unequal-length or bond-edit analogue mapping, side-chain chi, ring-shape metrics,
IMHB, SASA, polar exposure, biological-assembly expansion, energy or population.
An asymmetric-unit macromolecular contact is labeled as such; absence of such a
contact is not called experimental apo.



## Conformation V2 Descriptors

Development builds after v0.2.0 attach additional transparent descriptors to each detected cyclic-peptide conformation:

- `imhb`: intramolecular hydrogen-bond candidates inferred from observed heavy-atom distances. Hydrogens and protonation states are not added silently.
- `sasa`: deterministic neighbor-occlusion SASA proxy with total, polar, nonpolar, exposed donor/acceptor, and buried polar-atom counts. This is not a permeability model.
- `ring_shape`: rigid-transform-invariant backbone descriptors including principal moments, asphericity, planarity RMSD, compactness, and distance-matrix fingerprint.

Backbone comparisons now report deltas for Rg, IMHB count/fingerprint tokens, SASA, polar SASA, nonpolar SASA, ring planarity, and ring asphericity. Missing values remain `null`; unavailable descriptors are not converted to zero.
