from cyclicpepomega.comparison.modification import build_matched_analogue_pair


def _peptide(component, stereo, methyl, identity):
    annotation = {"stereochemistry": stereo, "n_methylated_backbone_n": methyl}
    return {
        "chemical_identity_id": identity,
        "residues": [{"name": component}],
        "residue_annotations": [annotation],
        "molecular_graph": {"nodes": [
            {"residue_index": 1, "atom_id": atom} for atom in ("N", "CA", "C")
        ]},
    }


def test_atom_mapped_stereo_change_is_explicit():
    pair = build_matched_analogue_pair(
        _peptide("PHE", "L", False, "left"),
        _peptide("DPH", "D", False, "right"),
    )
    assert pair.mapping_status == "candidate_multi_edit"
    assert {event.event_type for event in pair.modifications} == {"component_substitution", "stereochemistry_change"}
    assert ("CA", "CA") in pair.modifications[0].mapped_atoms
