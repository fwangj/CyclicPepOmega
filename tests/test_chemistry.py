from cyclicpepomega.chemistry.identity import n_methyl_atoms, residue_stereochemistry
from cyclicpepomega.core.schema import ComponentAtom, ComponentBond, ComponentDefinition


def test_ccd_stereo_and_n_methyl_annotation():
    component = ComponentDefinition("TEST", "D-peptide linking")
    component.atoms = {
        "N": ComponentAtom("N", "N"), "CA": ComponentAtom("CA", "C", "R"),
        "CM": ComponentAtom("CM", "C"),
        "H1": ComponentAtom("H1", "H"), "H2": ComponentAtom("H2", "H"), "H3": ComponentAtom("H3", "H"),
    }
    component.bonds = [
        ComponentBond("N", "CA", "sing"), ComponentBond("N", "CM", "sing"),
        ComponentBond("CM", "H1", "sing"), ComponentBond("CM", "H2", "sing"), ComponentBond("CM", "H3", "sing"),
    ]
    assert residue_stereochemistry(component)[0] == "D"
    assert n_methyl_atoms(component) == ["CM"]

