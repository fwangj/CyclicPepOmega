from cyclicpepomega.comparison.discovery import discover_analogue_candidates
from cyclicpepomega.core.schema import Atom, ComponentAtom, ComponentBond, ComponentDefinition, CovalentBond, Residue, ResidueKey, StructureData
from cyclicpepomega.detect import analyze_structure
from cyclicpepomega.recognition.interface import peptide_target_contacts


def _component(name):
    comp = ComponentDefinition(name, "L-peptide linking")
    comp.atoms = {atom: ComponentAtom(atom, element) for atom, element in {"N": "N", "CA": "C", "C": "C", "O": "O"}.items()}
    comp.bonds = [ComponentBond("N", "CA", "sing"), ComponentBond("CA", "C", "sing"), ComponentBond("C", "O", "doub")]
    return comp


def _residue(seq, x, chain="A", name="ALA"):
    key = ResidueKey(1, chain, str(seq))
    residue = Residue(key, name, "1", {
        "N": Atom("N", "N", (x, 0.0, 0.0)),
        "CA": Atom("CA", "C", (x + 0.5, 1.0, 0.0)),
        "C": Atom("C", "C", (x + 1.0, 0.0, 0.0)),
        "O": Atom("O", "O", (x + 1.2, 0.0, 0.0)),
    })
    return key, residue


def _cyclic_report(second_name="GLY", shift=0.0, structure_id="TEST"):
    k1, r1 = _residue(1, shift, name="ALA")
    k2, r2 = _residue(2, shift + 3.0, name=second_name)
    residues = {k1: r1, k2: r2}
    bonds = [
        CovalentBond(k1, "C", k2, "N", "test", "peptide"),
        CovalentBond(k2, "C", k1, "N", "test", "peptide"),
    ]
    data = StructureData(structure_id, residues, bonds, components={"ALA": _component("ALA"), "GLY": _component("GLY"), "SER": _component("SER")})
    return analyze_structure(data)


def test_interface_contacts_are_deterministic_and_typed():
    pk, peptide = _residue(1, 0.0, "A", "ALA")
    tk, target = _residue(1, 2.0, "B", "ASP")
    result = peptide_target_contacts("conf1", "pep1", [pk], [tk], {pk: peptide, tk: target})
    assert result["interaction_count"] > 0
    assert result["interface_id"].startswith("iface_")
    assert result["interactions"][0]["method_version"] == "cpo-interface-1"


def test_discover_analogue_candidates_reports_one_edit_pair():
    left = _cyclic_report("GLY", 0.0, "LEFT")
    right = _cyclic_report("SER", 1.0, "RIGHT")
    result = discover_analogue_candidates([left, right])
    assert result["candidate_count"] == 1
    candidate = result["candidates"][0]
    assert candidate["mapping_status"] == "candidate_single_edit"
    assert candidate["delta_sasa_angstrom2"] is not None
    assert "delta_planarity_rmsd_angstrom" in candidate["ring_shape_deltas"]
