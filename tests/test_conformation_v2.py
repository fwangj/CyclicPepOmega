import pytest

from cyclicpepomega.conformation.imhb import intramolecular_hydrogen_bonds
from cyclicpepomega.conformation.ring_shape import ring_shape_descriptors
from cyclicpepomega.conformation.sasa import solvent_exposure
from cyclicpepomega.core.schema import Atom, Residue, ResidueKey


def _residue(seq, x, name="ALA"):
    key = ResidueKey(1, "A", str(seq))
    atoms = {
        "N": Atom("N", "N", (x, 0.0, 0.0)),
        "CA": Atom("CA", "C", (x + 0.5, 1.0, 0.0)),
        "C": Atom("C", "C", (x + 1.0, 0.0, 0.0)),
        "O": Atom("O", "O", (x + 1.1, 0.0, 0.0)),
    }
    return key, Residue(key, name, "1", atoms)


def test_imhb_reports_heavy_atom_inference():
    k1, r1 = _residue(1, 0.0)
    k2, r2 = _residue(2, 8.0)
    k3, r3 = _residue(3, 2.0)
    r3.atoms["O"].xyz = (0.2, 0.0, 0.0)
    result = intramolecular_hydrogen_bonds([k1, k2, k3], {k1: r1, k2: r2, k3: r3}, min_sequence_separation=2)
    assert result["count"] > 0
    assert result["method"]["version"] == "cpo-imhb-1"
    assert "heavy_atom" in result["intramolecular_hydrogen_bonds"][0]["geometry_evidence"]


def test_sasa_reports_polar_and_nonpolar_components():
    keys = []
    residues = {}
    for index, x in enumerate((0.0, 4.0, 8.0), start=1):
        key, residue = _residue(index, x)
        keys.append(key)
        residues[key] = residue
    result = solvent_exposure(keys, residues)
    assert result["total_sasa_angstrom2"] > 0
    assert result["polar_sasa_angstrom2"] > 0
    assert result["nonpolar_sasa_angstrom2"] > 0
    assert result["method"]["version"] == "cpo-sasa-shrake-rupley-lite-1"


def test_ring_shape_is_translation_invariant():
    keys = []
    residues = {}
    shifted = {}
    for index, x in enumerate((0.0, 2.0, 4.0), start=1):
        key, residue = _residue(index, x)
        keys.append(key)
        residues[key] = residue
        _, other = _residue(index, x + 100.0)
        shifted[key] = other
    left = ring_shape_descriptors(keys, residues)
    right = ring_shape_descriptors(keys, shifted)
    assert left["principal_moments"] == pytest.approx(right["principal_moments"])
    assert left["planarity_rmsd_angstrom"] == pytest.approx(right["planarity_rmsd_angstrom"])
