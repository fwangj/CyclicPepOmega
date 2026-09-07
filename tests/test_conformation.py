import pytest

from cyclicpepomega.conformation.geometry import compare_backbones


def _instance(coords, names=("ALA", "GLY")):
    points = []
    for residue_index, residue_points in enumerate(coords, 1):
        for atom_name, xyz in zip(("N", "CA", "C"), residue_points):
            points.append({"residue_index": residue_index, "atom_name": atom_name, "xyz_angstrom": xyz})
    return {
        "residues": [{"name": name} for name in names],
        "conformation": {
            "conformation_id": "test",
            "backbone_coordinates": points,
            "torsions": [{"phi_degrees": 10.0, "psi_degrees": 20.0, "omega_degrees": 180.0} for _ in names],
        },
    }


def test_kabsch_removes_rigid_translation():
    left = _instance([[(0, 0, 0), (1, 0, 0), (1, 1, 0)], [(2, 1, 0), (2, 2, 0), (3, 2, 0)]])
    right = _instance([[tuple(v + 5 for v in p) for p in residue] for residue in [
        [(0, 0, 0), (1, 0, 0), (1, 1, 0)], [(2, 1, 0), (2, 2, 0), (3, 2, 0)]
    ]])
    result = compare_backbones(left, right)
    assert result["backbone_rmsd_angstrom"] == pytest.approx(0.0, abs=1e-10)
    assert result["mapped_backbone_atom_count"] == 6


def test_cyclic_permutation_is_evaluated():
    first = [(0, 0, 0), (1, 0, 0), (1, 1, 0)]
    second = [(2, 1, 0), (2, 2, 0), (3, 2, 0)]
    left = _instance([first, second], names=("ALA", "GLY"))
    right = _instance([second, first], names=("GLY", "ALA"))
    result = compare_backbones(left, right)
    assert result["backbone_rmsd_angstrom"] == pytest.approx(0.0, abs=1e-10)
    assert result["cyclic_permutation_shift"] == 1
