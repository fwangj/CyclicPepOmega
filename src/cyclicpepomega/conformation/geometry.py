from __future__ import annotations

import hashlib
import math
from typing import Iterable

import numpy as np

from ..core.schema import CovalentBond, Residue, ResidueKey
from .imhb import intramolecular_hydrogen_bonds
from .ring_shape import ring_shape_descriptors
from .sasa import solvent_exposure

BACKBONE_ATOMS = ("N", "CA", "C")


def _dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    b0 = -(b - a)
    b1 = c - b
    b2 = d - c
    norm = np.linalg.norm(b1)
    if norm == 0:
        return float("nan")
    b1 /= norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) == 0 or np.linalg.norm(w) == 0:
        return float("nan")
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _coord(residue: Residue, atom_name: str) -> np.ndarray | None:
    atom = residue.atom(atom_name)
    return np.asarray(atom.xyz, dtype=float) if atom else None


def is_backbone_cyclic(keys: list[ResidueKey], bonds: Iterable[CovalentBond]) -> bool:
    if not keys:
        return False
    first, last = keys[0], keys[-1]
    return any(
        {bond.left, bond.right} == {first, last}
        and {bond.left_atom.upper(), bond.right_atom.upper()} == {"C", "N"}
        for bond in bonds
    )


def backbone_torsions(
    keys: list[ResidueKey], residues: dict[ResidueKey, Residue], cyclic: bool
) -> list[dict[str, float | int | None]]:
    result = []
    n = len(keys)
    for i, key in enumerate(keys):
        prev_key = keys[(i - 1) % n] if cyclic or i > 0 else None
        next_key = keys[(i + 1) % n] if cyclic or i + 1 < n else None
        current = residues[key]
        phi = psi = omega = None
        if prev_key is not None:
            atoms = (_coord(residues[prev_key], "C"), _coord(current, "N"), _coord(current, "CA"), _coord(current, "C"))
            if all(atom is not None for atom in atoms):
                phi = _dihedral(*atoms)
        if next_key is not None:
            following = residues[next_key]
            atoms = (_coord(current, "N"), _coord(current, "CA"), _coord(current, "C"), _coord(following, "N"))
            if all(atom is not None for atom in atoms):
                psi = _dihedral(*atoms)
            atoms = (_coord(current, "CA"), _coord(current, "C"), _coord(following, "N"), _coord(following, "CA"))
            if all(atom is not None for atom in atoms):
                omega = _dihedral(*atoms)
        result.append({"residue_index": i + 1, "phi_degrees": phi, "psi_degrees": psi, "omega_degrees": omega})
    return result


def radius_of_gyration(keys: list[ResidueKey], residues: dict[ResidueKey, Residue]) -> float | None:
    coords = [atom.xyz for key in keys for atom in residues[key].atoms.values() if atom.element.upper() != "H"]
    if not coords:
        return None
    array = np.asarray(coords, dtype=float)
    center = array.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((array - center) ** 2, axis=1))))


def backbone_coordinates(keys: list[ResidueKey], residues: dict[ResidueKey, Residue]) -> list[dict[str, object]]:
    return [
        {"residue_index": i, "atom_name": atom_name, "xyz_angstrom": list(residues[key].atoms[atom_name].xyz)}
        for i, key in enumerate(keys, 1)
        for atom_name in BACKBONE_ATOMS
        if atom_name in residues[key].atoms
    ]


def describe_conformation(
    structure_id: str, identity_id: str, keys: list[ResidueKey],
    residues: dict[ResidueKey, Residue], bonds: list[CovalentBond],
    has_macromolecular_contact: bool = False,
) -> dict[str, object]:
    cyclic = is_backbone_cyclic(keys, bonds)
    key = f"{structure_id}|{keys[0].model}|{identity_id}|experimental"
    imhb = intramolecular_hydrogen_bonds(keys, residues)
    sasa = solvent_exposure(keys, residues)
    ring_shape = ring_shape_descriptors(keys, residues)
    return {
        "conformation_id": f"conf_{hashlib.sha256(key.encode()).hexdigest()[:24]}",
        "kind": "experimental_bound_asymmetric_unit_contact" if has_macromolecular_contact else "experimental_context_unresolved",
        "coordinate_model": keys[0].model,
        "backbone_cyclic": cyclic,
        "torsions": backbone_torsions(keys, residues, cyclic),
        "radius_of_gyration_angstrom": radius_of_gyration(keys, residues),
        "backbone_coordinates": backbone_coordinates(keys, residues),
        "imhb": imhb,
        "sasa": sasa,
        "ring_shape": ring_shape,
        "method": {"name": "cyclicpepomega_geometry", "version": "2"},
        "warnings": [
            "geometric (not mass-weighted) radius of gyration",
            *imhb.get("warnings", []),
            *sasa.get("warnings", []),
            *ring_shape.get("warnings", []),
        ]
    }


def _kabsch_rmsd(moving: np.ndarray, reference: np.ndarray) -> float:
    moving_centered = moving - moving.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)
    covariance = moving_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    fitted = moving_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((fitted - reference_centered) ** 2, axis=1))))


def compare_backbones(
    left: dict[str, object], right: dict[str, object], *, allow_component_substitutions: bool = False
) -> dict[str, object]:
    """Compare same-chemistry conformations across valid cyclic registrations."""
    left_res = left["residues"]
    right_res = right["residues"]
    if len(left_res) != len(right_res):
        raise ValueError("backbone comparison requires equal normalized residue counts")
    n = len(left_res)
    candidates = []
    left_names = [r["name"] for r in left_res]
    right_names = [r["name"] for r in right_res]
    left_coords = {(p["residue_index"] - 1, p["atom_name"]): p["xyz_angstrom"] for p in left["conformation"]["backbone_coordinates"]}
    right_coords = {(p["residue_index"] - 1, p["atom_name"]): p["xyz_angstrom"] for p in right["conformation"]["backbone_coordinates"]}
    for shift in range(n):
        mapping = [(i, (i + shift) % n) for i in range(n)]
        substitutions = sum(left_names[i] != right_names[j] for i, j in mapping)
        if substitutions and not allow_component_substitutions:
            continue
        atom_pairs = [
            (left_coords[(i, atom)], right_coords[(j, atom)])
            for i, j in mapping for atom in BACKBONE_ATOMS
            if (i, atom) in left_coords and (j, atom) in right_coords
        ]
        if len(atom_pairs) < 3:
            continue
        a, b = (np.asarray([pair[k] for pair in atom_pairs], dtype=float) for k in (0, 1))
        candidates.append((substitutions, _kabsch_rmsd(a, b), shift, len(atom_pairs), mapping))
    if not candidates:
        raise ValueError("no chemically compatible cyclic registration with at least three mapped backbone atoms")
    substitutions, rmsd, shift, atom_count, mapping = min(candidates)
    left_t = left["conformation"]["torsions"]
    right_t = right["conformation"]["torsions"]
    deltas = []
    for i, j in mapping:
        for name in ("phi_degrees", "psi_degrees", "omega_degrees"):
            x, y = left_t[i][name], right_t[j][name]
            if x is not None and y is not None:
                deltas.append(abs(((x - y + 180.0) % 360.0) - 180.0))
    torsion_rms = math.sqrt(sum(delta * delta for delta in deltas) / len(deltas)) if deltas else None
    left_conf = left["conformation"]
    right_conf = right["conformation"]
    left_sasa = left_conf.get("sasa", {})
    right_sasa = right_conf.get("sasa", {})
    left_ring = left_conf.get("ring_shape", {})
    right_ring = right_conf.get("ring_shape", {})
    left_imhb = set(str(left_conf.get("imhb", {}).get("fingerprint", "")).split(";")) - {""}
    right_imhb = set(str(right_conf.get("imhb", {}).get("fingerprint", "")).split(";")) - {""}
    return {
        "left_conformation_id": left_conf["conformation_id"],
        "right_conformation_id": right_conf["conformation_id"],
        "backbone_rmsd_angstrom": rmsd,
        "mapped_backbone_atom_count": atom_count,
        "cyclic_permutation_shift": shift,
        "component_substitution_count": substitutions,
        "torsion_rms_difference_degrees": torsion_rms,
        "mapped_torsion_count": len(deltas),
        "delta_radius_of_gyration_angstrom": _delta(left_conf.get("radius_of_gyration_angstrom"), right_conf.get("radius_of_gyration_angstrom")),
        "delta_imhb_count": _delta(left_conf.get("imhb", {}).get("count"), right_conf.get("imhb", {}).get("count")),
        "gained_imhb_fingerprint_tokens": sorted(right_imhb - left_imhb),
        "lost_imhb_fingerprint_tokens": sorted(left_imhb - right_imhb),
        "delta_total_sasa_angstrom2": _delta(left_sasa.get("total_sasa_angstrom2"), right_sasa.get("total_sasa_angstrom2")),
        "delta_polar_sasa_angstrom2": _delta(left_sasa.get("polar_sasa_angstrom2"), right_sasa.get("polar_sasa_angstrom2")),
        "delta_nonpolar_sasa_angstrom2": _delta(left_sasa.get("nonpolar_sasa_angstrom2"), right_sasa.get("nonpolar_sasa_angstrom2")),
        "delta_ring_planarity_rmsd_angstrom": _delta(left_ring.get("planarity_rmsd_angstrom"), right_ring.get("planarity_rmsd_angstrom")),
        "delta_ring_asphericity": _delta(left_ring.get("asphericity"), right_ring.get("asphericity")),
        "method": {"name": "Kabsch_same_component_cyclic_permutations", "version": "2"},
        "warnings": [
            "reverse sequence and side-chain symmetry mappings are not evaluated",
            *( ["backbone comparison includes mapped positions with component substitutions"] if substitutions else [] ),
        ],
    }


def _delta(left: object, right: object) -> float | None:
    if left is None or right is None:
        return None
    return float(right) - float(left)
