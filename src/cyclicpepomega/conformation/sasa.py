from __future__ import annotations

from math import pi, sqrt

from ..core.schema import Residue, ResidueKey

SASA_METHOD_VERSION = "cpo-sasa-shrake-rupley-lite-1"
POLAR_ELEMENTS = {"N", "O", "S", "P"}
NONPOLAR_ELEMENTS = {"C", "F", "CL", "BR", "I"}
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}


def solvent_exposure(
    keys: list[ResidueKey],
    residues: dict[ResidueKey, Residue],
    *,
    probe_radius: float = 1.4,
    exposure_threshold_fraction: float = 0.20,
) -> dict[str, object]:
    """Approximate SASA with a deterministic neighbor-occlusion proxy."""

    atoms = []
    for key in keys:
        for atom in residues[key].atoms.values():
            element = atom.element.upper()
            if element == "H":
                continue
            atoms.append((key, atom.name, element, atom.xyz))
    atom_rows = []
    unknown = []
    for key, atom_name, element, xyz in atoms:
        radius = VDW_RADII.get(element, 1.70)
        sphere = 4.0 * pi * (radius + probe_radius) ** 2
        occlusion = 0.0
        for other_key, other_name, other_element, other_xyz in atoms:
            if other_key == key and other_name == atom_name:
                continue
            cutoff = radius + VDW_RADII.get(other_element, 1.70) + 2.0 * probe_radius
            distance = sqrt(sum((xyz[i] - other_xyz[i]) ** 2 for i in range(3)))
            if distance < cutoff:
                occlusion += max(0.0, (cutoff - distance) / cutoff) * 0.08
        exposed_fraction = max(0.0, min(1.0, 1.0 - occlusion))
        sasa = sphere * exposed_fraction
        polarity = "polar" if element in POLAR_ELEMENTS else "nonpolar" if element in NONPOLAR_ELEMENTS else "unknown"
        if polarity == "unknown":
            unknown.append(f"{key.label()}:{atom_name}:{element}")
        atom_rows.append({
            "residue_id": key.label(),
            "atom_id": atom_name,
            "element": element,
            "polarity": polarity,
            "sasa_angstrom2": sasa,
            "exposed_fraction": exposed_fraction,
            "exposed": exposed_fraction >= exposure_threshold_fraction,
        })
    total = sum(row["sasa_angstrom2"] for row in atom_rows)
    polar = sum(row["sasa_angstrom2"] for row in atom_rows if row["polarity"] == "polar")
    nonpolar = sum(row["sasa_angstrom2"] for row in atom_rows if row["polarity"] == "nonpolar")
    exposed_hbd = sum(1 for row in atom_rows if row["element"] in {"N", "O", "S"} and row["exposed"])
    exposed_hba = exposed_hbd
    buried_polar = sum(1 for row in atom_rows if row["polarity"] == "polar" and not row["exposed"])
    return {
        "method": {
            "name": "cyclicpepomega_sasa_proxy",
            "version": SASA_METHOD_VERSION,
            "probe_radius_angstrom": probe_radius,
            "policy": "deterministic neighbor-occlusion proxy; not a replacement for FreeSASA/NACCESS",
            "exposure_threshold_fraction": exposure_threshold_fraction,
        },
        "total_sasa_angstrom2": total,
        "polar_sasa_angstrom2": polar,
        "nonpolar_sasa_angstrom2": nonpolar,
        "exposed_hbd_count": exposed_hbd,
        "exposed_hba_count": exposed_hba,
        "buried_polar_atom_count": buried_polar,
        "atom_exposure": atom_rows,
        "unsupported_atom_classifications": unknown,
        "warnings": ["approximate SASA descriptor for comparative screening, not permeability prediction"],
    }
