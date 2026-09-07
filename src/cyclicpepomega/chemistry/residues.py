from __future__ import annotations

from ..core.schema import Residue, ResidueKey

CANONICAL_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def noncanonical_annotations(
    keys: list[ResidueKey], residues: dict[ResidueKey, Residue]
) -> list[dict[str, str]]:
    """Return name-based v1 annotations; unknown is not silently canonicalized."""
    return [
        {"original_id": key.label(), "component_id": residues[key].name}
        for key in keys
        if residues[key].name.upper() not in CANONICAL_RESIDUES
    ]
