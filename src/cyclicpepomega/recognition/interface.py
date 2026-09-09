from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from math import dist

from ..core.schema import Residue, ResidueKey
from .schema import InterfaceObservation, InteractionObservation

INTERFACE_METHOD_VERSION = "cpo-interface-1"
POLAR = {"N", "O", "S", "P"}
HYDROPHOBIC = {"C", "S"}
AROMATIC_NAMES = {"PHE", "TYR", "TRP", "HIS"}
POSITIVE_NAMES = {"ARG", "LYS", "HIS"}
NEGATIVE_NAMES = {"ASP", "GLU"}


def peptide_target_contacts(
    conformation_id: str,
    peptide_instance_id: str,
    peptide_keys: list[ResidueKey],
    target_keys: list[ResidueKey],
    residues: dict[ResidueKey, Residue],
    *,
    assembly_id: str | None = None,
    cutoff_angstrom: float = 5.0,
) -> dict[str, object]:
    """Measure transparent asymmetric-unit peptide-target contacts."""

    interactions: list[InteractionObservation] = []
    for peptide_key in peptide_keys:
        peptide_residue = residues[peptide_key]
        for target_key in target_keys:
            target_residue = residues[target_key]
            for peptide_atom in peptide_residue.atoms.values():
                for target_atom in target_residue.atoms.values():
                    distance = dist(peptide_atom.xyz, target_atom.xyz)
                    if distance > cutoff_angstrom:
                        continue
                    interaction_type = _classify_contact(
                        peptide_residue.name, peptide_atom.element,
                        target_residue.name, target_atom.element,
                        distance,
                    )
                    interactions.append(InteractionObservation(
                        interaction_type=interaction_type,
                        peptide_residue_id=peptide_key.label(),
                        peptide_atom_id=peptide_atom.name,
                        target_residue_id=target_key.label(),
                        target_atom_id=target_atom.name,
                        distance_angstrom=distance,
                        method_version=INTERFACE_METHOD_VERSION,
                    ))
    target_instance = ",".join(sorted({key.chain for key in target_keys})) or "none"
    digest = sha256(f"{conformation_id}|{peptide_instance_id}|{target_instance}|{assembly_id}".encode()).hexdigest()[:16]
    observation = InterfaceObservation(
        interface_id=f"iface_{digest}",
        conformation_id=conformation_id,
        peptide_instance_id=peptide_instance_id,
        target_instance_id=target_instance,
        assembly_id=assembly_id,
        interactions=interactions,
        fingerprint_version=INTERFACE_METHOD_VERSION,
        warnings=["contacts are asymmetric-unit unless assembly_id/provenance says otherwise"],
    )
    payload = asdict(observation)
    payload["interaction_count"] = len(interactions)
    payload["method"] = {"name": "cyclicpepomega_interface", "version": INTERFACE_METHOD_VERSION, "cutoff_angstrom": cutoff_angstrom}
    return payload


def _classify_contact(peptide_name: str, peptide_element: str, target_name: str, target_element: str, distance: float) -> str:
    pe = peptide_element.upper()
    te = target_element.upper()
    if pe in POLAR and te in POLAR and distance <= 3.5:
        return "hydrogen_bond_candidate"
    if {pe, te} <= HYDROPHOBIC:
        if peptide_name.upper() in AROMATIC_NAMES or target_name.upper() in AROMATIC_NAMES:
            return "aromatic_or_hydrophobic_contact"
        return "hydrophobic_contact"
    if (pe in POLAR or te in POLAR) and {peptide_name.upper(), target_name.upper()} & (POSITIVE_NAMES | NEGATIVE_NAMES):
        return "salt_bridge_or_polar_contact_candidate"
    return "van_der_waals_contact"
