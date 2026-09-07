from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class InteractionObservation:
    interaction_type: str
    peptide_residue_id: str
    peptide_atom_id: str
    target_residue_id: str
    target_atom_id: str
    distance_angstrom: float | None = None
    angle_degrees: float | None = None
    method_version: str = "uncomputed"


@dataclass
class InterfaceObservation:
    interface_id: str
    conformation_id: str
    peptide_instance_id: str
    target_instance_id: str
    assembly_id: str | None
    interactions: list[InteractionObservation] = field(default_factory=list)
    interface_sasa_angstrom2: float | None = None
    buried_surface_area_angstrom2: float | None = None
    fingerprint_version: str | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class RecognitionComparison:
    left_conformation_id: str
    right_conformation_id: str
    left_interface_id: str
    right_interface_id: str
    gained_interactions: list[str] = field(default_factory=list)
    lost_interactions: list[str] = field(default_factory=list)
    conserved_interactions: list[str] = field(default_factory=list)
