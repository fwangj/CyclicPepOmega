from __future__ import annotations

from .modification import build_matched_analogue_pair
from ..conformation.geometry import compare_backbones

DISCOVERY_METHOD_VERSION = "cpo-analogue-discovery-1"


def discover_analogue_candidates(
    reports: list[dict[str, object]],
    *,
    max_component_edits: int = 1,
) -> dict[str, object]:
    """Discover bounded candidate analogue pairs from analyzed public structures."""

    observations = []
    for report in reports:
        source_id = report.get("structure", {}).get("id", "unknown") if isinstance(report.get("structure"), dict) else "unknown"
        method = report.get("structure", {}).get("experimental_method") if isinstance(report.get("structure"), dict) else None
        for item in report.get("cyclic_peptides", []):
            observations.append({"source_id": source_id, "experimental_method": method, "peptide": item})
    candidates = []
    failures: dict[str, int] = {}
    for left_index, left in enumerate(observations):
        for right in observations[left_index + 1:]:
            left_peptide = left["peptide"]
            right_peptide = right["peptide"]
            if left_peptide["instance_id"] == right_peptide["instance_id"]:
                continue
            if len(left_peptide["residues"]) != len(right_peptide["residues"]):
                failures["unequal_monomer_count"] = failures.get("unequal_monomer_count", 0) + 1
                continue
            try:
                mapping = build_matched_analogue_pair(left_peptide, right_peptide)
            except ValueError as exc:
                failures[str(exc)] = failures.get(str(exc), 0) + 1
                continue
            edit_count = len(mapping.modifications)
            if edit_count == 0 or edit_count > max_component_edits:
                failures["not_one_edit_candidate"] = failures.get("not_one_edit_candidate", 0) + 1
                continue
            try:
                conformation = compare_backbones(left_peptide, right_peptide, allow_component_substitutions=True)
            except ValueError as exc:
                failures[str(exc)] = failures.get(str(exc), 0) + 1
                continue
            candidates.append({
                "left_peptide_id": left_peptide["chemical_identity_id"],
                "right_peptide_id": right_peptide["chemical_identity_id"],
                "left_observation_id": left_peptide["instance_id"],
                "right_observation_id": right_peptide["instance_id"],
                "pdb_observations": [left["source_id"], right["source_id"]],
                "experimental_method_context": [left["experimental_method"], right["experimental_method"]],
                "chemical_mapping": mapping.to_dict(),
                "modification_events": [event.__dict__ for event in mapping.modifications],
                "mapping_status": mapping.mapping_status,
                "backbone_rmsd_angstrom": conformation.get("backbone_rmsd_angstrom"),
                "torsion_rms_difference_degrees": conformation.get("torsion_rms_difference_degrees"),
                "delta_rg_angstrom": conformation.get("delta_radius_of_gyration_angstrom"),
                "delta_imhb_count": conformation.get("delta_imhb_count"),
                "delta_sasa_angstrom2": conformation.get("delta_total_sasa_angstrom2"),
                "delta_polar_sasa_angstrom2": conformation.get("delta_polar_sasa_angstrom2"),
                "ring_shape_deltas": {
                    "delta_planarity_rmsd_angstrom": conformation.get("delta_ring_planarity_rmsd_angstrom"),
                    "delta_asphericity": conformation.get("delta_ring_asphericity"),
                },
                "warnings": [
                    "machine-generated candidate; not manually accepted as a matched analogue pair",
                    "differences are associated with deposited public structures and do not establish causality",
                    *mapping.warnings,
                    *conformation.get("warnings", []),
                ],
            })
    return {
        "schema_version": "0.3-dev",
        "method": {"name": "cyclicpepomega_analogue_discovery", "version": DISCOVERY_METHOD_VERSION, "max_component_edits": max_component_edits},
        "input_observation_count": len(observations),
        "candidate_count": len(candidates),
        "candidates": candidates,
        "failure_categories": failures,
    }
