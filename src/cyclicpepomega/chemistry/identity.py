from __future__ import annotations

import hashlib
import json
from ..core.schema import ComponentDefinition, CovalentBond, Residue, ResidueKey, StructureData
from .residues import CANONICAL_RESIDUES

IDENTITY_SCHEMA_VERSION = "cpo-chemical-identity-1"


def _component_signature(component: ComponentDefinition | None) -> dict[str, object]:
    if component is None:
        return {"definition": "missing"}
    return {
        "type": component.component_type,
        "atoms": sorted(
            (a.atom_id, a.element, a.stereo, a.aromatic)
            for a in component.atoms.values() if a.element.upper() != "H"
        ),
        "bonds": sorted(
            (min(b.atom_id_1, b.atom_id_2), max(b.atom_id_1, b.atom_id_2), b.order, b.aromatic, b.stereo)
            for b in component.bonds
            if component.atoms.get(b.atom_id_1, None) is None
            or component.atoms.get(b.atom_id_2, None) is None
            or (component.atoms[b.atom_id_1].element.upper() != "H" and component.atoms[b.atom_id_2].element.upper() != "H")
        ),
    }


def residue_stereochemistry(component: ComponentDefinition | None) -> tuple[str, str]:
    """Return L/D/unknown and its chemical evidence, without name heuristics."""
    component_type = (component.component_type or "").lower() if component else ""
    if "d-peptide" in component_type:
        return "D", "CCD chem_comp.type"
    if "l-peptide" in component_type:
        return "L", "CCD chem_comp.type"
    return "unknown", "CCD type does not specify L/D"


def n_methyl_atoms(component: ComponentDefinition | None) -> list[str]:
    """Find methyl carbons directly bonded to backbone N in the CCD graph."""
    if component is None or "N" not in component.atoms:
        return []
    neighbors: dict[str, set[str]] = {}
    for bond in component.bonds:
        neighbors.setdefault(bond.atom_id_1, set()).add(bond.atom_id_2)
        neighbors.setdefault(bond.atom_id_2, set()).add(bond.atom_id_1)
    result = []
    for atom_id in neighbors.get("N", set()) - {"CA", "C"}:
        atom = component.atoms.get(atom_id)
        if not atom or atom.element.upper() != "C":
            continue
        hydrogen_count = sum(
            component.atoms.get(n) is not None and component.atoms[n].element.upper() == "H"
            for n in neighbors.get(atom_id, set())
        )
        heavy_neighbors = sum(
            component.atoms.get(n) is None or component.atoms[n].element.upper() != "H"
            for n in neighbors.get(atom_id, set())
        )
        if hydrogen_count >= 2 and heavy_neighbors == 1:
            result.append(atom_id)
    return sorted(result)


def annotate_residue(index: int, key: ResidueKey, residue: Residue, component: ComponentDefinition | None) -> dict[str, object]:
    stereo, stereo_evidence = residue_stereochemistry(component)
    methyl_atoms = n_methyl_atoms(component)
    return {
        "internal_index": index,
        "original_id": key.label(),
        "label_asym_id": residue.label_asym_id,
        "label_seq_id": residue.label_seq_id,
        "auth_asym_id": residue.auth_asym_id or key.chain,
        "auth_seq_id": residue.auth_seq_id or key.seq,
        "component_id": residue.name,
        "component_type": component.component_type if component else None,
        "canonical_amino_acid": residue.name.upper() in CANONICAL_RESIDUES,
        "stereochemistry": stereo,
        "stereochemistry_evidence": stereo_evidence,
        "n_methylated_backbone_n": bool(methyl_atoms),
        "n_methyl_carbon_atoms": methyl_atoms,
        "annotation_source": component.source if component else "missing_CCD_definition",
    }


def modification_annotations(annotations: list[dict[str, object]]) -> list[dict[str, object]]:
    result = []
    for item in annotations:
        common = {
            "residue_index": item["internal_index"],
            "component_id": item["component_id"],
            "evidence_source": item["annotation_source"],
        }
        if not item["canonical_amino_acid"]:
            result.append({**common, "type": "noncanonical_amino_acid", "atoms": []})
        if item["stereochemistry"] == "D":
            result.append({**common, "type": "d_amino_acid", "atoms": ["CA"], "evidence": item["stereochemistry_evidence"]})
        if item["n_methylated_backbone_n"]:
            result.append({
                **common, "type": "backbone_n_methylation",
                "atoms": ["N", *item["n_methyl_carbon_atoms"]],
                "evidence": "CCD atom/bond graph",
            })
    return result


def atom_level_graph(
    data: StructureData, ordered: list[ResidueKey], bonds: list[CovalentBond],
    auxiliary: list[ResidueKey] | None = None,
) -> dict[str, object]:
    """Serializable CCD + inter-residue molecular graph with explicit evidence."""
    nodes = []
    graph_bonds = []
    auxiliary = auxiliary or []
    graph_keys = [(str(i), key) for i, key in enumerate(ordered, 1)] + [
        (f"aux{i}", key) for i, key in enumerate(auxiliary, 1)
    ]
    for residue_index, key in graph_keys:
        residue = data.residues[key]
        definition = data.components.get(residue.name)
        atom_ids = set(residue.atoms)
        if definition:
            atom_ids.update(definition.atoms)
        for atom_id in sorted(atom_ids):
            ccd_atom = definition.atoms.get(atom_id) if definition else None
            observed_atom = residue.atoms.get(atom_id)
            nodes.append({
                "node_id": f"{residue_index}:{atom_id}",
                "residue_index": int(residue_index) if residue_index.isdigit() else None,
                "auxiliary_component_index": residue_index if not residue_index.isdigit() else None,
                "component_id": residue.name,
                "atom_id": atom_id,
                "element": ccd_atom.element if ccd_atom else observed_atom.element if observed_atom else None,
                "stereo": ccd_atom.stereo if ccd_atom else None,
                "observed": observed_atom is not None,
                "source": "embedded_pdbx_chem_comp" if ccd_atom else "atom_site_only",
            })
        if definition:
            for bond in definition.bonds:
                graph_bonds.append({
                    "atom_1": f"{residue_index}:{bond.atom_id_1}",
                    "atom_2": f"{residue_index}:{bond.atom_id_2}",
                    "order": bond.order,
                    "kind": "intra_component",
                    "evidence_source": definition.source,
                    "evidence_status": "deposited_dictionary",
                })
    indices = {key: label for label, key in graph_keys}
    for bond in bonds:
        if bond.left not in indices or bond.right not in indices:
            continue
        graph_bonds.append({
            "atom_1": f"{indices[bond.left]}:{bond.left_atom}",
            "atom_2": f"{indices[bond.right]}:{bond.right_atom}",
            "order": "single",
            "kind": bond.kind,
            "evidence_source": bond.source,
            "evidence_status": bond.evidence_status,
            "provenance_id": bond.provenance_id,
        })
    return {
        "schema_version": "cpo-atom-graph-1",
        "nodes": nodes,
        "bonds": graph_bonds,
    }


def chemical_identity(
    data: StructureData,
    ordered: list[ResidueKey],
    bonds: list[CovalentBond],
    annotations: list[dict[str, object]],
    cyclic: bool,
    auxiliary: list[ResidueKey] | None = None,
) -> tuple[str, dict[str, object]]:
    """Build a model/PDB-independent, cyclic-registration-invariant identity."""
    n = len(ordered)
    auxiliary = auxiliary or []
    all_keys = ordered + auxiliary
    original_index = {key: i for i, key in enumerate(all_keys)}
    residue_payload = [
        {
            "component_id": data.residues[key].name,
            "component_graph": _component_signature(data.components.get(data.residues[key].name)),
            "stereochemistry": annotations[i]["stereochemistry"],
            "n_methylated": annotations[i]["n_methylated_backbone_n"],
        }
        for i, key in enumerate(ordered)
    ]
    auxiliary_payload = [
        {
            "component_id": data.residues[key].name,
            "component_graph": _component_signature(data.components.get(data.residues[key].name)),
        }
        for key in auxiliary
    ]
    crosslinks = []
    for bond in bonds:
        if bond.left not in original_index or bond.right not in original_index:
            continue
        atom_pair = {bond.left_atom.upper(), bond.right_atom.upper()}
        normalized_kind = (
            "peptide" if atom_pair == {"C", "N"}
            else "disulfide" if atom_pair == {"SG"}
            else "covalent"
        )
        crosslinks.append((original_index[bond.left], bond.left_atom, original_index[bond.right], bond.right_atom, normalized_kind))

    orders = []
    base = list(range(n))
    orientations = (base, list(reversed(base)))
    for orientation in orientations:
        rotations = range(n) if cyclic else range(1)
        for shift in rotations:
            orders.append(orientation[shift:] + orientation[:shift])
    candidates = []
    for order in orders:
        remap = {old: new for new, old in enumerate(order)}
        remap.update({n + i: n + i for i in range(len(auxiliary))})
        payload = {
            "version": IDENTITY_SCHEMA_VERSION,
            "residues": [residue_payload[i] for i in order],
            "auxiliary_components": auxiliary_payload,
            "inter_residue_bonds": sorted(
                (remap[i], atom_i, remap[j], atom_j, kind)
                if (remap[i], atom_i) <= (remap[j], atom_j)
                else (remap[j], atom_j, remap[i], atom_i, kind)
                for i, atom_i, j, atom_j, kind in crosslinks
            ),
        }
        candidates.append((json.dumps(payload, sort_keys=True, separators=(",", ":")), payload))
    canonical_json, canonical_payload = min(candidates)
    digest = hashlib.sha256(canonical_json.encode()).hexdigest()
    return f"pep_{digest[:24]}", canonical_payload
