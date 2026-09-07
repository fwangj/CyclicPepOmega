from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass(frozen=True)
class AtomMapping:
    left_residue_index: int
    right_residue_index: int
    atom_pairs: tuple[tuple[str, str], ...]
    mapping_basis: str


@dataclass(frozen=True)
class ModificationEvent:
    event_type: str
    left_residue_index: int
    right_residue_index: int
    before: dict[str, object]
    after: dict[str, object]
    mapped_atoms: tuple[tuple[str, str], ...]
    evidence_status: str


@dataclass
class MatchedAnaloguePair:
    left_identity_id: str
    right_identity_id: str
    residue_mapping: list[tuple[int, int]]
    atom_mappings: list[AtomMapping]
    modifications: list[ModificationEvent]
    mapping_status: str
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def build_matched_analogue_pair(left: dict[str, object], right: dict[str, object]) -> MatchedAnaloguePair:
    """Create an explicit candidate mapping; it does not assert medicinal-chemistry relatedness."""
    left_residues = left["residues"]
    right_residues = right["residues"]
    if len(left_residues) != len(right_residues):
        raise ValueError("v1 analogue mapping requires equal normalized monomer counts")
    n = len(left_residues)
    left_annotations = left["residue_annotations"]
    right_annotations = right["residue_annotations"]
    # Select the cyclic registration with the fewest component/stereo/N-methyl
    # edits. This proposes a mapping only; downstream curation must accept it.
    candidates = []
    for shift in range(n):
        mapping = [(i, (i + shift) % n) for i in range(n)]
        cost = sum(
            left_residues[i]["name"] != right_residues[j]["name"]
            or left_annotations[i]["stereochemistry"] != right_annotations[j]["stereochemistry"]
            or left_annotations[i]["n_methylated_backbone_n"] != right_annotations[j]["n_methylated_backbone_n"]
            for i, j in mapping
        )
        candidates.append((cost, shift, mapping))
    _, shift, mapping = min(candidates)
    atom_mappings = []
    events = []
    for left_i, right_i in mapping:
        left_nodes = {
            node["atom_id"] for node in left["molecular_graph"]["nodes"]
            if node.get("residue_index") == left_i + 1 and str(node.get("element", "")).upper() != "H"
        }
        right_nodes = {
            node["atom_id"] for node in right["molecular_graph"]["nodes"]
            if node.get("residue_index") == right_i + 1 and str(node.get("element", "")).upper() != "H"
        }
        common_atoms = tuple((atom, atom) for atom in sorted(left_nodes & right_nodes))
        atom_mappings.append(AtomMapping(left_i + 1, right_i + 1, common_atoms, "same CCD atom_id intersection"))
        before, after = left_annotations[left_i], right_annotations[right_i]
        changes = []
        if left_residues[left_i]["name"] != right_residues[right_i]["name"]:
            changes.append("component_substitution")
        if before["stereochemistry"] != after["stereochemistry"]:
            changes.append("stereochemistry_change")
        if before["n_methylated_backbone_n"] != after["n_methylated_backbone_n"]:
            changes.append("backbone_n_methylation_change")
        for event_type in changes:
            events.append(ModificationEvent(
                event_type, left_i + 1, right_i + 1,
                {"component_id": left_residues[left_i]["name"], "stereochemistry": before["stereochemistry"], "n_methylated": before["n_methylated_backbone_n"]},
                {"component_id": right_residues[right_i]["name"], "stereochemistry": after["stereochemistry"], "n_methylated": after["n_methylated_backbone_n"]},
                common_atoms, "proposed_from_CCD_atom_ids",
            ))
    status = "candidate_single_edit" if len(events) == 1 else "candidate_multi_edit" if events else "exact_identity"
    return MatchedAnaloguePair(
        left["chemical_identity_id"], right["chemical_identity_id"],
        [(i + 1, j + 1) for i, j in mapping], atom_mappings, events, status,
        [
            f"cyclic registration shift={shift}",
            "mapping is a transparent candidate and requires chemical curation before matched-pair mining",
            "bond-change and unequal-length mappings are not implemented",
        ],
    )
