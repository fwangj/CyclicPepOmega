from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass

from .core.schema import CovalentBond, Residue, ResidueKey, StructureData, atom_distance, bond_dict
from .chemistry.residues import noncanonical_annotations

PEPTIDE_CN_MAX = 1.9
DISULFIDE_SS_MAX = 2.35


@dataclass
class DetectionConfig:
    min_residues: int = 2
    # v1 operational scope. This is a secondary peptide/protein boundary, never
    # a cyclicity criterion: a positive still requires a covalent ring closure.
    max_peptide_residues: int = 60
    max_ambiguous_residues: int = 80
    infer_disulfides: bool = True


def _ordered_pair(a: ResidueKey, b: ResidueKey) -> tuple[ResidueKey, ResidueKey]:
    return (a, b) if a < b else (b, a)


def _residue_sort_key(key: ResidueKey) -> tuple[object, ...]:
    try:
        seq: object = (0, int(key.seq))
    except ValueError:
        seq = (1, key.seq)
    return key.model, key.chain, seq, key.insertion


def reconstruct_bonds(data: StructureData, config: DetectionConfig) -> list[CovalentBond]:
    """Combine deposited links with conservative coordinate-based inference."""
    bonds = list(data.explicit_bonds)
    seen = {(b.left, b.left_atom, b.right, b.right_atom) for b in bonds}
    peptide = [r for r in data.residues.values() if r.peptide_like]

    # Infer C(i)-N(j) only within an entity/chain. Distance is deliberately
    # conservative; ambiguous close contacts are retained in the report later.
    for left in peptide:
        c = left.atom("C")
        if c is None:
            continue
        for right in peptide:
            if left.key == right.key or left.key.model != right.key.model:
                continue
            if left.key.chain != right.key.chain and left.entity_id != right.entity_id:
                continue
            n = right.atom("N")
            if n and atom_distance(c, n) <= PEPTIDE_CN_MAX:
                key = (left.key, "C", right.key, "N")
                reverse = (right.key, "N", left.key, "C")
                if key not in seen and reverse not in seen:
                    bonds.append(CovalentBond(left.key, "C", right.key, "N", "geometry", "peptide"))
                    seen.add(key)

    if config.infer_disulfides:
        sulfur = [(r, r.atom("SG")) for r in peptide if r.atom("SG")]
        for i, (left, a) in enumerate(sulfur):
            for right, b in sulfur[i + 1 :]:
                if left.key.model == right.key.model and atom_distance(a, b) <= DISULFIDE_SS_MAX:
                    key = (left.key, "SG", right.key, "SG")
                    reverse = (right.key, "SG", left.key, "SG")
                    if key not in seen and reverse not in seen:
                        bonds.append(CovalentBond(left.key, "SG", right.key, "SG", "geometry", "disulfide"))
                        seen.add(key)
    return bonds


def _components(nodes: set[ResidueKey], edges: list[CovalentBond]) -> list[set[ResidueKey]]:
    graph: dict[ResidueKey, set[ResidueKey]] = defaultdict(set)
    for edge in edges:
        if edge.left in nodes and edge.right in nodes:
            graph[edge.left].add(edge.right)
            graph[edge.right].add(edge.left)
    remaining = set(nodes)
    result = []
    while remaining:
        root = min(remaining)
        comp, queue = set(), deque([root])
        while queue:
            node = queue.popleft()
            if node in comp:
                continue
            comp.add(node)
            queue.extend(graph[node] - comp)
        remaining -= comp
        result.append(comp)
    return result


def _classify(edges: list[CovalentBond], residues: dict[ResidueKey, Residue]) -> tuple[str, list[CovalentBond]]:
    closure = []
    peptide_edges = []
    for edge in edges:
        atom_pair = {edge.left_atom.upper(), edge.right_atom.upper()}
        if atom_pair == {"C", "N"}:
            peptide_edges.append(edge)
        else:
            closure.append(edge)
    nodes = {e.left for e in edges} | {e.right for e in edges}
    # A closed backbone has one incoming and one outgoing peptide bond/residue.
    incoming = defaultdict(int)
    outgoing = defaultdict(int)
    for edge in peptide_edges:
        if edge.left_atom.upper() == "C":
            outgoing[edge.left] += 1
            incoming[edge.right] += 1
        else:
            outgoing[edge.right] += 1
            incoming[edge.left] += 1
    if nodes and all(incoming[n] >= 1 and outgoing[n] >= 1 for n in nodes):
        ordered = sorted(nodes, key=_residue_sort_key)
        first, last = ordered[0], ordered[-1]
        closures = [e for e in peptide_edges if (e.left == last and e.right == first) or (e.right == last and e.left == first)]
        return "head-to-tail", closures or peptide_edges
    disulfides = [e for e in closure if e.kind == "disulfide" or {e.left_atom.upper(), e.right_atom.upper()} == {"SG"}]
    if disulfides:
        return "disulfide-constrained", disulfides
    for edge in closure:
        backbone_left = edge.left_atom.upper() in {"N", "C"}
        backbone_right = edge.right_atom.upper() in {"N", "C"}
        if backbone_left != backbone_right:
            return "head-to-side-chain", [edge]
    if closure:
        return "side-chain-to-side-chain", closure
    return "noncyclic", []


def analyze_structure(data: StructureData, config: DetectionConfig | None = None) -> dict[str, object]:
    config = config or DetectionConfig()
    bonds = reconstruct_bonds(data, config)
    peptide_nodes = {key for key, residue in data.residues.items() if residue.peptide_like}
    peptide_bonds = [b for b in bonds if b.left in peptide_nodes and b.right in peptide_nodes]
    components = _components(peptide_nodes, peptide_bonds)
    peptides = []
    ambiguous = []
    entities = []
    for component in components:
        if len(component) < config.min_residues:
            continue
        edges = [b for b in peptide_bonds if b.left in component and b.right in component]
        kind, closure = _classify(edges, data.residues)
        ordered = sorted(component, key=_residue_sort_key)
        noncanonical = noncanonical_annotations(ordered, data.residues)
        disulfides = [b for b in edges if b.kind == "disulfide" or {b.left_atom.upper(), b.right_atom.upper()} == {"SG"}]
        component_chains = {k.chain for k in component}
        component_atoms = [atom for k in component for atom in data.residues[k].atoms.values()]
        target_chains = set()
        for key, residue in data.residues.items():
            if key in component or key.model != ordered[0].model or key.chain in component_chains:
                continue
            nucleotide_like = "P" in residue.atoms and any(name in residue.atoms for name in ("C4'", "C4*"))
            if not (residue.peptide_like or nucleotide_like):
                continue
            if any(atom_distance(a, b) <= 5.0 for a in component_atoms for b in residue.atoms.values()):
                target_chains.add((key.chain, "nucleic_acid" if nucleotide_like else "protein_or_peptide"))
        entities.append({
            "model_number": ordered[0].model,
            "chain_ids": sorted(component_chains),
            "entity_ids": sorted({data.residues[k].entity_id for k in component if data.residues[k].entity_id}),
            "residue_count": len(component),
            "has_detected_cycle": kind != "noncyclic",
        })
        payload = {
            "instance_id": f"{data.structure_id}:model{ordered[0].model}:{ordered[0].chain}:{ordered[0].seq}-{ordered[-1].seq}",
            "model_number": ordered[0].model,
            "chains": sorted({k.chain for k in component}),
            "entity_ids": sorted({data.residues[k].entity_id for k in component if data.residues[k].entity_id}),
            "residue_count": len(component),
            "residues": [
                {"original_id": k.label(), "name": data.residues[k].name, "internal_index": i + 1}
                for i, k in enumerate(ordered)
            ],
            "sequence": [data.residues[k].name for k in ordered],
            "noncanonical_residues": noncanonical,
            "cyclization_type": kind,
            "cyclization_bonds": [bond_dict(b) for b in closure],
            "disulfide_bonds": [bond_dict(b) for b in disulfides],
            "bound_targets": [
                {"chain_id": chain, "macromolecule_type": kind, "assignment": "contact-within-5A"}
                for chain, kind in sorted(target_chains)
            ],
            "warnings": (["noncanonical residue chemistry requires chem_comp-level annotation"] if noncanonical else []),
            "all_covalent_bonds": [bond_dict(b) for b in edges],
        }
        if kind != "noncyclic" and len(component) > config.max_peptide_residues:
            # Avoid classifying ordinary disulfide-rich proteins as cyclic
            # peptides. Preserve the entity inventory but not as a candidate.
            continue
        if kind == "noncyclic":
            if len(component) <= config.max_ambiguous_residues:
                ambiguous.append({**payload, "reason": "peptide-like covalent component has no detected ring closure"})
        else:
            peptides.append(payload)
    return {
        "schema_version": "0.1.0",
        "structure": {"id": data.structure_id, **data.metadata},
        "peptide_like_entities": entities,
        "cyclic_peptides": peptides,
        "ambiguous_or_noncyclic_candidates": ambiguous,
        "summary": {"cyclic_peptide_count": len(peptides), "candidate_count": len(peptides) + len(ambiguous)},
    }
