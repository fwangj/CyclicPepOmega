from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass

from .core.schema import CovalentBond, Residue, ResidueKey, StructureData, atom_distance, bond_dict
from .chemistry.residues import noncanonical_annotations
from .chemistry.identity import atom_level_graph, annotate_residue, chemical_identity, modification_annotations
from .chemistry.scope import disulfide_entity_scope
from .conformation.geometry import describe_conformation

PEPTIDE_CN_MAX = 1.9
DISULFIDE_SS_MAX = 2.35


@dataclass
class DetectionConfig:
    min_residues: int = 2
    max_ambiguous_residues: int = 80
    min_macrocycle_coverage: float = 0.5
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

    # entity_poly_seq establishes ordinary adjacency. Coordinates identify the
    # observed atom endpoints, but distance is not the evidence source here.
    by_polymer: dict[tuple[int, str, str], list[Residue]] = defaultdict(list)
    for residue in peptide:
        if residue.entity_id and residue.label_seq_id:
            by_polymer[(residue.key.model, residue.entity_id, residue.label_asym_id or residue.key.chain)].append(residue)
    for group in by_polymer.values():
        try:
            ordered_group = sorted(group, key=lambda residue: int(residue.label_seq_id or ""))
        except ValueError:
            continue
        for left, right in zip(ordered_group, ordered_group[1:]):
            if int(right.label_seq_id or "0") - int(left.label_seq_id or "0") != 1:
                continue
            key = (left.key, "C", right.key, "N")
            reverse = (right.key, "N", left.key, "C")
            if key not in seen and reverse not in seen:
                bonds.append(CovalentBond(
                    left.key, "C", right.key, "N", "polymer_connectivity",
                    "peptide", f"entity_poly_seq:{left.entity_id}:{left.label_seq_id}-{right.label_seq_id}",
                    "deposited_polymer_sequence",
                ))
                seen.add(key)

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
                    bonds.append(CovalentBond(left.key, "C", right.key, "N", "coordinate_inferred", "peptide", None, "inferred"))
                    seen.add(key)

    if config.infer_disulfides:
        sulfur = [(r, r.atom("SG")) for r in peptide if r.atom("SG")]
        for i, (left, a) in enumerate(sulfur):
            for right, b in sulfur[i + 1 :]:
                if left.key.model == right.key.model and atom_distance(a, b) <= DISULFIDE_SS_MAX:
                    key = (left.key, "SG", right.key, "SG")
                    reverse = (right.key, "SG", left.key, "SG")
                    if key not in seen and reverse not in seen:
                        bonds.append(CovalentBond(left.key, "SG", right.key, "SG", "coordinate_inferred", "disulfide", None, "inferred"))
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


def _atom_element(edge: CovalentBond, left: bool, residues: dict[ResidueKey, Residue]) -> str:
    key, name = (edge.left, edge.left_atom) if left else (edge.right, edge.right_atom)
    atom = residues[key].atom(name)
    return atom.element.upper() if atom else name.strip("0123456789")[:1].upper()


def _classify(
    edges: list[CovalentBond], residues: dict[ResidueKey, Residue],
    min_coverage: float, entity_lengths: dict[str, int] | None = None,
) -> tuple[str, list[CovalentBond], list[str], float | None]:
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
        selected = closures or peptide_edges
        return "head-to-tail", selected, [], 1.0
    disulfides = [e for e in closure if e.kind == "disulfide" or {e.left_atom.upper(), e.right_atom.upper()} == {"SG"}]
    ordered = sorted(nodes, key=_residue_sort_key)
    index = {key: i for i, key in enumerate(ordered)}
    entity_ids = {residues[key].entity_id for key in nodes if residues[key].entity_id}
    denominator = len(ordered)
    if len(entity_ids) == 1 and entity_lengths:
        denominator = entity_lengths.get(next(iter(entity_ids)), denominator) or denominator
    coverages = []
    for edge in closure:
        try:
            span = abs(int(edge.left.seq) - int(edge.right.seq)) + 1
        except ValueError:
            span = abs(index[edge.left] - index[edge.right]) + 1
        coverages.append(span / denominator)
    coverage = max(coverages, default=None)
    if closure and len(disulfides) == len(closure) and coverage is not None and coverage < min_coverage:
        return "ambiguous_unresolved", closure, [
            "internal crosslink does not span enough of the observed polymer to establish peptide-macrocycle scope"
        ], coverage
    if disulfides:
        if len(disulfides) == len(closure):
            return "disulfide-constrained", disulfides, [], coverage
    if len(closure) > 1:
        return "multiple-crosslink peptide", closure, [], coverage
    for edge in closure:
        left_element = _atom_element(edge, True, residues)
        right_element = _atom_element(edge, False, residues)
        elements = {left_element, right_element}
        backbone_left = edge.left_atom.upper() in {"N", "C", "O"}
        backbone_right = edge.right_atom.upper() in {"N", "C", "O"}
        if backbone_left != backbone_right:
            return "head-to-side-chain", [edge], [], coverage
        if elements == {"C", "N"}:
            return "side-chain-to-side-chain lactam", [edge], [], coverage
        if elements == {"C", "S"}:
            return "thioether/sulfur crosslink", [edge], [], coverage
        if elements == {"C", "O"}:
            return "depsipeptide/ester macrocycle", [edge], [], coverage
    if closure:
        return "other covalent macrocycle", closure, [], coverage
    return "linear peptide", [], [], None


def analyze_structure(data: StructureData, config: DetectionConfig | None = None) -> dict[str, object]:
    config = config or DetectionConfig()
    bonds = reconstruct_bonds(data, config)
    peptide_nodes = {key for key, residue in data.residues.items() if residue.peptide_like}
    attached: dict[ResidueKey, set[ResidueKey]] = defaultdict(set)
    for bond in bonds:
        if bond.left in peptide_nodes and bond.right not in peptide_nodes:
            attached[bond.right].add(bond.left)
        elif bond.right in peptide_nodes and bond.left not in peptide_nodes:
            attached[bond.left].add(bond.right)
    connector_nodes = {key for key, neighbors in attached.items() if len(neighbors) >= 2}
    graph_nodes = peptide_nodes | connector_nodes
    graph_bonds = [b for b in bonds if b.left in graph_nodes and b.right in graph_nodes]
    components = _components(graph_nodes, graph_bonds)
    peptides = []
    ambiguous = []
    entities = []
    for component in components:
        peptide_component = component & peptide_nodes
        if len(peptide_component) < config.min_residues:
            continue
        edges = [b for b in graph_bonds if b.left in component and b.right in component]
        kind, closure, classification_warnings, macrocycle_coverage = _classify(
            edges, data.residues, config.min_macrocycle_coverage,
            {entity_id: len(entity.sequence_components) for entity_id, entity in data.entities.items()},
        )
        ordered = sorted(peptide_component, key=_residue_sort_key)
        auxiliary = sorted(component - peptide_nodes, key=lambda key: (data.residues[key].name, _residue_sort_key(key)))
        scope_evidence = {"scope": "chemically_confirmed", "evidence": "non-disulfide macrocycle topology", "source": "cpo-cyclization-1"}
        if kind == "disulfide-constrained":
            entity_ids_for_scope = {data.residues[key].entity_id for key in ordered if data.residues[key].entity_id}
            entity_id_for_scope = next(iter(entity_ids_for_scope)) if len(entity_ids_for_scope) == 1 else None
            scope_evidence = disulfide_entity_scope(data.structure_id, entity_id_for_scope)
            if scope_evidence["scope"] != "cyclic_peptide":
                kind = "ambiguous_unresolved"
                classification_warnings.append(
                    "disulfide macrocycle detected, but public entity semantics do not establish cyclic-peptide scope"
                )
        noncanonical = noncanonical_annotations(ordered, data.residues)
        residue_annotations = [
            annotate_residue(i + 1, key, data.residues[key], data.components.get(data.residues[key].name))
            for i, key in enumerate(ordered)
        ]
        identity_id, identity_payload = chemical_identity(
            data, ordered, edges, residue_annotations,
            kind not in {"linear peptide", "ambiguous_unresolved"},
            auxiliary,
        )
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
            "residue_count": len(peptide_component),
            "has_detected_cycle": kind not in {"linear peptide", "ambiguous_unresolved"},
        })
        payload = {
            "instance_id": f"{data.structure_id}:model{ordered[0].model}:{ordered[0].chain}:{ordered[0].seq}-{ordered[-1].seq}",
            "chemical_identity_id": identity_id,
            "chemical_identity_version": identity_payload["version"],
            "molecular_graph": atom_level_graph(data, ordered, edges, auxiliary),
            "covalent_connector_components": [
                {"component_id": data.residues[key].name, "original_id": key.label()}
                for key in auxiliary
            ],
            "model_number": ordered[0].model,
            "chains": sorted({k.chain for k in component}),
            "entity_ids": sorted({data.residues[k].entity_id for k in component if data.residues[k].entity_id}),
            "residue_count": len(ordered),
            "residues": [
                {"original_id": k.label(), "name": data.residues[k].name, "internal_index": i + 1}
                for i, k in enumerate(ordered)
            ],
            "sequence": [data.residues[k].name for k in ordered],
            "noncanonical_residues": noncanonical,
            "residue_annotations": residue_annotations,
            "modifications": modification_annotations(residue_annotations),
            "cyclization_type": kind,
            "cyclization_bonds": [bond_dict(b) for b in closure],
            "cyclization": {
                "classification": kind,
                "closure_atoms": [bond_dict(b) for b in closure],
                "evidence_status": (
                    "deposited" if closure and all(b.evidence_status == "deposited" for b in closure)
                    else "coordinate_inferred" if closure and all(b.evidence_status == "inferred" for b in closure)
                    else "mixed" if closure else "none"
                ),
                "macrocycle_coverage": macrocycle_coverage,
                "warnings": classification_warnings,
            },
            "peptide_scope": scope_evidence,
            "disulfide_bonds": [bond_dict(b) for b in disulfides],
            "bound_targets": [
                {"chain_id": chain, "macromolecule_type": kind, "assignment": "contact-within-5A"}
                for chain, kind in sorted(target_chains)
            ],
            "warnings": classification_warnings + (["one or more component definitions are missing"] if any(data.residues[k].name not in data.components for k in ordered) else []),
            "all_covalent_bonds": [bond_dict(b) for b in edges],
        }
        payload["conformation"] = describe_conformation(
            data.structure_id, identity_id, ordered, data.residues, edges, bool(target_chains)
        )
        if kind in {"linear peptide", "ambiguous_unresolved"}:
            if len(peptide_component) <= config.max_ambiguous_residues:
                reason = "no detected ring closure" if kind == "linear peptide" else "crosslinked polymer has unresolved peptide/macrocycle scope"
                ambiguous.append({**payload, "reason": reason})
        else:
            peptides.append(payload)
    identity_ids = {x["chemical_identity_id"] for x in peptides}
    return {
        "schema_version": "0.2.0",
        "structure": {"id": data.structure_id, **data.metadata},
        "peptide_like_entities": entities,
        "cyclic_peptides": peptides,
        "ambiguous_or_noncyclic_candidates": ambiguous,
        "summary": {
            "chemical_identity_count": len(identity_ids),
            "cyclic_peptide_instance_count": len(peptides),
            "cyclic_peptide_count": len(peptides),
            "candidate_count": len(peptides) + len(ambiguous),
            "coordinate_model_count": len({key.model for key in data.residues}),
        },
    }
