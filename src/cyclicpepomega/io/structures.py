from __future__ import annotations

from pathlib import Path

from ..assembly import assembly_records_to_dicts, parse_mmcif_assembly_metadata
from ..core.schema import (
    Atom, ComponentAtom, ComponentBond, ComponentDefinition, CovalentBond,
    EntityDefinition, Residue, ResidueKey, StructureData,
)


def _clean(value: object) -> str | None:
    text = str(value).strip().strip("'").strip('"')
    return None if text in {"", ".", "?"} else text


def _table_rows(block: object, category: str) -> tuple[dict[str, int], object]:
    table = block.find_mmcif_category(category)
    tags = {str(tag).split(".")[-1]: index for index, tag in enumerate(table.tags)}
    return tags, table


def _value(row: object, tags: dict[str, int], name: str) -> str | None:
    return _clean(row[tags[name]]) if name in tags else None


def _element(atom_name: str, field: str = "") -> str:
    return (field.strip() or atom_name.strip().lstrip("0123456789")[:1]).upper()


def parse_pdb(path: str | Path, structure_id: str | None = None) -> StructureData:
    """Parse coordinates and deposited LINK/SSBOND records from legacy PDB."""
    path = Path(path)
    residues: dict[ResidueKey, Residue] = {}
    links_raw: list[tuple[str, str, str, str, str, str, str, str, str]] = []
    model = 1
    method = None
    resolution = None
    with path.open(errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                try: model = int(line[10:14])
                except ValueError: model += 1
            elif record in {"ATOM", "HETATM"}:
                alt = line[16:17].strip()
                if alt not in {"", "A", "1"}:
                    continue
                chain, seq, ins = line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip()
                key = ResidueKey(model, chain, seq, ins)
                residue = residues.setdefault(key, Residue(key, line[17:20].strip()))
                name = line[12:16].strip()
                try: xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError: continue
                residue.atoms.setdefault(name, Atom(name, _element(name, line[76:78]), xyz, alt))
            elif record == "LINK":
                links_raw.append((line[12:16].strip(), line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip(), line[42:46].strip(), line[51:52].strip() or "_", line[52:56].strip(), line[56:57].strip(), "LINK"))
            elif record == "SSBOND":
                links_raw.append(("SG", line[15:16].strip() or "_", line[17:21].strip(), line[21:22].strip(), "SG", line[29:30].strip() or "_", line[31:35].strip(), line[35:36].strip(), "SSBOND"))
            elif record == "EXPDTA": method = line[10:].strip()
            elif line.startswith("REMARK   2 RESOLUTION."):
                try: resolution = float(line[22:30])
                except ValueError: pass
    bonds = []
    for a1, c1, s1, i1, a2, c2, s2, i2, source in links_raw:
        k1, k2 = ResidueKey(1, c1, s1, i1), ResidueKey(1, c2, s2, i2)
        if k1 in residues and k2 in residues:
            kind = "disulfide" if source == "SSBOND" else "covalent"
            bonds.append(CovalentBond(k1, a1, k2, a2, source, kind))
    metadata = {"format": "pdb", "experimental_method": method, "resolution_angstrom": resolution}
    return StructureData((structure_id or path.stem).upper(), residues, bonds, metadata)


def parse_mmcif(path: str | Path, structure_id: str | None = None) -> StructureData:
    """Parse mmCIF with gemmi while preserving author residue identifiers."""
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("mmCIF parsing requires gemmi; install cyclicpepomega[gemmi]") from exc
    path = Path(path)
    document = gemmi.cif.read_file(str(path))
    block = document.sole_block()
    structure = gemmi.make_structure_from_block(block)
    residues: dict[ResidueKey, Residue] = {}
    atom_tags, atom_rows = _table_rows(block, "_atom_site.")
    label_lookup: dict[tuple[int, str, str, str], ResidueKey] = {}
    for row in atom_rows:
        alt = _value(row, atom_tags, "label_alt_id") or ""
        if alt not in {"", "A", "1"}:
            continue
        model = int(_value(row, atom_tags, "pdbx_PDB_model_num") or "1")
        label_asym = _value(row, atom_tags, "label_asym_id") or "_"
        auth_asym = _value(row, atom_tags, "auth_asym_id") or label_asym
        raw_label_seq = _value(row, atom_tags, "label_seq_id")
        label_seq = raw_label_seq or _value(row, atom_tags, "auth_seq_id") or "?"
        auth_seq = _value(row, atom_tags, "auth_seq_id") or label_seq
        insertion = _value(row, atom_tags, "pdbx_PDB_ins_code") or ""
        component_id = _value(row, atom_tags, "label_comp_id") or _value(row, atom_tags, "auth_comp_id") or "UNK"
        entity_id = _value(row, atom_tags, "label_entity_id")
        key = ResidueKey(model, auth_asym, auth_seq, insertion)
        residue = residues.setdefault(key, Residue(
            key, component_id, entity_id, {}, label_asym, label_seq, auth_asym, auth_seq
        ))
        label_atom = _value(row, atom_tags, "label_atom_id") or "?"
        auth_atom = _value(row, atom_tags, "auth_atom_id") or label_atom
        try:
            xyz = tuple(float(_value(row, atom_tags, axis) or "nan") for axis in ("Cartn_x", "Cartn_y", "Cartn_z"))
            occupancy = float(_value(row, atom_tags, "occupancy") or "1")
        except ValueError:
            continue
        element = _value(row, atom_tags, "type_symbol") or _element(auth_atom)
        residue.atoms.setdefault(auth_atom, Atom(auth_atom, element, xyz, alt, label_atom, auth_atom, occupancy))
        label_lookup[(model, label_asym, label_seq, component_id)] = key
        if raw_label_seq is None:
            label_lookup[(model, label_asym, "?", component_id)] = key

    components: dict[str, ComponentDefinition] = {}
    comp_tags, comp_rows = _table_rows(block, "_chem_comp.")
    for row in comp_rows:
        component_id = _value(row, comp_tags, "id")
        if component_id:
            components[component_id] = ComponentDefinition(
                component_id, _value(row, comp_tags, "type"),
                _value(row, comp_tags, "name"), _value(row, comp_tags, "formula"),
            )
    cca_tags, cca_rows = _table_rows(block, "_chem_comp_atom.")
    for row in cca_rows:
        component_id, atom_id = _value(row, cca_tags, "comp_id"), _value(row, cca_tags, "atom_id")
        if component_id and atom_id:
            definition = components.setdefault(component_id, ComponentDefinition(component_id))
            definition.atoms[atom_id] = ComponentAtom(
                atom_id, _value(row, cca_tags, "type_symbol") or "?",
                _value(row, cca_tags, "pdbx_stereo_config"),
                (_value(row, cca_tags, "pdbx_aromatic_flag") or "N").upper() == "Y",
            )
    ccb_tags, ccb_rows = _table_rows(block, "_chem_comp_bond.")
    for row in ccb_rows:
        component_id = _value(row, ccb_tags, "comp_id")
        atom1, atom2 = _value(row, ccb_tags, "atom_id_1"), _value(row, ccb_tags, "atom_id_2")
        if component_id and atom1 and atom2:
            definition = components.setdefault(component_id, ComponentDefinition(component_id))
            definition.bonds.append(ComponentBond(
                atom1, atom2, _value(row, ccb_tags, "value_order") or "unknown",
                (_value(row, ccb_tags, "pdbx_aromatic_flag") or "N").upper() == "Y",
                _value(row, ccb_tags, "pdbx_stereo_config"),
            ))

    entities: dict[str, EntityDefinition] = {}
    entity_tags, entity_rows = _table_rows(block, "_entity.")
    for row in entity_rows:
        entity_id = _value(row, entity_tags, "id")
        if entity_id:
            entities[entity_id] = EntityDefinition(
                entity_id, _value(row, entity_tags, "type"), None,
                _value(row, entity_tags, "pdbx_description"),
            )
    poly_tags, poly_rows = _table_rows(block, "_entity_poly.")
    for row in poly_rows:
        entity_id = _value(row, poly_tags, "entity_id")
        if entity_id:
            entities.setdefault(entity_id, EntityDefinition(entity_id)).polymer_type = _value(row, poly_tags, "type")
    seq_tags, seq_rows = _table_rows(block, "_entity_poly_seq.")
    for row in seq_rows:
        entity_id, monomer = _value(row, seq_tags, "entity_id"), _value(row, seq_tags, "mon_id")
        if entity_id and monomer:
            entities.setdefault(entity_id, EntityDefinition(entity_id)).sequence_components.append(monomer)

    # Read deposited covalent connections directly to preserve label IDs and the
    # connection record identifier. A deposited connection is replicated across
    # coordinate models only when both referenced residues exist in that model.
    bonds = []
    conn_tags, conn_rows = _table_rows(block, "_struct_conn.")
    models = sorted({key.model for key in residues})
    for row in conn_rows:
        connection_type = (_value(row, conn_tags, "conn_type_id") or "").lower()
        if not any(token in connection_type for token in ("covale", "disulf")):
            continue
        p1 = (_value(row, conn_tags, "ptnr1_label_asym_id") or "_", _value(row, conn_tags, "ptnr1_label_seq_id") or "?", _value(row, conn_tags, "ptnr1_label_comp_id") or "UNK")
        p2 = (_value(row, conn_tags, "ptnr2_label_asym_id") or "_", _value(row, conn_tags, "ptnr2_label_seq_id") or "?", _value(row, conn_tags, "ptnr2_label_comp_id") or "UNK")
        atom1 = _value(row, conn_tags, "ptnr1_label_atom_id") or "?"
        atom2 = _value(row, conn_tags, "ptnr2_label_atom_id") or "?"
        for model in models:
            k1, k2 = label_lookup.get((model, *p1)), label_lookup.get((model, *p2))
            if k1 and k2:
                kind = "disulfide" if "disulf" in connection_type else "covalent"
                bonds.append(CovalentBond(
                    k1, atom1, k2, atom2, "deposited_struct_conn", kind,
                    _value(row, conn_tags, "id"), "deposited",
                ))
    method = None
    if hasattr(structure, "info"):
        try:
            method = structure.info["_exptl.method"]
        except (KeyError, TypeError):
            pass
    metadata = {
        "format": "mmcif",
        "experimental_method": method,
        "resolution_angstrom": structure.resolution or None,
        "biological_assemblies": assembly_records_to_dicts(parse_mmcif_assembly_metadata(path)),
    }
    return StructureData(
        (structure_id or structure.name or path.stem).upper(), residues, bonds,
        metadata, components, entities,
    )


def parse_structure(path: str | Path, structure_id: str | None = None) -> StructureData:
    path = Path(path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return parse_mmcif(path, structure_id)
    return parse_pdb(path, structure_id)
