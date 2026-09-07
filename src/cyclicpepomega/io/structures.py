from __future__ import annotations

from pathlib import Path

from ..core.schema import Atom, CovalentBond, Residue, ResidueKey, StructureData


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
    structure = gemmi.read_structure(str(path))
    residues: dict[ResidueKey, Residue] = {}
    for model_index, model in enumerate(structure, 1):
        for chain in model:
            chain_id = chain.name or "_"
            for source_residue in chain:
                key = ResidueKey(model_index, chain_id, str(source_residue.seqid.num), source_residue.seqid.icode.strip())
                residue = Residue(key, source_residue.name, str(source_residue.entity_id) if source_residue.entity_id else None)
                for source_atom in source_residue:
                    alt = source_atom.altloc.strip("\x00 ")
                    if alt not in {"", "A", "1"}:
                        continue
                    pos = source_atom.pos
                    residue.atoms.setdefault(source_atom.name.strip(), Atom(source_atom.name.strip(), source_atom.element.name, (pos.x, pos.y, pos.z), alt))
                residues[key] = residue
    # gemmi resolves both _struct_conn and legacy LINK records to connections.
    bonds = []
    for conn in structure.connections:
        p1, p2 = conn.partner1, conn.partner2
        k1 = ResidueKey(1, p1.chain_name or "_", str(p1.res_id.seqid.num), p1.res_id.seqid.icode.strip())
        k2 = ResidueKey(1, p2.chain_name or "_", str(p2.res_id.seqid.num), p2.res_id.seqid.icode.strip())
        if k1 in residues and k2 in residues:
            kind = "disulfide" if "disulf" in str(conn.type).lower() else "covalent"
            bonds.append(CovalentBond(k1, p1.atom_name, k2, p2.atom_name, "struct_conn", kind))
    method = None
    if hasattr(structure, "info"):
        try:
            method = structure.info["_exptl.method"]
        except (KeyError, TypeError):
            pass
    metadata = {"format": "mmcif", "experimental_method": method, "resolution_angstrom": structure.resolution or None}
    return StructureData((structure_id or structure.name or path.stem).upper(), residues, bonds, metadata)


def parse_structure(path: str | Path, structure_id: str | None = None) -> StructureData:
    path = Path(path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return parse_mmcif(path, structure_id)
    return parse_pdb(path, structure_id)
