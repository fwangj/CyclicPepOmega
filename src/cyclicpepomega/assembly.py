from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

from .core.schema import StructureData

ASSEMBLY_METHOD_VERSION = "cpo-assembly-1"


@dataclass(frozen=True)
class AssemblyRecord:
    assembly_id: str
    details: str | None = None
    method_details: str | None = None
    oligomeric_details: str | None = None
    oligomeric_count: int | None = None
    asym_ids: tuple[str, ...] = ()
    operator_expression: str | None = None
    status: str = "metadata_only"
    warnings: tuple[str, ...] = ()


def parse_mmcif_assembly_metadata(path: str | Path) -> list[AssemblyRecord]:
    """Parse mmCIF biological assembly categories without generating coordinates."""

    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("assembly parsing requires gemmi; install cyclicpepomega[gemmi]") from exc
    block = gemmi.cif.read_file(str(path)).sole_block()
    gen_by_id: dict[str, list[tuple[tuple[str, ...], str | None]]] = {}
    try:
        gen = block.find_mmcif_category("_pdbx_struct_assembly_gen.")
        tags = {str(tag).split(".")[-1]: index for index, tag in enumerate(gen.tags)}
        for row in gen:
            assembly_id = _clean(row[tags.get("assembly_id")]) if "assembly_id" in tags else None
            asym = _clean(row[tags.get("asym_id_list")]) if "asym_id_list" in tags else None
            oper = _clean(row[tags.get("oper_expression")]) if "oper_expression" in tags else None
            if assembly_id:
                gen_by_id.setdefault(assembly_id, []).append((tuple(part.strip() for part in (asym or "").split(",") if part.strip()), oper))
    except Exception:
        gen_by_id = {}
    records: list[AssemblyRecord] = []
    try:
        table = block.find_mmcif_category("_pdbx_struct_assembly.")
        tags = {str(tag).split(".")[-1]: index for index, tag in enumerate(table.tags)}
        for row in table:
            assembly_id = _clean(row[tags.get("id")]) if "id" in tags else None
            if not assembly_id:
                continue
            generated = gen_by_id.get(assembly_id, [])
            records.append(AssemblyRecord(
                assembly_id=assembly_id,
                details=_clean(row[tags.get("details")]) if "details" in tags else None,
                method_details=_clean(row[tags.get("method_details")]) if "method_details" in tags else None,
                oligomeric_details=_clean(row[tags.get("oligomeric_details")]) if "oligomeric_details" in tags else None,
                oligomeric_count=_optional_int(_clean(row[tags.get("oligomeric_count")])) if "oligomeric_count" in tags else None,
                asym_ids=tuple(sorted({asym for group, _ in generated for asym in group})),
                operator_expression=";".join(op for _, op in generated if op) or None,
                warnings=("coordinate transform generation is not implemented in cpo-assembly-1",),
            ))
    except Exception:
        return []
    return records


def assembly_summary_for_structure(data: StructureData) -> dict[str, object]:
    records = data.metadata.get("biological_assemblies", [])
    return {
        "method": {"name": "cyclicpepomega_assembly", "version": ASSEMBLY_METHOD_VERSION},
        "assembly_count": len(records) if isinstance(records, list) else 0,
        "contact_context_policy": "contacts are asymmetric-unit unless biological assembly metadata is available and explicitly linked",
        "biological_assemblies": records,
        "warnings": ["assembly coordinate transforms are metadata-only in this version"],
    }


def assembly_records_to_dicts(records: list[AssemblyRecord]) -> list[dict[str, object]]:
    return [asdict(record) for record in records]


def _clean(value: object) -> str | None:
    text = str(value).strip().strip("'").strip('"')
    return None if text in {"", ".", "?", "None"} else text


def _optional_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None
