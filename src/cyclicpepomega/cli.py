from __future__ import annotations

import argparse
import json
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

from .detect import analyze_structure
from .io.structures import parse_structure
from .validation import validate_manifest
from .conformation.geometry import compare_backbones
from .comparison.modification import build_matched_analogue_pair


def _download_pdb(identifier: str, destination: Path) -> Path:
    identifier = identifier.upper()
    url = f"https://files.rcsb.org/download/{identifier}.cif"
    target = destination / f"{identifier}.cif"
    try:
        with urllib.request.urlopen(url, timeout=30) as response, target.open("wb") as handle:
            handle.write(response.read())
    except urllib.error.URLError as exc:
        raise RuntimeError(f"failed to download {identifier} from RCSB: {exc}") from exc
    return target


def _analyze(source: str) -> dict[str, object]:
    path = Path(source)
    if path.exists():
        return analyze_structure(parse_structure(path))
    if len(source) == 4 and source.isalnum():
        with tempfile.TemporaryDirectory(prefix="cycatlas-") as tmp:
            downloaded = _download_pdb(source, Path(tmp))
            return analyze_structure(parse_structure(downloaded, source))
    raise ValueError(f"not a local structure file or four-character PDB ID: {source}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cyclicpepomega", description="Detect cyclic peptides using deposited and inferred covalent connectivity")
    sub = parser.add_subparsers(dest="command", required=True)
    analyze = sub.add_parser("analyze", help="analyze a PDB ID or local PDB/mmCIF file")
    analyze.add_argument("source")
    analyze.add_argument("-o", "--output", type=Path, help="write JSON to a new file")
    batch = sub.add_parser("batch", help="analyze PDB IDs/files listed one per line")
    batch.add_argument("list_file", type=Path)
    batch.add_argument("-o", "--output", type=Path, help="write JSON array to a new file")
    validate = sub.add_parser("validate", help="run a public-structure validation manifest")
    validate.add_argument("manifest", type=Path)
    validate.add_argument("--structure-dir", type=Path, help="use cached <PDB_ID>.cif files from this directory")
    validate.add_argument("-o", "--output", type=Path, help="write validation summary to a new file")
    compare = sub.add_parser("compare", help="compare experimental conformations with the same chemical identity")
    compare.add_argument("source", help="first PDB ID or structure file")
    compare.add_argument("other", nargs="?", help="optional second PDB ID or structure file")
    compare.add_argument("-o", "--output", type=Path, help="write comparison JSON to a new file")
    map_mod = sub.add_parser("map-modification", help="propose an atom-mapped analogue change and compare conformations")
    map_mod.add_argument("source", help="first PDB ID or structure file")
    map_mod.add_argument("other", help="second PDB ID or structure file")
    map_mod.add_argument("-o", "--output", type=Path, help="write mapping JSON to a new file")
    return parser


def _write(payload: object, output: Path | None) -> None:
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if output:
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing output: {output}")
        output.write_text(text)
    else:
        sys.stdout.write(text)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "analyze":
            payload = _analyze(args.source)
        elif args.command == "batch":
            sources = [line.strip() for line in args.list_file.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
            payload = [_analyze(source) for source in sources]
        elif args.command == "validate":
            def analyze_validation(source: str) -> dict[str, object]:
                if args.structure_dir and len(source) == 4 and source.isalnum():
                    cached = args.structure_dir / f"{source.upper()}.cif"
                    if cached.exists():
                        return _analyze(str(cached))
                return _analyze(source)
            payload = validate_manifest(args.manifest, analyze_validation)
        elif args.command == "compare":
            left_report = _analyze(args.source)
            right_report = _analyze(args.other) if args.other else left_report
            left_by_identity: dict[str, list[dict[str, object]]] = {}
            right_by_identity: dict[str, list[dict[str, object]]] = {}
            for item in left_report["cyclic_peptides"]:
                left_by_identity.setdefault(item["chemical_identity_id"], []).append(item)
            for item in right_report["cyclic_peptides"]:
                right_by_identity.setdefault(item["chemical_identity_id"], []).append(item)
            shared = sorted(set(left_by_identity) & set(right_by_identity))
            comparisons = []
            for identity_id in shared:
                if args.other:
                    pairs = [(left_by_identity[identity_id][0], right_by_identity[identity_id][0])]
                else:
                    observations = left_by_identity[identity_id]
                    pairs = [(observations[0], observations[1])] if len(observations) >= 2 else []
                for left, right in pairs:
                    comparisons.append({
                        "chemical_identity_id": identity_id,
                        "left_instance_id": left["instance_id"],
                        "right_instance_id": right["instance_id"],
                        **compare_backbones(left, right),
                    })
            if not comparisons:
                raise ValueError("no pair of cyclic-peptide observations with a shared chemical identity")
            payload = {
                "schema_version": "0.2.0",
                "sources": [args.source] + ([args.other] if args.other else []),
                "comparisons": comparisons,
            }
        else:
            left_report, right_report = _analyze(args.source), _analyze(args.other)
            if not left_report["cyclic_peptides"] or not right_report["cyclic_peptides"]:
                raise ValueError("both sources must contain a detected cyclic peptide")
            left, right = left_report["cyclic_peptides"][0], right_report["cyclic_peptides"][0]
            mapping = build_matched_analogue_pair(left, right)
            payload = {
                "schema_version": "0.2.0",
                "sources": [args.source, args.other],
                "matched_analogue_candidate": mapping.to_dict(),
                "conformational_difference": compare_backbones(left, right, allow_component_substitutions=True),
                "readiness": "requires expert chemical curation before inclusion in a Modification-to-Conformation dataset",
            }
        _write(payload, args.output)
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"cyclicpepomega: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
