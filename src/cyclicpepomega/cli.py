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
    validate.add_argument("-o", "--output", type=Path, help="write validation summary to a new file")
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
        else:
            payload = validate_manifest(args.manifest, _analyze)
        _write(payload, args.output)
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"cyclicpepomega: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
