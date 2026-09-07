from __future__ import annotations

import json
from pathlib import Path
from typing import Callable


def validate_manifest(
    manifest_path: Path,
    analyze: Callable[[str], dict[str, object]],
) -> dict[str, object]:
    """Run transparent field-level regression checks from a public-data manifest."""
    manifest = json.loads(manifest_path.read_text())
    cases = []
    totals = {
        "passed": 0, "failed": 0, "false_positives": 0,
        "false_negatives": 0, "classification_errors": 0,
        "parsing_failures": 0,
    }
    for case in manifest["cases"]:
        result = {"pdb_id": case["pdb_id"], "status": "passed", "errors": []}
        try:
            report = analyze(case["source"])
            found = report["cyclic_peptides"]
            expected = case["expected"]
            if len(found) != expected["cyclic_instance_count"]:
                result["errors"].append(f"cyclic_instance_count: expected {expected['cyclic_instance_count']}, found {len(found)}")
                delta = len(found) - expected["cyclic_instance_count"]
                totals["false_positives" if delta > 0 else "false_negatives"] += abs(delta)
            for assertion in expected.get("contains", []):
                matches = [x for x in found if assertion.items() <= {
                    "chain_id": x["chains"][0] if len(x["chains"]) == 1 else None,
                    "length": x["residue_count"],
                    "cyclization_type": x["cyclization_type"],
                }.items()]
                if not matches:
                    result["errors"].append(f"missing expected cyclic peptide: {assertion}")
                    same_entity = [x for x in found if x["chains"] == [assertion["chain_id"]] and x["residue_count"] == assertion["length"]]
                    if same_entity:
                        totals["classification_errors"] += 1
            if expected.get("linear_candidate_chain"):
                chains = {c for x in report["ambiguous_or_noncyclic_candidates"] for c in x["chains"]}
                if expected["linear_candidate_chain"] not in chains:
                    result["errors"].append(f"missing linear candidate chain {expected['linear_candidate_chain']}")
        except Exception as exc:  # validation must record, not abort, parse failures
            result["status"] = "parsing_failure"
            result["errors"].append(f"{type(exc).__name__}: {exc}")
            totals["parsing_failures"] += 1
        if result["errors"] and result["status"] != "parsing_failure":
            result["status"] = "failed"
        totals["passed" if result["status"] == "passed" else "failed"] += 1
        cases.append(result)
    return {"schema_version": "0.1.0", "totals": totals, "cases": cases}
