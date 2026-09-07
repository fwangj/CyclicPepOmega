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
        "stereochemistry_errors": 0, "modification_annotation_errors": 0,
        "parsing_failures": 0,
    }
    for case in manifest["cases"]:
        result = {"pdb_id": case["pdb_id"], "status": "passed", "errors": []}
        try:
            report = analyze(case["source"])
            found = report["cyclic_peptides"]
            expected = case["expected"]
            expected_instances = expected.get("cyclic_instance_count", 0)
            if len(found) != expected_instances:
                result["errors"].append(f"cyclic_instance_count: expected {expected_instances}, found {len(found)}")
                delta = len(found) - expected_instances
                totals["false_positives" if delta > 0 else "false_negatives"] += abs(delta)
            expected_identities = expected.get("chemical_identity_count")
            found_identities = report["summary"].get("chemical_identity_count")
            if expected_identities is not None and found_identities != expected_identities:
                result["errors"].append(f"chemical_identity_count: expected {expected_identities}, found {found_identities}")
            expected_models = expected.get("coordinate_model_count")
            found_models = report["summary"].get("coordinate_model_count")
            if expected_models is not None and found_models != expected_models:
                result["errors"].append(f"coordinate_model_count: expected {expected_models}, found {found_models}")
            for assertion in expected.get("contains", []):
                core_assertion = {k: assertion[k] for k in ("chain_id", "length", "cyclization_type") if k in assertion}
                matches = [x for x in found if core_assertion.items() <= {
                    "chain_id": x["chains"][0] if len(x["chains"]) == 1 else None,
                    "length": x["residue_count"],
                    "cyclization_type": x["cyclization_type"],
                }.items()]
                if not matches:
                    result["errors"].append(f"missing expected cyclic peptide: {assertion}")
                    same_entity = [x for x in found if x["chains"] == [assertion.get("chain_id")] and x["residue_count"] == assertion.get("length")]
                    if same_entity:
                        totals["classification_errors"] += 1
                    continue
                item = matches[0]
                if "closure_atoms" in assertion:
                    observed = sorted(
                        sorted((bond["left_atom"], bond["right_atom"]))
                        for bond in item["cyclization_bonds"]
                    )
                    wanted = sorted(sorted(pair) for pair in assertion["closure_atoms"])
                    if observed != wanted:
                        result["errors"].append(f"closure atoms mismatch for {assertion.get('chain_id')}: expected {wanted}, found {observed}")
                        totals["classification_errors"] += 1
                if "d_residue_count" in assertion:
                    observed = sum(a["stereochemistry"] == "D" for a in item["residue_annotations"])
                    if observed != assertion["d_residue_count"]:
                        result["errors"].append(f"D residue count mismatch: expected {assertion['d_residue_count']}, found {observed}")
                        totals["stereochemistry_errors"] += 1
                if "n_methylated_count" in assertion:
                    observed = sum(a["n_methylated_backbone_n"] for a in item["residue_annotations"])
                    if observed != assertion["n_methylated_count"]:
                        result["errors"].append(f"N-methyl count mismatch: expected {assertion['n_methylated_count']}, found {observed}")
                        totals["modification_annotation_errors"] += 1
                if "noncanonical_count" in assertion:
                    observed = sum(not a["canonical_amino_acid"] for a in item["residue_annotations"])
                    if observed != assertion["noncanonical_count"]:
                        result["errors"].append(f"noncanonical count mismatch: expected {assertion['noncanonical_count']}, found {observed}")
                        totals["modification_annotation_errors"] += 1
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
    return {"schema_version": "0.2.0", "totals": totals, "cases": cases}
