import json

from cyclicpepomega.validation import validate_manifest


def test_validation_manifest_records_mismatch(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"cases": [{
        "pdb_id": "TEST", "source": "unused",
        "expected": {"cyclic_instance_count": 1},
    }]}))
    report = validate_manifest(manifest, lambda _: {
        "cyclic_peptides": [], "ambiguous_or_noncyclic_candidates": [],
        "summary": {},
    })
    assert report["totals"]["failed"] == 1
    assert report["totals"]["false_negatives"] == 1
    assert report["cases"][0]["status"] == "failed"
