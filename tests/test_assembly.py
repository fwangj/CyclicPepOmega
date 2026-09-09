from cyclicpepomega.assembly import assembly_summary_for_structure
from cyclicpepomega.core.schema import StructureData


def test_assembly_summary_preserves_metadata_only_status():
    data = StructureData("TEST", {}, metadata={"biological_assemblies": [{"assembly_id": "1", "status": "metadata_only"}]})
    summary = assembly_summary_for_structure(data)
    assert summary["assembly_count"] == 1
    assert summary["method"]["version"] == "cpo-assembly-1"
    assert "asymmetric-unit" in summary["contact_context_policy"]
