from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files


@lru_cache(maxsize=1)
def _registry() -> dict[tuple[str, str], dict[str, str]]:
    path = files("cyclicpepomega.data").joinpath("entity_scope_v1.json")
    payload = json.loads(path.read_text())
    return {
        (item["pdb_id"].upper(), item["entity_id"]): item
        for item in payload["annotations"]
    }


def disulfide_entity_scope(structure_id: str, entity_id: str | None) -> dict[str, str]:
    item = _registry().get((structure_id.upper(), entity_id or ""))
    if item:
        return item
    return {
        "scope": "unresolved",
        "evidence": "pure disulfide topology cannot distinguish peptide from disulfide-rich protein",
        "source": "cpo-entity-scope-1 policy",
    }
