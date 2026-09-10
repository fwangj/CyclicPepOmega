# CyclicPepOmega Agent Instructions

This repository follows the shared development policy maintained in `fwangj/OmegaCanon`. This local file is authoritative for routine Codex work; runtime access to OmegaCanon is not required.

## Sandbox Autonomy Principle

**Use maximum practical autonomy inside the disposable sandbox and least privilege outside it.**

Work autonomously on routine and reversible repository-local tasks. Without routine user confirmation, Codex may install packages and scientific software, create project environments, run shell commands, compile software, download clearly licensed public test/scientific data, run tests and benchmarks, inspect logs, diagnose failures, retry failed steps, and choose reasonable alternative implementations.

Do not stop at the first ordinary sandbox-local failure. Diagnose and attempt safe alternatives first. Prefer a minimal reproducible base environment with autonomous on-demand installation.

Ask primarily when an action crosses the sandbox boundary, requires sensitive credentials or unavailable information, creates significant external cost, affects protected persistent resources, or is destructive and not safely reversible.

## Persistence and Git

Treat the Codex Cloud sandbox as disposable. Create frequent meaningful commits. Verify that valuable work reaches GitHub or another approved persistent store before abandoning a sandbox. A sandbox-local commit alone is not proof of remote persistence.

Never force-push, rewrite published history, delete remote branches/tags, or perform repository administration unless explicitly authorized.

## Data, licensing, and secrets

Never commit credentials, tokens, passwords, private keys, sensitive client information, proprietary application data, downloaded structures that should remain external, or inappropriate large generated datasets. Preserve provenance and compatible licensing for public scientific data and dependencies.

Broad sandbox autonomy does not grant broad external authority. Do not access local Windows computers, home Linux systems, NAS storage, production services, or unrelated repositories unless explicitly authorized.

## CyclicPepOmega scope

Preserve the public-data-driven scientific scope and the distinction between implemented behavior, development-stage features, proposed capabilities, and validated evidence. Do not fabricate experimental values, causal conclusions, validation results, or performance claims. Maintain deterministic tests and chemical/structural provenance.

Use only public structures, public annotations, synthetic fixtures, and repository content. Keep downloaded coordinate caches, generated analysis outputs, local environments, and private working reports outside Git. Preserve the frozen `v0.2.0` tag.

## Setup and verification

From the repository root, use Python 3.10 or newer:

```bash
python -m pip install -e '.[gemmi,test]'
pytest
```

The maintained validation manifests require public PDB mmCIF files named `<PDB_ID>.cif` in a caller-supplied cache:

```bash
cyclicpepomega validate validation/public_set.json --structure-dir /path/to/public/mmcif-cache
cyclicpepomega validate validation/public_set_v2.json --structure-dir /path/to/public/mmcif-cache
```

Do not commit downloaded coordinates or generated JSON reports, and do not initiate full-PDB ingestion unless explicitly assigned. Add focused tests for behavior changes and run the full test suite before handoff. Executable scripts must include an appropriate shebang.

Keep missing and ambiguous states explicit. Do not equate asymmetric-unit contacts with biological interfaces, absent contacts with an apo state, or approximate IMHB/SASA/contact calculations with validated energetic, affinity, permeability, or selectivity predictions. Candidate analogue mappings require independent scientific review.
