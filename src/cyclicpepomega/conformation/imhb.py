from __future__ import annotations

from dataclasses import asdict, dataclass
from math import dist

from ..core.schema import Residue, ResidueKey

IMHB_METHOD_VERSION = "cpo-imhb-1"
DONOR_ELEMENTS = {"N", "O", "S"}
ACCEPTOR_ELEMENTS = {"N", "O", "S"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class IntramolecularHydrogenBond:
    donor_residue_id: str
    donor_residue_index: int
    donor_atom_id: str
    acceptor_residue_id: str
    acceptor_residue_index: int
    acceptor_atom_id: str
    donor_acceptor_distance_angstrom: float
    geometry_evidence: str
    atom_class: str
    warning: str | None = None

    def fingerprint_token(self) -> str:
        return f"{self.donor_residue_index}:{self.donor_atom_id}->{self.acceptor_residue_index}:{self.acceptor_atom_id}"


def classify_atom_pair(donor_atom: str, acceptor_atom: str) -> str:
    donor_backbone = donor_atom.upper() in BACKBONE_ATOMS
    acceptor_backbone = acceptor_atom.upper() in BACKBONE_ATOMS
    if donor_backbone and acceptor_backbone:
        return "backbone/backbone"
    if donor_backbone or acceptor_backbone:
        return "backbone/side-chain"
    return "side-chain/side-chain"


def intramolecular_hydrogen_bonds(
    keys: list[ResidueKey],
    residues: dict[ResidueKey, Residue],
    *,
    heavy_atom_cutoff: float = 3.5,
    min_sequence_separation: int = 2,
) -> dict[str, object]:
    """Infer cyclic-peptide IMHBs from heavy atoms without adding hydrogens."""

    observations: list[IntramolecularHydrogenBond] = []
    for left_index, left_key in enumerate(keys, start=1):
        left_residue = residues[left_key]
        for right_index, right_key in enumerate(keys, start=1):
            if abs(left_index - right_index) < min_sequence_separation:
                continue
            right_residue = residues[right_key]
            for donor in left_residue.atoms.values():
                if donor.element.upper() not in DONOR_ELEMENTS:
                    continue
                for acceptor in right_residue.atoms.values():
                    if acceptor.element.upper() not in ACCEPTOR_ELEMENTS:
                        continue
                    if donor.name == acceptor.name and left_key == right_key:
                        continue
                    distance = dist(donor.xyz, acceptor.xyz)
                    if distance <= heavy_atom_cutoff:
                        observations.append(
                            IntramolecularHydrogenBond(
                                left_key.label(), left_index, donor.name,
                                right_key.label(), right_index, acceptor.name,
                                distance,
                                "heavy_atom_distance_inference_no_hydrogens",
                                classify_atom_pair(donor.name, acceptor.name),
                                None if donor.element.upper() != acceptor.element.upper() else "same-element donor/acceptor assignment is ambiguous",
                            )
                        )
    tokens = sorted({item.fingerprint_token() for item in observations})
    return {
        "method": {
            "name": "cyclicpepomega_imhb",
            "version": IMHB_METHOD_VERSION,
            "heavy_atom_cutoff_angstrom": heavy_atom_cutoff,
            "hydrogen_policy": "hydrogens are not added; observations are heavy-atom geometric inferences unless hydrogens are present in deposited coordinates",
        },
        "intramolecular_hydrogen_bonds": [asdict(item) for item in observations],
        "fingerprint": ";".join(tokens),
        "count": len(observations),
        "warnings": [
            "protonation states are not assigned by CyclicPepOmega",
            "donor/acceptor roles for noncanonical atoms are conservative element-based candidates",
        ],
    }
