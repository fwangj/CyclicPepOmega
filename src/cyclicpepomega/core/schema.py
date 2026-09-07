from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import dist


@dataclass(frozen=True, order=True)
class ResidueKey:
    model: int
    chain: str
    seq: str
    insertion: str = ""

    def label(self) -> str:
        return f"{self.chain}:{self.seq}{self.insertion}"


@dataclass
class Atom:
    name: str
    element: str
    xyz: tuple[float, float, float]
    altloc: str = ""


@dataclass
class Residue:
    key: ResidueKey
    name: str
    entity_id: str | None = None
    atoms: dict[str, Atom] = field(default_factory=dict)

    def atom(self, name: str) -> Atom | None:
        return self.atoms.get(name)

    @property
    def peptide_like(self) -> bool:
        # Chemistry-first proxy: an alpha-amino-acid-like backbone, independent
        # of residue/component names.
        return all(name in self.atoms for name in ("N", "CA", "C"))


@dataclass(frozen=True)
class CovalentBond:
    left: ResidueKey
    left_atom: str
    right: ResidueKey
    right_atom: str
    source: str
    kind: str = "covalent"


@dataclass
class StructureData:
    structure_id: str
    residues: dict[ResidueKey, Residue]
    explicit_bonds: list[CovalentBond] = field(default_factory=list)
    metadata: dict[str, object] = field(default_factory=dict)


def atom_distance(a: Atom, b: Atom) -> float:
    return dist(a.xyz, b.xyz)


def bond_dict(bond: CovalentBond) -> dict[str, object]:
    value = asdict(bond)
    value["left"] = bond.left.label()
    value["right"] = bond.right.label()
    return value

