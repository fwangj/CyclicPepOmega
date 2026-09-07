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
    label_name: str | None = None
    auth_name: str | None = None
    occupancy: float | None = None


@dataclass(frozen=True)
class ComponentAtom:
    atom_id: str
    element: str
    stereo: str | None = None
    aromatic: bool = False


@dataclass(frozen=True)
class ComponentBond:
    atom_id_1: str
    atom_id_2: str
    order: str
    aromatic: bool = False
    stereo: str | None = None


@dataclass
class ComponentDefinition:
    component_id: str
    component_type: str | None = None
    name: str | None = None
    formula: str | None = None
    atoms: dict[str, ComponentAtom] = field(default_factory=dict)
    bonds: list[ComponentBond] = field(default_factory=list)
    source: str = "embedded_pdbx_chem_comp"


@dataclass
class EntityDefinition:
    entity_id: str
    entity_type: str | None = None
    polymer_type: str | None = None
    description: str | None = None
    sequence_components: list[str] = field(default_factory=list)


@dataclass
class Residue:
    key: ResidueKey
    name: str
    entity_id: str | None = None
    atoms: dict[str, Atom] = field(default_factory=dict)
    label_asym_id: str | None = None
    label_seq_id: str | None = None
    auth_asym_id: str | None = None
    auth_seq_id: str | None = None

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
    provenance_id: str | None = None
    evidence_status: str = "observed"


@dataclass
class StructureData:
    structure_id: str
    residues: dict[ResidueKey, Residue]
    explicit_bonds: list[CovalentBond] = field(default_factory=list)
    metadata: dict[str, object] = field(default_factory=dict)
    components: dict[str, ComponentDefinition] = field(default_factory=dict)
    entities: dict[str, EntityDefinition] = field(default_factory=dict)


def atom_distance(a: Atom, b: Atom) -> float:
    return dist(a.xyz, b.xyz)


def bond_dict(bond: CovalentBond) -> dict[str, object]:
    value = asdict(bond)
    value["left"] = bond.left.label()
    value["right"] = bond.right.label()
    return value
