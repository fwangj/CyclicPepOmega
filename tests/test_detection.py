import unittest

from cyclicpepomega.detect import DetectionConfig, analyze_structure, reconstruct_bonds
from cyclicpepomega.core.schema import Atom, CovalentBond, Residue, ResidueKey, StructureData


def residue(number, x=0.0, sulfur=False):
    key = ResidueKey(1, "A", str(number))
    atoms = {
        "N": Atom("N", "N", (x, 0.0, 0.0)),
        "CA": Atom("CA", "C", (x + 0.5, 1.0, 0.0)),
        "C": Atom("C", "C", (x + 1.0, 0.0, 0.0)),
    }
    if sulfur:
        atoms["SG"] = Atom("SG", "S", (x, 2.0, 0.0))
    return Residue(key, "CYS" if sulfur else "UNK", "1", atoms)


class DetectionTests(unittest.TestCase):
    def test_polymer_connectivity_does_not_join_entity_copies(self):
        left = residue(1, 0)
        right = residue(2, 4)
        copy_left = residue(1, 20)
        copy_right = residue(2, 24)
        for item, chain in ((left, "A"), (right, "A"), (copy_left, "B"), (copy_right, "B")):
            old = item.key
            item.key = ResidueKey(1, chain, old.seq)
            item.entity_id = "1"
            item.label_asym_id = chain
            item.label_seq_id = old.seq
        residues = {item.key: item for item in (left, right, copy_left, copy_right)}
        bonds = reconstruct_bonds(StructureData("COPIES", residues), DetectionConfig())
        self.assertEqual(sum(b.source == "polymer_connectivity" for b in bonds), 2)
        self.assertFalse(any(b.left.chain != b.right.chain for b in bonds))

    def test_head_to_tail_from_connectivity(self):
        residues = {r.key: r for r in (residue(1, 0), residue(2, 4), residue(3, 8))}
        keys = list(residues)
        bonds = [
            CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide"),
            CovalentBond(keys[1], "C", keys[2], "N", "struct_conn", "peptide"),
            CovalentBond(keys[2], "C", keys[0], "N", "struct_conn", "peptide"),
        ]
        report = analyze_structure(StructureData("TEST", residues, bonds))
        self.assertEqual(report["cyclic_peptides"][0]["cyclization_type"], "head-to-tail")
        self.assertEqual([r["original_id"] for r in report["cyclic_peptides"][0]["residues"]], ["A:1", "A:2", "A:3"])

    def test_disulfide_constrained(self):
        residues = {r.key: r for r in (residue(1, 0, True), residue(2, 4), residue(3, 8, True))}
        for item in residues.values():
            item.entity_id = "2"
        keys = list(residues)
        bonds = [
            CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide"),
            CovalentBond(keys[1], "C", keys[2], "N", "struct_conn", "peptide"),
            CovalentBond(keys[0], "SG", keys[2], "SG", "struct_conn", "disulfide"),
        ]
        report = analyze_structure(StructureData("1NPO", residues, bonds))
        self.assertEqual(report["cyclic_peptides"][0]["cyclization_type"], "disulfide-constrained")

    def test_unscoped_disulfide_polypeptide_is_ambiguous(self):
        residues = {r.key: r for r in (residue(1, 0, True), residue(2, 4), residue(3, 8, True))}
        keys = list(residues)
        bonds = [
            CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide"),
            CovalentBond(keys[1], "C", keys[2], "N", "struct_conn", "peptide"),
            CovalentBond(keys[0], "SG", keys[2], "SG", "struct_conn", "disulfide"),
        ]
        report = analyze_structure(StructureData("UNKNOWN", residues, bonds))
        self.assertEqual(report["summary"]["cyclic_peptide_count"], 0)
        self.assertEqual(report["ambiguous_or_noncyclic_candidates"][0]["cyclization_type"], "ambiguous_unresolved")

    def test_linear_candidate_is_reported(self):
        residues = {r.key: r for r in (residue(1, 0), residue(2, 4))}
        keys = list(residues)
        bonds = [CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide")]
        report = analyze_structure(StructureData("TEST", residues, bonds))
        self.assertEqual(report["summary"]["cyclic_peptide_count"], 0)
        self.assertEqual(len(report["ambiguous_or_noncyclic_candidates"]), 1)


if __name__ == "__main__":
    unittest.main()
