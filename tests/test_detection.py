import unittest

from cyclicpepomega.detect import analyze_structure
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
        keys = list(residues)
        bonds = [
            CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide"),
            CovalentBond(keys[1], "C", keys[2], "N", "struct_conn", "peptide"),
            CovalentBond(keys[0], "SG", keys[2], "SG", "struct_conn", "disulfide"),
        ]
        report = analyze_structure(StructureData("TEST", residues, bonds))
        self.assertEqual(report["cyclic_peptides"][0]["cyclization_type"], "disulfide-constrained")

    def test_linear_candidate_is_reported(self):
        residues = {r.key: r for r in (residue(1, 0), residue(2, 4))}
        keys = list(residues)
        bonds = [CovalentBond(keys[0], "C", keys[1], "N", "struct_conn", "peptide")]
        report = analyze_structure(StructureData("TEST", residues, bonds))
        self.assertEqual(report["summary"]["cyclic_peptide_count"], 0)
        self.assertEqual(len(report["ambiguous_or_noncyclic_candidates"]), 1)


if __name__ == "__main__":
    unittest.main()
