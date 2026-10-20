"""Lossless migration tests: global substitutions must corrupt these fixtures."""
import unittest

from cationize_protein.rename import rename_parm7, rename_pdb

MAPPING = {"ASL": "ASC", "GLS": "GLC", "NSL": "ASX", "LCL": "LET"}

PDB = (
    b"REMARK ASL GLS NSL LCL must remain verbatim\r\n"
    b"ATOM      1  CA  ASL A   1      11.123  22.456  -3.789  1.00 20.00           C  \r\n"
    b"HETATM    2  C1  GLS B   7      12.123  23.456  -4.789  1.00 20.00           C  \n"
    b"TER       3      ASL A   1\n"
    b"END\n"
)
PARM7 = (
    b"%VERSION VERSION_STAMP = V0001.000\r\n"
    b"%FLAG TITLE\r\n%FORMAT(20a4)\r\nASL GLS NSL LCL title retained\r\n"
    b"%FLAG ATOM_NAME\n%FORMAT(20a4)\nASL C1  \n"
    b"%FLAG CHARGE\n%FORMAT(5E16.8)\n  1.82223000E+01 -1.82223000E+01\n"
    b"%FLAG RESIDUE_LABEL\r\n%FORMAT(20a4)\r\nASL GLS NSL LCL ALA \r\n"
    b"%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1       2       3       4       5\n"
    b"%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n       0       3       1\n"
)


class RenameBytesTests(unittest.TestCase):
    def test_pdb_only_atom_and_hetatm_resname_columns_change(self):
        expected = PDB.replace(b" CA  ASL", b" CA  ASC").replace(b" C1  GLS", b" C1  GLC")
        self.assertEqual(rename_pdb(PDB, MAPPING), expected)

    def test_parm7_only_residue_label_section_changes(self):
        expected = PARM7.replace(b"ASL GLS NSL LCL ALA ", b"ASC GLC ASX LET ALA ")
        self.assertEqual(rename_parm7(PARM7, MAPPING), expected)

    def test_both_formats_are_exactly_reversible(self):
        inverse = {new: old for old, new in MAPPING.items()}
        for data, rename in [(PDB, rename_pdb), (PARM7, rename_parm7)]:
            with self.subTest(format=rename.__name__):
                self.assertEqual(rename(rename(data, MAPPING), inverse), data)

    def test_existing_target_label_is_rejected_as_irreversible(self):
        for data, rename in [(PDB.replace(b" C1  GLS", b" C1  ASC"), rename_pdb),
                             (PARM7.replace(b"LCL ALA ", b"LCL ASC "), rename_parm7)]:
            with self.subTest(format=rename.__name__):
                with self.assertRaisesRegex(ValueError, "collision"):
                    rename(data, MAPPING)

    def test_noninjective_mapping_is_rejected(self):
        with self.assertRaises(ValueError):
            rename_pdb(PDB, {"ASL": "ASC", "GLS": "ASC"})

    def test_width_changing_or_invalid_names_are_rejected(self):
        for mapping in ({"ASL": "A"}, {"ASL": "A C"}, {"ASL": "ABCDE"}):
            with self.subTest(mapping=mapping):
                with self.assertRaises(ValueError):
                    rename_pdb(PDB, mapping)

    def test_missing_duplicate_or_malformed_residue_section_is_rejected(self):
        for data in (PARM7.replace(b"%FLAG RESIDUE_LABEL", b"%FLAG OTHER"),
                     PARM7 + b"%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nASL \n",
                     PARM7.replace(b"%FORMAT(20a4)\r\nASL GLS NSL", b"%FORMAT(10I8)\r\nASL GLS NSL")):
            with self.subTest(data=data):
                with self.assertRaises(ValueError):
                    rename_parm7(data, MAPPING)

    def test_generic_mapping_can_exchange_existing_names(self):
        self.assertEqual(rename_pdb(PDB, {"ASL": "GLS", "GLS": "ASL"}),
                         PDB.replace(b" CA  ASL", b" CA  GLS").replace(b" C1  GLS", b" C1  ASL"))


class ExportSeriesTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        from pathlib import Path
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "source"
        self.level = self.source / "25"
        self.level.mkdir(parents=True)
        (self.level / "protein.parm7").write_bytes(PARM7)
        (self.level / "minimized.pdb").write_bytes(PDB)
        (self.level / "minimized.rst7").write_bytes(b"Original coordinates\n2\n  1.0  2.0  3.0\n")
        self.output = self.source / "structures_named"

    def test_copies_are_reversible_and_sources_restarts_are_preserved(self):
        import hashlib
        import json
        from cationize_protein.rename import export_series
        original = {p.name: p.read_bytes() for p in self.level.iterdir()}
        report = export_series(self.source, self.output, MAPPING, levels=("25",))
        self.assertTrue((self.output / "25" / "protein.parm7").exists())
        self.assertEqual((self.output / "25" / "minimized.pdb").read_bytes(),
                         PDB.replace(b" CA  ASL", b" CA  ASC").replace(b" C1  GLS", b" C1  GLC"))
        self.assertTrue((self.output / "25" / "minimized.rst7").is_symlink())
        self.assertEqual((self.output / "25" / "minimized.rst7").resolve(),
                         (self.level / "minimized.rst7").resolve())
        self.assertEqual({p.name: p.read_bytes() for p in self.level.iterdir()}, original)
        self.assertTrue(report["original_inputs_unchanged"])
        record = report["cases"]["25"]["files"]["protein.parm7"]
        self.assertEqual(record["source_sha256"], hashlib.sha256(PARM7).hexdigest())
        self.assertTrue(record["validation"]["exact_inverse"])
        self.assertTrue(record["validation"]["all_other_bytes_identical"])
        self.assertEqual(json.loads((self.output / "migration_manifest.json").read_text()), report)
        self.assertEqual(report["cases"]["25"]["coordinates"]["sha256"],
                         hashlib.sha256(original["minimized.rst7"]).hexdigest())

    def test_optional_solvated_topology_is_exported_when_present(self):
        from cationize_protein.rename import export_series
        (self.level / "system.parm7").write_bytes(PARM7)
        export_series(self.source, self.output, MAPPING, levels=("25",))
        self.assertTrue((self.output / "25" / "system.parm7").exists())
        self.assertEqual((self.output / "25" / "system.parm7").read_bytes(),
                         PARM7.replace(b"ASL GLS NSL LCL ALA ", b"ASC GLC ASX LET ALA "))

    def test_existing_output_is_never_overwritten(self):
        from cationize_protein.rename import export_series
        self.output.mkdir()
        marker = self.output / "keep"
        marker.write_bytes(b"unchanged")
        with self.assertRaises(FileExistsError):
            export_series(self.source, self.output, MAPPING, levels=("25",))
        self.assertEqual(marker.read_bytes(), b"unchanged")

    def test_missing_source_leaves_no_partial_exports_or_temporary_files(self):
        from cationize_protein.rename import export_series
        (self.level / "minimized.pdb").unlink()
        with self.assertRaises(FileNotFoundError):
            export_series(self.source, self.output, MAPPING, levels=("25",))
        self.assertFalse(self.output.exists())
        self.assertEqual(sorted(p.name for p in self.source.iterdir()), ["25"])

    def test_traversal_and_duplicate_levels_are_rejected(self):
        from cationize_protein.rename import export_series
        for levels in (("../elsewhere",), ("25", "25"), ()):
            with self.subTest(levels=levels):
                with self.assertRaises(ValueError):
                    export_series(self.source, self.output, MAPPING, levels=levels)

    def test_broken_symlink_output_is_refused_without_creating_its_target(self):
        from cationize_protein.rename import export_series
        dangling = self.root / "must_not_be_created"
        self.output.symlink_to(dangling, target_is_directory=True)
        with self.assertRaises(FileExistsError):
            export_series(self.source, self.output, MAPPING, levels=("25",))
        self.assertFalse(dangling.exists())
        self.assertTrue(self.output.is_symlink())


if __name__ == "__main__":
    unittest.main()
