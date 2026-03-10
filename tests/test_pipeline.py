#!/usr/bin/env python3
"""Automated test suite for the Vina + Verdict pipeline.

Runs against smoke_test/ data with synthetic fixtures.
No external dependencies beyond the standard library + project packages.

Usage:
    python tests/test_pipeline.py              # run all tests
    python tests/test_pipeline.py -v           # verbose
    python tests/test_pipeline.py TestBootstrap # run one test class
"""
import csv
import math
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.residue_utils import normalize_residue_id, parse_residue_set


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _csv_fields(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        _ = next(reader, None)
        return list(reader.fieldnames) if reader.fieldnames else []


def _write_csv(path: Path, rows: list, fieldnames: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


SMOKE_CONFIG = PROJECT_ROOT / "smoke_test" / "config.yaml"
SMOKE_OUTPUT = PROJECT_ROOT / "smoke_test" / "output" / "smoke_vina"


# ---------------------------------------------------------------------------
# Fixtures: synthetic pose table with contact residues for richer testing
# ---------------------------------------------------------------------------

FIXTURE_POSE_FIELDS = [
    "receptor_id", "ligand_id", "pose_rank", "affinity",
    "rmsd_lb", "rmsd_ub", "centroid_x", "centroid_y", "centroid_z",
    "raw_pose_file", "pocket_id", "contact_residues", "n_contact_residues",
]

def _make_fixture_poses() -> list:
    """Create synthetic poses with contact residues for full pipeline testing."""
    rows = []
    # Receptor 1: two tight poses (same pocket) + one outlier
    for lig, aff, cx, cy, cz, contacts in [
        ("lig_a", -8.0, 10.0, 10.0, 10.0, "ALA744;GLY752;LEU831"),
        ("lig_a", -7.5, 10.5, 10.2, 10.1, "ALA744;GLY752"),
        ("lig_b", -7.8, 10.3, 10.1, 10.0, "ALA744;GLY752;LEU831;ASP835"),
        ("lig_b", -6.0, 30.0, 30.0, 30.0, "PHE899;TRP950"),
    ]:
        rows.append({
            "receptor_id": "R_raw", "ligand_id": lig,
            "pose_rank": len([r for r in rows if r["receptor_id"] == "R_raw" and r["ligand_id"] == lig]) + 1,
            "affinity": aff, "rmsd_lb": 0, "rmsd_ub": 0,
            "centroid_x": cx, "centroid_y": cy, "centroid_z": cz,
            "raw_pose_file": "dummy.pdbqt", "pocket_id": "", "contact_residues": contacts,
            "n_contact_residues": len(contacts.split(";")),
        })
    # Receptor 2: similar pocket to R_raw P001
    for lig, aff, cx, cy, cz, contacts in [
        ("lig_a", -8.2, 10.2, 10.3, 10.1, "ALA744;GLY752;LEU831"),
        ("lig_b", -7.0, 10.8, 10.5, 10.3, "ALA744;GLY752"),
        ("lig_b", -6.5, 40.0, 40.0, 40.0, "ILE960;VAL970"),
    ]:
        rows.append({
            "receptor_id": "R_md1", "ligand_id": lig,
            "pose_rank": len([r for r in rows if r["receptor_id"] == "R_md1" and r["ligand_id"] == lig]) + 1,
            "affinity": aff, "rmsd_lb": 0, "rmsd_ub": 0,
            "centroid_x": cx, "centroid_y": cy, "centroid_z": cz,
            "raw_pose_file": "dummy.pdbqt", "pocket_id": "", "contact_residues": contacts,
            "n_contact_residues": len(contacts.split(";")),
        })
    return rows


def _make_fixture_config(tmpdir: Path, experimental: bool = True) -> Path:
    """Write a minimal JSON config pointing to tmpdir."""
    import json
    cfg_dict = {
        "project_name": "test_run",
        "output_root": str(tmpdir),
        "max_workers": 1,
        "mode": "blind",
        "receptors": [
            {"id": "R_raw", "source_type": "raw", "pdb": "dummy.pdb",
             "pdbqt": "dummy.pdbqt", "chain": "A"},
            {"id": "R_md1", "source_type": "md", "pdb": "dummy.pdb",
             "pdbqt": "dummy.pdbqt", "chain": "A"},
        ],
        "ligands": [
            {"id": "lig_a", "pdbqt": "dummy.pdbqt"},
            {"id": "lig_b", "pdbqt": "dummy.pdbqt"},
        ],
        "postprocess": {
            "parse_results": False,
            "extract_contacts": False,
            "contact_cutoff": 4.0,
            "cluster_pockets": False,
            "pocket_cutoff": 4.0,
            "merge_by_residue": False,
            "merge_jaccard": 0.3,
            "merge_overlap": 0.5,
            "keep_chain": False,
            "summarize_pockets": False,
            "compare_pockets": False,
            "comparison_centroid_cutoff": 15.0,
            "extract_ppi_residues": False,
            "generate_report": False,
        },
        "bootstrap": {
            "n_replicates": 10,
            "sample_fraction": 0.8,
            "seed": 42,
        },
    }
    if experimental:
        cfg_dict["experimental"] = {
            "known_binding_residues": [744, 752, 831, 835],
            "known_non_binding_residues": [899, 950, 960, 970],
            "source": "test_fixture",
        }
    cfg = tmpdir / "config.json"
    cfg.write_text(json.dumps(cfg_dict, indent=2), encoding="utf-8")
    return cfg


# ===================================================================
# Test 1: Summarize — uncertainty columns
# ===================================================================

class TestSummarize(unittest.TestCase):
    """B1: centroid_spread_A, affinity_std, affinity_iqr in pocket table."""

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp(prefix="test_summarize_"))
        self.outdir = self.tmpdir / "test_run"
        self.outdir.mkdir()
        self.cfg = _make_fixture_config(self.tmpdir)
        # Write pose table
        _write_csv(self.outdir / "vina_pose_table.csv",
                   _make_fixture_poses(), FIXTURE_POSE_FIELDS)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_uncertainty_columns_present(self):
        from egfr_pipeline.vina.summarize import summarize_from_config
        # First: cluster (assign pocket_id)
        from egfr_pipeline.vina.cluster import cluster_pose_table
        cluster_pose_table(str(self.cfg))

        pocket_csv, drug_csv = summarize_from_config(str(self.cfg))
        self.assertTrue(pocket_csv.exists())
        fields = _csv_fields(pocket_csv)
        for col in ["centroid_spread_A", "affinity_std", "affinity_iqr"]:
            self.assertIn(col, fields, f"Missing column: {col}")

    def test_spread_values_reasonable(self):
        from egfr_pipeline.vina.cluster import cluster_pose_table
        from egfr_pipeline.vina.summarize import summarize_from_config
        cluster_pose_table(str(self.cfg))
        pocket_csv, _ = summarize_from_config(str(self.cfg))
        rows = _read_csv(pocket_csv)
        self.assertGreater(len(rows), 0, "No pocket rows")
        for r in rows:
            spread = float(r["centroid_spread_A"])
            self.assertGreaterEqual(spread, 0, "Negative spread")
            # Single-pose pocket may have spread=0
            if int(r["n_pose"]) >= 2:
                aff_std = float(r["affinity_std"])
                self.assertGreaterEqual(aff_std, 0)


# ===================================================================
# Test 2: Bootstrap
# ===================================================================

class TestBootstrap(unittest.TestCase):
    """A1-A5: Bootstrap resampling pipeline."""

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp(prefix="test_bootstrap_"))
        self.outdir = self.tmpdir / "test_run"
        self.outdir.mkdir()
        self.cfg = _make_fixture_config(self.tmpdir)
        # Cluster first to get pocket_ids
        from egfr_pipeline.vina.cluster import cluster_pose_table
        _write_csv(self.outdir / "vina_pose_table.csv",
                   _make_fixture_poses(), FIXTURE_POSE_FIELDS)
        cluster_pose_table(str(self.cfg))

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_bootstrap_from_config(self):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        out = bootstrap_from_config(str(self.cfg), n_replicates=5)
        self.assertTrue(out.exists())
        rows = _read_csv(out)
        self.assertGreater(len(rows), 0, "No bootstrap rows")

    def test_bootstrap_fields(self):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config, BOOTSTRAP_FIELDS
        out = bootstrap_from_config(str(self.cfg), n_replicates=5)
        fields = _csv_fields(out)
        self.assertEqual(fields, BOOTSTRAP_FIELDS)

    def test_pocket_exists_frac_range(self):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        out = bootstrap_from_config(str(self.cfg), n_replicates=10)
        rows = _read_csv(out)
        for r in rows:
            frac = float(r["pocket_exists_frac"])
            self.assertGreaterEqual(frac, 0.0)
            self.assertLessEqual(frac, 1.0)

    def test_bootstrap_reproducible(self):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        out1 = bootstrap_from_config(str(self.cfg), n_replicates=5, seed=99,
                                     output_path=str(self.outdir / "bs1.csv"))
        out2 = bootstrap_from_config(str(self.cfg), n_replicates=5, seed=99,
                                     output_path=str(self.outdir / "bs2.csv"))
        rows1 = _read_csv(out1)
        rows2 = _read_csv(out2)
        self.assertEqual(len(rows1), len(rows2))
        for r1, r2 in zip(rows1, rows2):
            self.assertEqual(r1["pocket_exists_frac"], r2["pocket_exists_frac"])


# ===================================================================
# Test 3: Experimental Priors
# ===================================================================

class TestExperimentalPriors(unittest.TestCase):
    """C2-C5: Experimental correlation + verdict integration."""

    def test_parse_residue_ranges(self):
        from egfr_pipeline.verdict import _parse_residue_ranges
        self.assertEqual(_parse_residue_ranges("744, 752"), {744, 752})
        self.assertEqual(_parse_residue_ranges("831-835"), {831, 832, 833, 834, 835})
        self.assertEqual(_parse_residue_ranges("744, 831-833"), {744, 831, 832, 833})
        self.assertEqual(_parse_residue_ranges(""), set())
        self.assertEqual(_parse_residue_ranges("abc"), set())

    def test_compute_experimental_correlation(self):
        from egfr_pipeline.verdict import compute_experimental_correlation
        pockets = [
            {"receptor_id": "R1", "pocket_id": "P01",
             "union_contact_residues": "ALA744;GLY752;LEU831"},
            {"receptor_id": "R1", "pocket_id": "P02",
             "union_contact_residues": "PHE899;TRP950"},
        ]
        binding = {744, 752, 831}
        non_binding = {899, 950}

        result = compute_experimental_correlation(pockets, binding, non_binding)

        # P01 hits 3/3 binding, 0 non-binding
        p01 = result[("R1", "P01")]
        self.assertAlmostEqual(p01["exp_sensitivity"], 1.0, places=2)
        self.assertEqual(p01["exp_false_pos"], 0)
        self.assertGreater(p01["exp_enrichment"], 1.0)

        # P02 hits 0 binding, 2 non-binding
        p02 = result[("R1", "P02")]
        self.assertAlmostEqual(p02["exp_sensitivity"], 0.0, places=2)
        self.assertEqual(p02["exp_false_pos"], 2)

    def test_empty_experimental_data(self):
        from egfr_pipeline.verdict import compute_experimental_correlation
        result = compute_experimental_correlation([], set(), set())
        self.assertEqual(result, {})

    def test_score_pocket_exp_tags(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {
            "best_affinity": "-7.5", "n_pose": "3", "n_ligand": "2",
            "union_contact_residues": "ALA744;GLY752",
        }
        exp = {"exp_sensitivity": 0.5, "exp_hit_count": 2,
               "exp_false_pos": 0, "exp_enrichment": 3.0}
        total, verdict, reasons, v, p, c = score_pocket(
            pocket, None, [], DEFAULT_THRESHOLDS, False, exp_correlation=exp,
        )
        # Exp tags are informational — should appear in reasons
        reason_text = "; ".join(reasons)
        self.assertIn("exp_hit=", reason_text)
        self.assertIn("exp_enriched", reason_text)
        self.assertNotIn("exp_fp", reason_text)  # no false positives

    def test_score_pocket_without_exp(self):
        """Existing behavior: no exp_correlation → no exp tags."""
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {
            "best_affinity": "-7.5", "n_pose": "3", "n_ligand": "2",
            "union_contact_residues": "ALA744",
        }
        total, verdict, reasons, v, p, c = score_pocket(
            pocket, None, [], DEFAULT_THRESHOLDS, False,
        )
        reason_text = "; ".join(reasons)
        self.assertNotIn("exp_", reason_text)

    def test_rank_impact_labels(self):
        from egfr_pipeline.verdict import _compute_exp_rank_impact
        self.assertEqual(_compute_exp_rank_impact(
            {"exp_sensitivity": 0.4, "exp_enrichment": 2.5, "exp_false_pos": 0}
        ), "supports")
        self.assertEqual(_compute_exp_rank_impact(
            {"exp_sensitivity": 0.2, "exp_enrichment": 1.8, "exp_false_pos": 0}
        ), "consistent")
        self.assertEqual(_compute_exp_rank_impact(
            {"exp_sensitivity": 0, "exp_enrichment": 0, "exp_false_pos": 3}
        ), "contradicts")
        self.assertEqual(_compute_exp_rank_impact(
            {"exp_sensitivity": 0, "exp_enrichment": 0, "exp_false_pos": 0}
        ), "neutral")
        self.assertEqual(_compute_exp_rank_impact(None), "")


# ===================================================================
# Test 4: Compare — bootstrap CI column
# ===================================================================

class TestCompareBootstrapCI(unittest.TestCase):
    """B5: centroid_dist_bootstrap_ci in comparison table."""

    def test_ci_column_present(self):
        from egfr_pipeline.vina.compare import COMPARISON_FIELDS
        self.assertIn("centroid_dist_bootstrap_ci", COMPARISON_FIELDS)

    def test_ci_computation(self):
        from egfr_pipeline.vina.compare import compare_all_pockets
        pockets = [
            {"receptor_id": "R1", "pocket_id": "P01",
             "centroid_x": "10", "centroid_y": "10", "centroid_z": "10",
             "union_contact_residues": "ALA744", "best_affinity": "-7",
             "n_pose": "3", "n_ligand": "1"},
            {"receptor_id": "R2", "pocket_id": "P01",
             "centroid_x": "12", "centroid_y": "10", "centroid_z": "10",
             "union_contact_residues": "ALA744", "best_affinity": "-7.5",
             "n_pose": "4", "n_ligand": "1"},
        ]
        bootstrap_idx = {
            ("R1", "P01"): {"centroid_std_A": "0.5"},
            ("R2", "P01"): {"centroid_std_A": "0.8"},
        }
        results = compare_all_pockets(pockets, [], centroid_cutoff=0,
                                      bootstrap_index=bootstrap_idx)
        self.assertEqual(len(results), 1)
        ci = results[0]["centroid_dist_bootstrap_ci"]
        self.assertNotEqual(ci, "")
        lo, hi = ci.split("-")
        self.assertLess(float(lo), 2.0)
        self.assertGreater(float(hi), 2.0)

    def test_ci_empty_without_bootstrap(self):
        from egfr_pipeline.vina.compare import compare_all_pockets
        pockets = [
            {"receptor_id": "R1", "pocket_id": "P01",
             "centroid_x": "10", "centroid_y": "10", "centroid_z": "10",
             "union_contact_residues": "", "best_affinity": "-7",
             "n_pose": "2", "n_ligand": "1"},
            {"receptor_id": "R2", "pocket_id": "P01",
             "centroid_x": "12", "centroid_y": "10", "centroid_z": "10",
             "union_contact_residues": "", "best_affinity": "-7",
             "n_pose": "2", "n_ligand": "1"},
        ]
        results = compare_all_pockets(pockets, [], centroid_cutoff=0)
        self.assertEqual(results[0]["centroid_dist_bootstrap_ci"], "")


# ===================================================================
# Test 5: Verdict — pocket_stability from bootstrap
# ===================================================================

class TestVerdictStability(unittest.TestCase):
    """B3: pocket_stability column in valid_sites.csv."""

    def test_verdict_fields_include_stability(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        self.assertIn("pocket_stability", VERDICT_FIELDS)

    def test_verdict_fields_include_exp(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        for f in ["exp_sensitivity", "exp_specificity", "exp_enrichment", "exp_rank_impact"]:
            self.assertIn(f, VERDICT_FIELDS)


# ===================================================================
# Test 6: Schema consistency (validate.py)
# ===================================================================

class TestSchemaConsistency(unittest.TestCase):
    """C7/A7: All field lists match validate.py EXPECTED_SCHEMAS."""

    def test_verdict_fields_match(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        self.assertEqual(VERDICT_FIELDS, EXPECTED_SCHEMAS["valid_sites.csv"])

    def test_comparison_fields_match(self):
        from egfr_pipeline.vina.compare import COMPARISON_FIELDS
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        self.assertEqual(COMPARISON_FIELDS, EXPECTED_SCHEMAS["vina_pocket_comparison.csv"])

    def test_bootstrap_fields_match(self):
        from egfr_pipeline.vina.bootstrap import BOOTSTRAP_FIELDS
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        self.assertEqual(BOOTSTRAP_FIELDS, EXPECTED_SCHEMAS["vina_pocket_bootstrap.csv"])

    def test_pocket_table_fields_match(self):
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        expected = EXPECTED_SCHEMAS["vina_pocket_table.csv"]
        # Check the field list used in summarize_from_config via source
        for col in ["centroid_spread_A", "affinity_std", "affinity_iqr"]:
            self.assertIn(col, expected)

    def test_agreement_fields_match(self):
        from egfr_pipeline.verdict import AGREEMENT_FIELDS
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        self.assertEqual(AGREEMENT_FIELDS, EXPECTED_SCHEMAS["cross_method_agreement.csv"])

    def test_consensus_fields_match(self):
        from egfr_pipeline.verdict import CONSENSUS_FIELDS
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        self.assertEqual(CONSENSUS_FIELDS, EXPECTED_SCHEMAS["vina_consensus_sites.csv"])


# ===================================================================
# Test 7: End-to-end integration (cluster → summarize → bootstrap → verdict)
# ===================================================================

class TestEndToEnd(unittest.TestCase):
    """Full pipeline on synthetic fixtures."""

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp(prefix="test_e2e_"))
        self.outdir = self.tmpdir / "test_run"
        self.outdir.mkdir()
        self.cfg = _make_fixture_config(self.tmpdir)
        _write_csv(self.outdir / "vina_pose_table.csv",
                   _make_fixture_poses(), FIXTURE_POSE_FIELDS)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_full_pipeline(self):
        from egfr_pipeline.vina.cluster import cluster_pose_table
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        from egfr_pipeline.verdict import generate_verdict

        # Step 1: Cluster
        cluster_pose_table(str(self.cfg))
        pose_rows = _read_csv(self.outdir / "vina_pose_table.csv")
        pockets_assigned = sum(1 for r in pose_rows if r.get("pocket_id"))
        self.assertGreater(pockets_assigned, 0, "No pockets assigned")

        # Step 2: Summarize
        pocket_csv, drug_csv = summarize_from_config(str(self.cfg))
        pocket_rows = _read_csv(pocket_csv)
        self.assertGreater(len(pocket_rows), 0, "No pocket summary rows")

        # Step 3: Compare
        compare_csv = compare_from_config(str(self.cfg))
        self.assertTrue(compare_csv.exists())

        # Step 4: Bootstrap
        bootstrap_csv = bootstrap_from_config(str(self.cfg), n_replicates=5)
        bootstrap_rows = _read_csv(bootstrap_csv)
        self.assertGreater(len(bootstrap_rows), 0, "No bootstrap rows")

        # Step 5: Verdict (with experimental priors)
        agreement_csv, verdict_csv = generate_verdict(str(self.cfg))
        self.assertTrue(verdict_csv.exists())
        verdict_rows = _read_csv(verdict_csv)
        self.assertGreater(len(verdict_rows), 0, "No verdict rows")

        # Verify new columns exist
        fields = _csv_fields(verdict_csv)
        for col in ["exp_sensitivity", "exp_specificity", "exp_enrichment",
                     "exp_rank_impact", "pocket_stability"]:
            self.assertIn(col, fields, f"Missing in valid_sites: {col}")

        # Verify experimental data was populated (config has known residues)
        has_exp = any(r.get("exp_sensitivity") not in ("", None) for r in verdict_rows)
        self.assertTrue(has_exp, "Experimental priors not populated")

        # Verify pocket_stability was populated (bootstrap ran)
        has_stability = any(r.get("pocket_stability") not in ("", None) for r in verdict_rows)
        self.assertTrue(has_stability, "pocket_stability not populated from bootstrap")

        # Verify verdicts are valid labels
        for r in verdict_rows:
            self.assertIn(r["verdict"], ("STRONG", "MODERATE", "WEAK"))

    def test_verdict_without_experimental(self):
        """Verdict works when experimental section is empty."""
        from egfr_pipeline.vina.cluster import cluster_pose_table
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict

        # Recreate config without experimental
        self.cfg = _make_fixture_config(self.tmpdir, experimental=False)
        _write_csv(self.outdir / "vina_pose_table.csv",
                   _make_fixture_poses(), FIXTURE_POSE_FIELDS)

        cluster_pose_table(str(self.cfg))
        summarize_from_config(str(self.cfg))
        compare_from_config(str(self.cfg))
        agreement_csv, verdict_csv = generate_verdict(str(self.cfg))

        verdict_rows = _read_csv(verdict_csv)
        self.assertGreater(len(verdict_rows), 0)
        # exp fields should be empty
        for r in verdict_rows:
            self.assertEqual(r.get("exp_sensitivity", ""), "")

    def test_verdict_without_bootstrap(self):
        """Verdict works when bootstrap CSV doesn't exist."""
        from egfr_pipeline.vina.cluster import cluster_pose_table
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict

        cluster_pose_table(str(self.cfg))
        summarize_from_config(str(self.cfg))
        compare_from_config(str(self.cfg))
        # Skip bootstrap intentionally
        agreement_csv, verdict_csv = generate_verdict(str(self.cfg))

        verdict_rows = _read_csv(verdict_csv)
        self.assertGreater(len(verdict_rows), 0)
        # stability should be empty
        for r in verdict_rows:
            self.assertEqual(r.get("pocket_stability", ""), "")


# ===================================================================
# Test 8: Regression — existing smoke output still valid
# ===================================================================

class TestSmokeRegression(unittest.TestCase):
    """Verify smoke_test/output files haven't broken."""

    @unittest.skipUnless(SMOKE_OUTPUT.exists(), "smoke_test/output not found")
    def test_smoke_pose_table_readable(self):
        rows = _read_csv(SMOKE_OUTPUT / "vina_pose_table.csv")
        self.assertEqual(len(rows), 12, "Expected 12 smoke poses")

    @unittest.skipUnless(SMOKE_OUTPUT.exists(), "smoke_test/output not found")
    def test_smoke_pocket_table_has_uncertainty_header(self):
        fields = _csv_fields(SMOKE_OUTPUT / "vina_pocket_table.csv")
        for col in ["centroid_spread_A", "affinity_std", "affinity_iqr"]:
            self.assertIn(col, fields)

    @unittest.skipUnless(
        (SMOKE_OUTPUT / "valid_sites.csv").exists(),
        "valid_sites.csv not in smoke output"
    )
    def test_smoke_verdict_verdicts(self):
        rows = _read_csv(SMOKE_OUTPUT / "valid_sites.csv")
        for r in rows:
            self.assertIn(r["verdict"], ("STRONG", "MODERATE", "WEAK"))


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    unittest.main(verbosity=2)
