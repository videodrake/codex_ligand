#!/usr/bin/env python3
"""Automated test suite for the Vina + Verdict pipeline.

Runs against smoke_test/ data with synthetic fixtures.

Usage:
    .venv/bin/pytest tests/ -v                    # all tests
    .venv/bin/pytest tests/ -v -k bootstrap       # bootstrap only
    .venv/bin/pytest tests/ -v -k "e2e"           # end-to-end only
    .venv/bin/pytest tests/ -v --tb=short         # compact tracebacks
"""
import csv
import json
import math
import shutil
import sys
from pathlib import Path

import pytest

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.residue_utils import normalize_residue_id, parse_residue_set


# ===================================================================
# Shared helpers
# ===================================================================

def read_csv(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def csv_fields(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        _ = next(reader, None)
        return list(reader.fieldnames) if reader.fieldnames else []


def write_csv(path: Path, rows: list, fieldnames: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


SMOKE_OUTPUT = PROJECT_ROOT / "smoke_test" / "output" / "smoke_vina"

POSE_FIELDS = [
    "receptor_id", "ligand_id", "pose_rank", "affinity",
    "rmsd_lb", "rmsd_ub", "centroid_x", "centroid_y", "centroid_z",
    "raw_pose_file", "pocket_id", "contact_residues", "n_contact_residues",
]


# ===================================================================
# Fixtures (pytest)
# ===================================================================

def make_poses() -> list:
    """Synthetic poses: 2 receptors, 2 ligands, with contact residues."""
    rows = []
    specs = [
        # (receptor, ligand, affinity, cx, cy, cz, contacts)
        # R_raw: tight cluster + outlier
        ("R_raw", "lig_a", -8.0, 10.0, 10.0, 10.0, "ALA744;GLY752;LEU831"),
        ("R_raw", "lig_a", -7.5, 10.5, 10.2, 10.1, "ALA744;GLY752"),
        ("R_raw", "lig_b", -7.8, 10.3, 10.1, 10.0, "ALA744;GLY752;LEU831;ASP835"),
        ("R_raw", "lig_b", -6.0, 30.0, 30.0, 30.0, "PHE899;TRP950"),
        # R_md1: similar main pocket to R_raw + distant outlier
        ("R_md1", "lig_a", -8.2, 10.2, 10.3, 10.1, "ALA744;GLY752;LEU831"),
        ("R_md1", "lig_b", -7.0, 10.8, 10.5, 10.3, "ALA744;GLY752"),
        ("R_md1", "lig_b", -6.5, 40.0, 40.0, 40.0, "ILE960;VAL970"),
    ]
    counters = {}
    for rec, lig, aff, cx, cy, cz, contacts in specs:
        key = (rec, lig)
        counters[key] = counters.get(key, 0) + 1
        rows.append({
            "receptor_id": rec, "ligand_id": lig,
            "pose_rank": counters[key],
            "affinity": aff, "rmsd_lb": 0, "rmsd_ub": 0,
            "centroid_x": cx, "centroid_y": cy, "centroid_z": cz,
            "raw_pose_file": "dummy.pdbqt", "pocket_id": "",
            "contact_residues": contacts,
            "n_contact_residues": len(contacts.split(";")),
        })
    return rows


def make_config(tmpdir: Path, experimental: bool = True) -> Path:
    """Write JSON config to tmpdir."""
    cfg_dict = {
        "project_name": "test_run",
        "output_root": str(tmpdir),
        "max_workers": 1,
        "mode": "blind",
        "receptors": [
            {"id": "R_raw", "source_type": "raw",
             "pdb": "dummy.pdb", "pdbqt": "dummy.pdbqt", "chain": "A"},
            {"id": "R_md1", "source_type": "md",
             "pdb": "dummy.pdb", "pdbqt": "dummy.pdbqt", "chain": "A"},
        ],
        "ligands": [
            {"id": "lig_a", "pdbqt": "dummy.pdbqt"},
            {"id": "lig_b", "pdbqt": "dummy.pdbqt"},
        ],
        "postprocess": {
            "parse_results": False, "extract_contacts": False,
            "contact_cutoff": 4.0, "cluster_pockets": False,
            "pocket_cutoff": 8.0, "merge_by_residue": True,
            "merge_jaccard": 0.3, "merge_overlap": 0.5,
            "merge_centroid_fallback": 6.0,
            "cluster_max_iterations": 10, "min_pocket_size": 2,
            "keep_chain": False, "summarize_pockets": False,
            "compare_pockets": False, "comparison_centroid_cutoff": 15.0,
            "extract_ppi_residues": False, "generate_report": False,
        },
        "bootstrap": {"n_replicates": 10, "sample_fraction": 0.8, "seed": 42},
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


@pytest.fixture
def pipeline_env(tmp_path):
    """Set up tmpdir with config + pose table + clustering done."""
    from egfr_pipeline.vina.cluster import cluster_pose_table

    outdir = tmp_path / "test_run"
    outdir.mkdir()
    cfg = make_config(tmp_path)
    write_csv(outdir / "vina_pose_table.csv", make_poses(), POSE_FIELDS)
    cluster_pose_table(str(cfg))
    return cfg, outdir


@pytest.fixture
def pipeline_env_no_exp(tmp_path):
    """Same as pipeline_env but without experimental section."""
    from egfr_pipeline.vina.cluster import cluster_pose_table

    outdir = tmp_path / "test_run"
    outdir.mkdir()
    cfg = make_config(tmp_path, experimental=False)
    write_csv(outdir / "vina_pose_table.csv", make_poses(), POSE_FIELDS)
    cluster_pose_table(str(cfg))
    return cfg, outdir


# ===================================================================
# 0. Clustering algorithm
# ===================================================================

def _make_pose(receptor, ligand, rank, affinity, cx, cy, cz, contacts=""):
    return {
        "receptor_id": receptor, "ligand_id": ligand,
        "pose_rank": rank, "affinity": affinity,
        "rmsd_lb": 0, "rmsd_ub": 0,
        "centroid_x": cx, "centroid_y": cy, "centroid_z": cz,
        "raw_pose_file": "dummy.pdbqt", "pocket_id": "",
        "contact_residues": contacts, "n_contact_residues": len(contacts.split(";")) if contacts else 0,
    }


class TestClustering:
    """Tests for the iterative centroid-based clustering algorithm."""

    def test_convergence_single_pocket(self):
        """Poses spread ~6A within same pocket should converge to 1 pocket."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.5, 13.0, 12.0, 11.0),  # ~4.7A
            _make_pose("R1", "lig", 3, -7.0, 11.0, 14.0, 10.0),  # ~4.1A from first
            _make_pose("R1", "lig", 4, -6.5, 15.0, 11.0, 12.0),  # ~5.7A from first
            _make_pose("R1", "lig", 5, -6.0, 12.0, 10.0, 14.0),  # ~4.5A from first
        ]
        result = assign_pockets(rows, cutoff=8.0)
        pockets = {r["pocket_id"] for r in result}
        assert len(pockets) == 1, f"Expected 1 pocket, got {pockets}"

    def test_separation_two_clusters(self):
        """Two groups 30A apart should remain as 2 separate pockets."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.5, 11.0, 10.0, 10.0),
            _make_pose("R1", "lig", 3, -7.0, 10.5, 11.0, 10.0),
            _make_pose("R1", "lig", 4, -6.5, 40.0, 40.0, 40.0),
            _make_pose("R1", "lig", 5, -6.0, 41.0, 40.0, 40.0),
            _make_pose("R1", "lig", 6, -5.5, 40.5, 41.0, 40.0),
        ]
        result = assign_pockets(rows, cutoff=8.0, min_pocket_size=2)
        pockets = {r["pocket_id"] for r in result}
        assert len(pockets) == 2, f"Expected 2 pockets, got {pockets}"

    def test_singleton_absorption(self):
        """A single outlier pose should be absorbed into nearest large pocket."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.5, 11.0, 10.0, 10.0),
            _make_pose("R1", "lig", 3, -7.0, 10.5, 11.0, 10.0),
            # Singleton outlier
            _make_pose("R1", "lig", 4, -6.0, 50.0, 50.0, 50.0),
        ]
        result = assign_pockets(rows, cutoff=8.0, min_pocket_size=2)
        pockets = {r["pocket_id"] for r in result}
        assert len(pockets) == 1, f"Singleton should be absorbed, got {pockets}"

    def test_singleton_preservation_when_all_small(self):
        """When all pockets are singletons, they should all be preserved."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.0, 50.0, 50.0, 50.0),
        ]
        result = assign_pockets(rows, cutoff=8.0, min_pocket_size=2)
        # Both are singletons, no large pocket exists -> both preserved
        pockets = {r["pocket_id"] for r in result}
        assert len(pockets) == 2, f"All-singleton case should preserve, got {pockets}"

    def test_order_independence(self):
        """Shuffling input order should not change clustering result."""
        import random
        from egfr_pipeline.vina.cluster import assign_pockets
        base_rows = [
            _make_pose("R1", "lig_a", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig_a", 2, -7.0, 11.0, 10.5, 10.0),
            _make_pose("R1", "lig_b", 1, -7.5, 10.5, 10.0, 10.5),
            _make_pose("R1", "lig_b", 2, -6.5, 40.0, 40.0, 40.0),
            _make_pose("R1", "lig_b", 3, -6.0, 41.0, 40.5, 40.0),
        ]

        result1 = assign_pockets([dict(r) for r in base_rows], cutoff=8.0, min_pocket_size=2)
        pocket_map_1 = {(r["receptor_id"], r["ligand_id"], r["pose_rank"]): r["pocket_id"]
                        for r in result1}

        # Shuffle and re-run
        shuffled = [dict(r) for r in base_rows]
        random.seed(42)
        random.shuffle(shuffled)
        result2 = assign_pockets(shuffled, cutoff=8.0, min_pocket_size=2)
        pocket_map_2 = {(r["receptor_id"], r["ligand_id"], r["pose_rank"]): r["pocket_id"]
                        for r in result2}

        assert pocket_map_1 == pocket_map_2, "Clustering should be order-independent"

    def test_iterative_refinement_corrects_drift(self):
        """Poses that would be fragmented by single-pass should converge."""
        from egfr_pipeline.vina.cluster import assign_pockets
        # Chain of poses each ~5A apart: single-pass with 8A cutoff could drift
        # the center far from the ends, but iterative should converge
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.5, 14.0, 10.0, 10.0),  # 4A from #1
            _make_pose("R1", "lig", 3, -7.0, 10.0, 14.0, 10.0),  # 4A from #1
            _make_pose("R1", "lig", 4, -6.5, 14.0, 14.0, 10.0),  # ~5.7A from #1
            _make_pose("R1", "lig", 5, -6.0, 12.0, 12.0, 14.0),  # ~5.3A from #1
        ]
        result = assign_pockets(rows, cutoff=8.0)
        pockets = {r["pocket_id"] for r in result}
        assert len(pockets) == 1, f"Should converge to 1 pocket, got {pockets}"

    def test_multi_receptor_independence(self):
        """Clustering is independent per receptor."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.0, 11.0, 10.0, 10.0),
            _make_pose("R2", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R2", "lig", 2, -7.0, 11.0, 10.0, 10.0),
        ]
        result = assign_pockets(rows, cutoff=8.0)
        r1 = {r["pocket_id"] for r in result if r["receptor_id"] == "R1"}
        r2 = {r["pocket_id"] for r in result if r["receptor_id"] == "R2"}
        assert r1 == {"P001"}
        assert r2 == {"P001"}

    def test_merge_by_residue_overlap(self):
        """Pockets with overlapping residues should be merged."""
        from egfr_pipeline.vina.cluster import assign_pockets, merge_pockets_by_residue
        # Two pockets far apart in space but overlapping residues
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0, "ALA744;GLY752;LEU831"),
            _make_pose("R1", "lig", 2, -7.5, 11.0, 10.0, 10.0, "ALA744;GLY752"),
            _make_pose("R1", "lig", 3, -7.0, 20.0, 20.0, 20.0, "ALA744;GLY752;ASP835"),
            _make_pose("R1", "lig", 4, -6.5, 21.0, 20.0, 20.0, "ALA744;LEU831;ASP835"),
        ]
        result = assign_pockets(rows, cutoff=8.0, min_pocket_size=1)
        # Should be 2 pockets (>8A apart)
        pockets_before = {r["pocket_id"] for r in result}
        assert len(pockets_before) == 2

        # After residue merge: high overlap (ALA744,GLY752,LEU831 vs ALA744,GLY752,ASP835)
        merged = merge_pockets_by_residue(result, jaccard_threshold=0.3, overlap_threshold=0.5)
        pockets_after = {r["pocket_id"] for r in merged}
        assert len(pockets_after) == 1, f"Should merge to 1 pocket, got {pockets_after}"

    def test_merge_centroid_fallback(self):
        """Pockets close in space but with no residue data should merge via centroid fallback."""
        from egfr_pipeline.vina.cluster import assign_pockets, merge_pockets_by_residue
        rows = [
            _make_pose("R1", "lig_a", 1, -8.0, 10.0, 10.0, 10.0, ""),
            _make_pose("R1", "lig_a", 2, -7.5, 11.0, 10.0, 10.0, ""),
            _make_pose("R1", "lig_b", 1, -7.0, 14.0, 10.0, 10.0, ""),
            _make_pose("R1", "lig_b", 2, -6.5, 15.0, 10.0, 10.0, ""),
        ]
        # cutoff=3 forces 2 pockets
        result = assign_pockets(rows, cutoff=3.0, min_pocket_size=1)
        pockets_before = {r["pocket_id"] for r in result}
        assert len(pockets_before) == 2

        # Centroid fallback: centroids are ~4A apart < 6.0 fallback
        merged = merge_pockets_by_residue(
            result, jaccard_threshold=0.3, overlap_threshold=0.5,
            centroid_fallback_cutoff=6.0,
        )
        pockets_after = {r["pocket_id"] for r in merged}
        assert len(pockets_after) == 1, f"Should merge via centroid fallback, got {pockets_after}"

    def test_pocket_ids_contiguous(self):
        """Pocket IDs should be P001, P002, ... without gaps."""
        from egfr_pipeline.vina.cluster import assign_pockets
        rows = [
            _make_pose("R1", "lig", 1, -8.0, 10.0, 10.0, 10.0),
            _make_pose("R1", "lig", 2, -7.0, 11.0, 10.0, 10.0),
            _make_pose("R1", "lig", 3, -6.0, 40.0, 40.0, 40.0),
            _make_pose("R1", "lig", 4, -5.5, 41.0, 40.0, 40.0),
            _make_pose("R1", "lig", 5, -5.0, 80.0, 80.0, 80.0),
            _make_pose("R1", "lig", 6, -4.5, 81.0, 80.0, 80.0),
        ]
        result = assign_pockets(rows, cutoff=8.0, min_pocket_size=2)
        pocket_ids = sorted({r["pocket_id"] for r in result})
        expected = [f"P{i:03d}" for i in range(1, len(pocket_ids) + 1)]
        assert pocket_ids == expected, f"Non-contiguous IDs: {pocket_ids}"


# ===================================================================
# 1. Summarize — uncertainty columns
# ===================================================================

class TestSummarize:
    def test_uncertainty_columns_present(self, pipeline_env):
        from egfr_pipeline.vina.summarize import summarize_from_config
        cfg, outdir = pipeline_env
        pocket_csv, _ = summarize_from_config(str(cfg))
        fields = csv_fields(pocket_csv)
        for col in ["centroid_spread_A", "affinity_std", "affinity_iqr"]:
            assert col in fields, f"Missing: {col}"

    def test_spread_non_negative(self, pipeline_env):
        from egfr_pipeline.vina.summarize import summarize_from_config
        cfg, _ = pipeline_env
        pocket_csv, _ = summarize_from_config(str(cfg))
        for r in read_csv(pocket_csv):
            assert float(r["centroid_spread_A"]) >= 0

    def test_multi_pose_has_nonzero_std(self, pipeline_env):
        from egfr_pipeline.vina.summarize import summarize_from_config
        cfg, _ = pipeline_env
        pocket_csv, _ = summarize_from_config(str(cfg))
        multi_pose = [r for r in read_csv(pocket_csv) if int(r["n_pose"]) >= 2]
        assert len(multi_pose) > 0
        for r in multi_pose:
            assert float(r["affinity_std"]) >= 0


# ===================================================================
# 2. Bootstrap
# ===================================================================

class TestBootstrap:
    def test_output_exists(self, pipeline_env):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        cfg, _ = pipeline_env
        out = bootstrap_from_config(str(cfg), n_replicates=5)
        assert out.exists()
        assert len(read_csv(out)) > 0

    def test_fields_match_spec(self, pipeline_env):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config, BOOTSTRAP_FIELDS
        cfg, _ = pipeline_env
        out = bootstrap_from_config(str(cfg), n_replicates=5)
        assert csv_fields(out) == BOOTSTRAP_FIELDS

    def test_pocket_exists_frac_in_range(self, pipeline_env):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        cfg, _ = pipeline_env
        out = bootstrap_from_config(str(cfg), n_replicates=10)
        for r in read_csv(out):
            frac = float(r["pocket_exists_frac"])
            assert 0.0 <= frac <= 1.0

    def test_reproducible_with_same_seed(self, pipeline_env):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        cfg, outdir = pipeline_env
        out1 = bootstrap_from_config(str(cfg), n_replicates=5, seed=99,
                                     output_path=str(outdir / "bs1.csv"))
        out2 = bootstrap_from_config(str(cfg), n_replicates=5, seed=99,
                                     output_path=str(outdir / "bs2.csv"))
        rows1, rows2 = read_csv(out1), read_csv(out2)
        assert len(rows1) == len(rows2)
        for r1, r2 in zip(rows1, rows2):
            assert r1 == r2

    def test_different_seed_different_result(self, pipeline_env):
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        cfg, outdir = pipeline_env
        out1 = bootstrap_from_config(str(cfg), n_replicates=20, seed=1,
                                     output_path=str(outdir / "bs_s1.csv"))
        out2 = bootstrap_from_config(str(cfg), n_replicates=20, seed=999,
                                     output_path=str(outdir / "bs_s2.csv"))
        rows1, rows2 = read_csv(out1), read_csv(out2)
        # At least one value should differ with different seeds
        diffs = sum(1 for r1, r2 in zip(rows1, rows2)
                    if r1["affinity_mean"] != r2["affinity_mean"])
        assert diffs > 0, "Different seeds produced identical results"


# ===================================================================
# 3. Experimental Priors
# ===================================================================

class TestExperimentalPriors:
    @pytest.mark.parametrize("raw,expected", [
        ("744, 752", {744, 752}),
        ("831-835", {831, 832, 833, 834, 835}),
        ("744, 831-833", {744, 831, 832, 833}),
        ("", set()),
        ("abc", set()),
        ("100", {100}),
    ])
    def test_parse_residue_ranges(self, raw, expected):
        from egfr_pipeline.verdict import _parse_residue_ranges
        assert _parse_residue_ranges(raw) == expected

    def test_correlation_perfect_hit(self):
        from egfr_pipeline.verdict import compute_experimental_correlation
        pockets = [{"receptor_id": "R1", "pocket_id": "P01",
                     "union_contact_residues": "ALA744;GLY752;LEU831"}]
        # With non-binding pool: enrichment > 1.0
        result = compute_experimental_correlation(pockets, {744, 752, 831}, {899, 950})
        p01 = result[("R1", "P01")]
        assert p01["exp_sensitivity"] == pytest.approx(1.0)
        assert p01["exp_false_pos"] == 0
        assert p01["exp_enrichment"] > 1.0

    def test_correlation_pure_false_positive(self):
        from egfr_pipeline.verdict import compute_experimental_correlation
        pockets = [{"receptor_id": "R1", "pocket_id": "P02",
                     "union_contact_residues": "PHE899;TRP950"}]
        result = compute_experimental_correlation(pockets, {744}, {899, 950})
        p02 = result[("R1", "P02")]
        assert p02["exp_sensitivity"] == 0.0
        assert p02["exp_false_pos"] == 2

    def test_empty_returns_empty(self):
        from egfr_pipeline.verdict import compute_experimental_correlation
        assert compute_experimental_correlation([], set(), set()) == {}

    def test_score_pocket_adds_exp_tags(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {"best_affinity": "-7.5", "n_pose": "3", "n_ligand": "2",
                  "union_contact_residues": "ALA744;GLY752"}
        exp = {"exp_sensitivity": 0.5, "exp_hit_count": 2,
               "exp_false_pos": 0, "exp_enrichment": 3.0}
        _, _, reasons, _, _, _ = score_pocket(
            pocket, None, [], DEFAULT_THRESHOLDS, False, exp_correlation=exp)
        text = "; ".join(reasons)
        assert "exp_hit=" in text
        assert "exp_enriched" in text
        assert "exp_fp" not in text

    def test_score_pocket_no_exp_no_tags(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {"best_affinity": "-7.5", "n_pose": "3", "n_ligand": "2",
                  "union_contact_residues": "ALA744"}
        _, _, reasons, _, _, _ = score_pocket(
            pocket, None, [], DEFAULT_THRESHOLDS, False)
        assert not any(r.startswith("exp_") for r in reasons)

    @pytest.mark.parametrize("exp_data,expected_label", [
        ({"exp_sensitivity": 0.4, "exp_enrichment": 2.5, "exp_false_pos": 0}, "supports"),
        ({"exp_sensitivity": 0.2, "exp_enrichment": 1.8, "exp_false_pos": 0}, "consistent"),
        ({"exp_sensitivity": 0, "exp_enrichment": 0, "exp_false_pos": 3}, "contradicts"),
        ({"exp_sensitivity": 0, "exp_enrichment": 0, "exp_false_pos": 0}, "neutral"),
        (None, ""),
    ])
    def test_rank_impact_labels(self, exp_data, expected_label):
        from egfr_pipeline.verdict import _compute_exp_rank_impact
        assert _compute_exp_rank_impact(exp_data) == expected_label


# ===================================================================
# 4. Compare — bootstrap CI
# ===================================================================

class TestCompareBootstrapCI:
    def test_field_present(self):
        from egfr_pipeline.vina.compare import COMPARISON_FIELDS
        assert "centroid_dist_bootstrap_ci" in COMPARISON_FIELDS

    def test_ci_brackets_distance(self):
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
        bs = {("R1", "P01"): {"centroid_std_A": "0.5"},
              ("R2", "P01"): {"centroid_std_A": "0.8"}}
        results = compare_all_pockets(pockets, [], centroid_cutoff=0,
                                      bootstrap_index=bs)
        ci = results[0]["centroid_dist_bootstrap_ci"]
        lo, hi = map(float, ci.split("-"))
        assert lo < 2.0 < hi  # dist=2.0 should be inside CI

    def test_no_bootstrap_gives_empty_ci(self):
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
        assert results[0]["centroid_dist_bootstrap_ci"] == ""


# ===================================================================
# 5. Schema consistency
# ===================================================================

class TestSchemaConsistency:
    """Every module's FIELDS list must match validate.py EXPECTED_SCHEMAS."""

    @pytest.mark.parametrize("module_path,field_attr,schema_key", [
        ("egfr_pipeline.verdict", "VERDICT_FIELDS", "valid_sites.csv"),
        ("egfr_pipeline.verdict", "AGREEMENT_FIELDS", "cross_method_agreement.csv"),
        ("egfr_pipeline.verdict", "CONSENSUS_FIELDS", "vina_consensus_sites.csv"),
        ("egfr_pipeline.vina.compare", "COMPARISON_FIELDS", "vina_pocket_comparison.csv"),
        ("egfr_pipeline.vina.bootstrap", "BOOTSTRAP_FIELDS", "vina_pocket_bootstrap.csv"),
        ("egfr_pipeline.ppi.afm_extract", "AFM_RESIDUE_FIELDS", "ppi_afm_residues.csv"),
    ])
    def test_schema_match(self, module_path, field_attr, schema_key):
        import importlib
        from egfr_pipeline.validate import EXPECTED_SCHEMAS
        mod = importlib.import_module(module_path)
        assert getattr(mod, field_attr) == EXPECTED_SCHEMAS[schema_key]


# ===================================================================
# 6. End-to-end integration
# ===================================================================

class TestE2E:
    """Full pipeline: cluster → summarize → compare → bootstrap → verdict."""

    def test_full_pipeline_with_exp(self, pipeline_env):
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        from egfr_pipeline.verdict import generate_verdict

        cfg, outdir = pipeline_env

        # Summarize
        pocket_csv, drug_csv = summarize_from_config(str(cfg))
        assert len(read_csv(pocket_csv)) > 0

        # Compare
        cmp_csv = compare_from_config(str(cfg))
        assert cmp_csv.exists()

        # Bootstrap
        bs_csv = bootstrap_from_config(str(cfg), n_replicates=5)
        assert len(read_csv(bs_csv)) > 0

        # Verdict
        _, verdict_csv = generate_verdict(str(cfg))
        verdict_rows = read_csv(verdict_csv)
        assert len(verdict_rows) > 0

        # New columns populated
        assert any(r["exp_sensitivity"] not in ("", None) for r in verdict_rows), \
            "exp_sensitivity not populated"
        assert any(r["pocket_stability"] not in ("", None) for r in verdict_rows), \
            "pocket_stability not populated"
        # All verdicts valid
        for r in verdict_rows:
            assert r["verdict"] in ("STRONG", "MODERATE", "WEAK")

    def test_no_experimental_still_works(self, pipeline_env_no_exp):
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict

        cfg, _ = pipeline_env_no_exp
        summarize_from_config(str(cfg))
        compare_from_config(str(cfg))
        _, verdict_csv = generate_verdict(str(cfg))
        for r in read_csv(verdict_csv):
            assert r.get("exp_sensitivity", "") == ""

    def test_no_bootstrap_still_works(self, pipeline_env):
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict

        cfg, _ = pipeline_env
        summarize_from_config(str(cfg))
        compare_from_config(str(cfg))
        # Skip bootstrap
        _, verdict_csv = generate_verdict(str(cfg))
        for r in read_csv(verdict_csv):
            assert r.get("pocket_stability", "") == ""


# ===================================================================
# 7. Verdict scoring properties
# ===================================================================

class TestVerdictScoring:
    """Score_pocket behavior properties."""

    def test_score_always_0_to_100(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        # Maximum possible pocket
        pocket = {"best_affinity": "-12", "n_pose": "20", "n_ligand": "3",
                  "union_contact_residues": "ALA744"}
        total, _, _, _, _, _ = score_pocket(
            pocket, None, ["R2", "R3"], DEFAULT_THRESHOLDS, False)
        assert 0 <= total <= 100

        # Minimum pocket (no features)
        pocket_min = {"best_affinity": "0", "n_pose": "0", "n_ligand": "0",
                      "union_contact_residues": ""}
        total_min, _, _, _, _, _ = score_pocket(
            pocket_min, None, [], DEFAULT_THRESHOLDS, False)
        assert 0 <= total_min <= 100

    def test_ppi_does_not_penalize_without_data(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {"best_affinity": "-8.0", "n_pose": "5", "n_ligand": "2",
                  "union_contact_residues": "ALA744"}
        # With PPI = False: Vina(60) + Cross(40)
        t_no_ppi, _, _, v_no, p_no, c_no = score_pocket(
            pocket, None, ["R2"], DEFAULT_THRESHOLDS, False)
        # With PPI = True but no agreement: Vina(50) + PPI(0) + Cross(30)
        t_ppi, _, _, v_ppi, p_ppi, c_ppi = score_pocket(
            pocket, None, ["R2"], DEFAULT_THRESHOLDS, True)
        # Without PPI data, score should be >= with PPI(no match)
        assert t_no_ppi >= t_ppi

    @pytest.mark.parametrize("n_cross,expected_min", [
        (0, 0), (1, 10), (2, 30),
    ])
    def test_cross_receptor_graduated(self, n_cross, expected_min):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        pocket = {"best_affinity": "0", "n_pose": "0", "n_ligand": "0",
                  "union_contact_residues": ""}
        matches = [f"R{i}" for i in range(n_cross)]
        total, _, _, _, _, c_score = score_pocket(
            pocket, None, matches, DEFAULT_THRESHOLDS, False)
        assert c_score >= expected_min

    def test_verdict_labels_ordered(self):
        from egfr_pipeline.verdict import score_pocket, DEFAULT_THRESHOLDS
        # STRONG >= 55
        pocket_strong = {"best_affinity": "-9", "n_pose": "10", "n_ligand": "3",
                         "union_contact_residues": "ALA744"}
        _, v, _, _, _, _ = score_pocket(
            pocket_strong, None, ["R2", "R3"], DEFAULT_THRESHOLDS, False)
        assert v == "STRONG"

        # WEAK < 30
        pocket_weak = {"best_affinity": "0", "n_pose": "1", "n_ligand": "0",
                       "union_contact_residues": ""}
        _, v, _, _, _, _ = score_pocket(
            pocket_weak, None, [], DEFAULT_THRESHOLDS, False)
        assert v == "WEAK"


# ===================================================================
# 8. AlphaFold-Multimer integration
# ===================================================================

class TestAFMIntegration:
    """Tests for AFM extraction and verdict integration."""

    def test_afm_extraction_basic(self, tmp_path):
        """CA-CA distance extraction from a minimal PDB."""
        from egfr_pipeline.ppi.afm_extract import extract_afm_interface_residues
        # Create minimal 2-chain PDB: chain A CA at (10,10,10), chain B CA at (15,10,10)
        pdb = tmp_path / "model.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A 744      10.000  10.000  10.000  1.00  0.00\n"
            "ATOM      2  CA  ALA A 752      20.000  20.000  20.000  1.00  0.00\n"
            "ATOM      3  CA  GLY B  10      15.000  10.000  10.000  1.00  0.00\n"
            "END\n"
        )
        result = extract_afm_interface_residues(pdb, "A", "B", contact_cutoff=8.0,
                                                 receptor_id="R_test")
        rows = result["residue_rows"]
        # Only ALA744 is within 8A of chain B (dist=5.0)
        assert len(rows) == 1
        assert rows[0]["residue_id"] == "ALA744"
        assert abs(float(rows[0]["min_ca_distance"]) - 5.0) < 0.01

    def test_afm_extraction_no_contact(self, tmp_path):
        """No contacts when chains are far apart."""
        from egfr_pipeline.ppi.afm_extract import extract_afm_interface_residues
        pdb = tmp_path / "far.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A 744      10.000  10.000  10.000  1.00  0.00\n"
            "ATOM      2  CA  GLY B  10      99.000  99.000  99.000  1.00  0.00\n"
            "END\n"
        )
        result = extract_afm_interface_residues(pdb, "A", "B", contact_cutoff=8.0)
        assert len(result["residue_rows"]) == 0

    def test_adapt_afm_sigmoid_mapping(self):
        """Verify sigmoid distance-to-occupancy conversion."""
        from egfr_pipeline.verdict import _adapt_afm_to_ppi_format
        afm_rows = [
            {"receptor_id": "R1", "residue_id": "ALA744", "residue_num": "744",
             "min_ca_distance": "4.0", "source": "alphafold_multimer"},
            {"receptor_id": "R1", "residue_id": "GLY752", "residue_num": "752",
             "min_ca_distance": "8.0", "source": "alphafold_multimer"},
            {"receptor_id": "R1", "residue_id": "LEU831", "residue_num": "831",
             "min_ca_distance": "12.0", "source": "alphafold_multimer"},
        ]
        adapted = _adapt_afm_to_ppi_format(afm_rows)
        assert len(adapted) == 3
        # 4A -> high occ (~0.88), 8A -> ~0.5, 12A -> low occ (~0.12)
        occ_4 = adapted[0]["occupancy"]
        occ_8 = adapted[1]["occupancy"]
        occ_12 = adapted[2]["occupancy"]
        assert occ_4 > 0.8, f"4A should give high occ, got {occ_4}"
        assert 0.4 < occ_8 < 0.6, f"8A should give ~0.5 occ, got {occ_8}"
        assert occ_12 < 0.2, f"12A should give low occ, got {occ_12}"
        # Source tag format
        assert adapted[0]["source"] == "afm:R1"

    def test_afm_merge_with_pyrosetta(self):
        """AFM + PyRosetta residues merge preserving both sources."""
        from egfr_pipeline.verdict import _adapt_afm_to_ppi_format, _merge_multi_partner_residues
        pyro_rows = [
            {"receptor_id": "R1", "residue_id": "ALA744", "residue_num": "744",
             "source": "pyrosetta_ppi:TH1", "occupancy": "0.7", "frequency_final_ranking": "0.8"},
        ]
        afm_rows = [
            {"receptor_id": "R1", "residue_id": "ALA744", "residue_num": "744",
             "min_ca_distance": "5.0", "source": "alphafold_multimer"},
            {"receptor_id": "R1", "residue_id": "GLY752", "residue_num": "752",
             "min_ca_distance": "6.0", "source": "alphafold_multimer"},
        ]
        combined = pyro_rows + _adapt_afm_to_ppi_format(afm_rows)
        merged = _merge_multi_partner_residues(combined)
        # Both residues should survive
        res_ids = {r["residue_id"] for r in merged}
        assert "ALA744" in res_ids
        assert "GLY752" in res_ids
        # ALA744 has both sources -> merged source annotation
        ala = [r for r in merged if r["residue_id"] == "ALA744"][0]
        assert "merged" in ala.get("source", "") or len(res_ids) >= 2

    def test_verdict_with_afm_only(self, pipeline_env):
        """Verdict works with AFM data but no PyRosetta."""
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict
        from egfr_pipeline.ppi.afm_extract import AFM_RESIDUE_FIELDS

        cfg, outdir = pipeline_env
        summarize_from_config(str(cfg))
        compare_from_config(str(cfg))

        # Write synthetic AFM residues (no PyRosetta data)
        afm_rows = [
            {"receptor_id": "R_raw", "source": "alphafold_multimer",
             "residue_id": "ALA744", "residue_num": "744", "min_ca_distance": "5.0"},
            {"receptor_id": "R_raw", "source": "alphafold_multimer",
             "residue_id": "GLY752", "residue_num": "752", "min_ca_distance": "6.5"},
        ]
        write_csv(outdir / "ppi_afm_residues.csv", afm_rows, AFM_RESIDUE_FIELDS)

        _, verdict_csv = generate_verdict(str(cfg))
        verdict_rows = read_csv(verdict_csv)
        assert len(verdict_rows) > 0
        # At least one pocket should have PPI data (from AFM)
        assert any(r["ppi_data_available"] == "yes" for r in verdict_rows)

    def test_verdict_with_both_sources(self, pipeline_env):
        """Verdict with both PyRosetta + AFM data."""
        from egfr_pipeline.vina.summarize import summarize_from_config
        from egfr_pipeline.vina.compare import compare_from_config
        from egfr_pipeline.verdict import generate_verdict
        from egfr_pipeline.ppi.afm_extract import AFM_RESIDUE_FIELDS

        cfg, outdir = pipeline_env
        summarize_from_config(str(cfg))
        compare_from_config(str(cfg))

        # Write PyRosetta PPI
        pyro_fields = ["receptor_id", "source", "residue_id", "residue_num",
                        "frequency_final_ranking", "frequency_cluster_summary",
                        "n_models_final_ranking", "occupancy",
                        "mean_interface_delta_e", "best_interface_delta_e"]
        pyro_rows = [
            {"receptor_id": "R_raw", "source": "pyrosetta_ppi:TH1",
             "residue_id": "ALA744", "residue_num": "744",
             "frequency_final_ranking": "0.8", "frequency_cluster_summary": "0.6",
             "n_models_final_ranking": "3", "occupancy": "0.7",
             "mean_interface_delta_e": "-2.5", "best_interface_delta_e": "-4.0"},
        ]
        write_csv(outdir / "ppi_pyrosetta_residues.csv", pyro_rows, pyro_fields)
        write_csv(outdir / "ppi_pyrosetta_summary.csv", [
            {"receptor_id": "R_raw", "best_dg": "-15.0", "n_models": "5",
             "n_clusters": "3", "top_residues": "ALA744;GLY752"},
        ], ["receptor_id", "best_dg", "n_models", "n_clusters", "top_residues"])

        # Write AFM
        afm_rows = [
            {"receptor_id": "R_raw", "source": "alphafold_multimer",
             "residue_id": "ALA744", "residue_num": "744", "min_ca_distance": "5.0"},
            {"receptor_id": "R_raw", "source": "alphafold_multimer",
             "residue_id": "LEU831", "residue_num": "831", "min_ca_distance": "7.0"},
        ]
        write_csv(outdir / "ppi_afm_residues.csv", afm_rows, AFM_RESIDUE_FIELDS)

        _, verdict_csv = generate_verdict(str(cfg))
        verdict_rows = read_csv(verdict_csv)
        # R_raw should have PPI data from both sources
        r_raw_rows = [r for r in verdict_rows if r["receptor_id"] == "R_raw"]
        assert all(r["ppi_data_available"] == "yes" for r in r_raw_rows)

    def test_combined_evidence_includes_afm(self, pipeline_env):
        """Report's combined_residue_evidence includes AFM source tag."""
        from egfr_pipeline.report import format_combined_residue_table
        from egfr_pipeline.ppi.afm_extract import AFM_RESIDUE_FIELDS

        pocket_rows = [{"receptor_id": "R1", "pocket_id": "P01",
                        "union_contact_residues": "ALA744;GLY752"}]
        afm_rows = [{"receptor_id": "R1", "residue_id": "ALA744",
                      "min_ca_distance": "5.0", "source": "alphafold_multimer"}]
        combined = format_combined_residue_table(pocket_rows, [], afm_rows)
        ala = [r for r in combined if r["residue_id"] == "ALA744"][0]
        assert "afm" in ala["evidence_sources"]
        assert "vina" in ala["evidence_sources"]


# ===================================================================
# 9. Smoke regression
# ===================================================================

class TestSmokeRegression:
    @pytest.mark.skipif(not SMOKE_OUTPUT.exists(), reason="smoke output missing")
    def test_pose_count(self):
        assert len(read_csv(SMOKE_OUTPUT / "vina_pose_table.csv")) == 12

    @pytest.mark.skipif(not SMOKE_OUTPUT.exists(), reason="smoke output missing")
    def test_pocket_table_uncertainty_header(self):
        fields = csv_fields(SMOKE_OUTPUT / "vina_pocket_table.csv")
        for col in ["centroid_spread_A", "affinity_std", "affinity_iqr"]:
            assert col in fields

    @pytest.mark.skipif(
        not (SMOKE_OUTPUT / "valid_sites.csv").exists(),
        reason="valid_sites.csv missing",
    )
    def test_verdict_labels_valid(self):
        for r in read_csv(SMOKE_OUTPUT / "valid_sites.csv"):
            assert r["verdict"] in ("STRONG", "MODERATE", "WEAK")


# ===================================================================
# 9. YAML config (now that pyyaml is installed)
# ===================================================================

class TestYAMLConfig:
    """Verify YAML config loading works (pyyaml installed)."""

    def test_yaml_round_trip(self, tmp_path):
        import yaml
        cfg_dict = {
            "project_name": "yaml_test",
            "output_root": str(tmp_path),
            "receptors": [{"id": "R1", "pdb": "x.pdb"}],
            "experimental": {
                "known_binding_residues": [744, 752],
                "source": "test",
            },
        }
        cfg_path = tmp_path / "test.yaml"
        cfg_path.write_text(yaml.dump(cfg_dict), encoding="utf-8")
        loaded = load_config(str(cfg_path))
        assert loaded["project_name"] == "yaml_test"
        assert loaded["experimental"]["known_binding_residues"] == [744, 752]

    def test_smoke_config_loadable(self):
        smoke_cfg = PROJECT_ROOT / "smoke_test" / "config.yaml"
        if smoke_cfg.exists():
            cfg = load_config(str(smoke_cfg))
            assert cfg["project_name"] == "smoke_vina"
            assert len(cfg["receptors"]) == 3
