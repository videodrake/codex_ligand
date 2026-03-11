"""Comprehensive tests for Phase 3 pipeline (TG 3.0–3.7).

Tests cover:
  - Unit tests for saturation/budget/diversity logic
  - Schema validation for all CSV outputs
  - Cross-module data flow (field name consistency)
  - Edge cases (empty inputs, single pocket, no dockable pockets)
  - End-to-end pipeline integration
"""

import csv
import json
import math
import shutil
from pathlib import Path
from typing import Dict, List

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _load_csv(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _csv_columns(path: Path) -> List[str]:
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        return next(reader)


# ---------------------------------------------------------------------------
# Fixtures: synthetic Phase 2 → Phase 3 test data
# ---------------------------------------------------------------------------

PHASE3_POCKET_REFERENCE_HEADER = (
    "receptor_id,pocket_id,centroid_x,centroid_y,centroid_z,"
    "box_size_x,box_size_y,box_size_z,relationship_class,"
    "hotspot_overlap_count,hotspot_overlap_fraction,"
    "centroid_distance_to_patch_A,druggability_confidence,"
    "overall_druggability_tier,best_proposal_score,"
    "best_druggability_score,normalized_druggability_score,"
    "n_sources,sources,multi_source_support,state_class,"
    "n_states_matched,n_residues,residue_ids,mean_volume_A3,"
    "docking_priority,priority_basis"
)


@pytest.fixture
def tmp_phase3(tmp_path):
    """Create full Phase 3 test environment with multi-state multi-pocket data."""
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    ligands = tmp_path / "ligands"
    receptors = tmp_path / "receptors"
    phase2.mkdir()
    phase3.mkdir()
    ligands.mkdir()
    receptors.mkdir()

    # Phase 2 export (3 pockets across 2 states)
    _write(phase2 / "phase3_candidate_pocket_reference.csv", [
        PHASE3_POCKET_REFERENCE_HEADER,
        # 3GT8_raw PKT01: tier_1 orthosteric → primary
        "3GT8_raw,3GT8_raw_PKT01,44.5,32.3,18.6,14.0,13.5,13.0,"
        "orthosteric_candidate,4,0.571,132.5,high,tier_1,0.82,0.71,0.85,"
        "2,fpocket;p2rank,True,state_robust_pocket,2,4,"
        "LEU838;ASP855;ILE857;GLU866,890.0,primary,tier_1 + orthosteric_candidate",
        # 3GT8_raw PKT02: tier_2 orthosteric → primary
        "3GT8_raw,3GT8_raw_PKT02,28.5,45.5,12.6,11.0,11.5,12.0,"
        "orthosteric_candidate,3,0.429,124.4,medium,tier_2,0.56,0.44,0.50,"
        "1,fpocket,False,state_specific_pocket,1,3,"
        "LEU819;ALA822;VAL834,570.0,primary,tier_2 + orthosteric_candidate",
        # 3GT8_raw PKT03: tier_3 low_relevance → skip
        "3GT8_raw,3GT8_raw_PKT03,60.8,25.2,31.0,11.3,10.6,10.4,"
        "low_relevance_candidate,0,0.000,150.0,low,tier_3,0.35,0.20,0.10,"
        "1,fpocket,False,state_specific_pocket,1,2,"
        "ARG748;GLY750,300.0,skip,tier_3 + low_relevance (not recommended)",
        # 3GT8_cl38_48 PKT01: tier_2 rim → secondary
        "3GT8_cl38_48,3GT8_cl38_48_PKT01,45.0,33.0,19.0,13.0,12.5,12.0,"
        "rim_candidate,1,0.143,135.0,medium,tier_2,0.75,0.62,0.55,"
        "1,fpocket,False,state_specific_pocket,1,3,"
        "LEU838;ASP855;GLU866,850.0,secondary,tier_2 + rim_candidate",
    ])

    # Ligand files
    for lig in ["lig_A", "lig_B"]:
        (ligands / f"{lig}_ligand.sdf").write_text("mock SDF\n")

    # Receptor files
    for rec in ["3GT8_raw", "3GT8_cl38_48"]:
        (receptors / f"{rec}.pdb").write_text("mock PDB\n")

    return {
        "phase2": phase2,
        "phase3": phase3,
        "ligands": ligands,
        "receptors": receptors,
    }


@pytest.fixture
def tmp_phase3_single_pocket(tmp_path):
    """Phase 3 with only 1 dockable pocket (edge case)."""
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    ligands = tmp_path / "ligands"
    receptors = tmp_path / "receptors"
    phase2.mkdir()
    phase3.mkdir()
    ligands.mkdir()
    receptors.mkdir()

    _write(phase2 / "phase3_candidate_pocket_reference.csv", [
        PHASE3_POCKET_REFERENCE_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,44.5,32.3,18.6,14.0,13.5,13.0,"
        "orthosteric_candidate,4,0.571,132.5,high,tier_1,0.82,0.71,0.85,"
        "1,fpocket,False,state_specific_pocket,1,4,"
        "LEU838;ASP855;ILE857;GLU866,890.0,primary,tier_1 + orthosteric_candidate",
    ])

    (ligands / "lig_A_ligand.sdf").write_text("mock\n")
    (receptors / "3GT8_raw.pdb").write_text("mock\n")

    return {
        "phase2": phase2, "phase3": phase3,
        "ligands": ligands, "receptors": receptors,
    }


@pytest.fixture
def tmp_phase3_no_dockable(tmp_path):
    """Phase 3 where all pockets are skip (edge case)."""
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    ligands = tmp_path / "ligands"
    receptors = tmp_path / "receptors"
    phase2.mkdir()
    phase3.mkdir()
    ligands.mkdir()
    receptors.mkdir()

    _write(phase2 / "phase3_candidate_pocket_reference.csv", [
        PHASE3_POCKET_REFERENCE_HEADER,
        "3GT8_raw,3GT8_raw_PKT03,60.8,25.2,31.0,11.3,10.6,10.4,"
        "low_relevance_candidate,0,0.000,150.0,low,tier_3,0.35,0.20,0.10,"
        "1,fpocket,False,state_specific_pocket,1,2,"
        "ARG748;GLY750,300.0,skip,tier_3 + low_relevance (not recommended)",
    ])

    (ligands / "lig_A_ligand.sdf").write_text("mock\n")
    (receptors / "3GT8_raw.pdb").write_text("mock\n")

    return {
        "phase2": phase2, "phase3": phase3,
        "ligands": ligands, "receptors": receptors,
    }


# ---------------------------------------------------------------------------
# Helper: run full Phase 3 pipeline on fixture data
# ---------------------------------------------------------------------------

def _run_phase3_pipeline(env):
    """Run TG 3.0 → 3.7 on a test environment. Returns output dir."""
    from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
    from egfr_pipeline.phase3.job_construction import run_job_construction
    from egfr_pipeline.phase3.budget_policy import run_budget_policy
    from egfr_pipeline.phase3.run_diverse_docking import run_setup
    from egfr_pipeline.phase3.pose_attribution import run_pose_attribution
    from egfr_pipeline.phase3.diversity_validation import run_diversity_validation
    from egfr_pipeline.phase3.phase4_export import run_phase4_export
    from egfr_pipeline.phase3.review_report import run_review_report

    run_pocket_reference_ingestion(env["phase2"], env["phase3"])
    run_job_construction(env["phase3"], env["ligands"], env["receptors"])
    run_budget_policy(env["phase3"])
    run_setup(env["phase3"])
    run_pose_attribution(env["phase3"])
    run_diversity_validation(env["phase3"])
    run_phase4_export(env["phase3"])
    run_review_report(env["phase3"])

    return env["phase3"]


# ===========================================================================
# TG 3.0: Pocket Reference Ingestion
# ===========================================================================

class TestPocketReferenceIngestion:
    """Tests for TG 3.0."""

    def test_normalized_csv_created(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        norm, val = run_pocket_reference_ingestion(tmp_phase3["phase2"], tmp_phase3["phase3"])
        assert norm.exists()
        assert val.exists()

    def test_normalized_row_count(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        run_pocket_reference_ingestion(tmp_phase3["phase2"], tmp_phase3["phase3"])
        rows = _load_csv(tmp_phase3["phase3"] / "phase3_candidate_reference_normalized.csv")
        assert len(rows) == 4  # 3 raw + 1 cl38_48

    def test_dockable_flag(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        run_pocket_reference_ingestion(tmp_phase3["phase2"], tmp_phase3["phase3"])
        rows = _load_csv(tmp_phase3["phase3"] / "phase3_candidate_reference_normalized.csv")
        dockable = [r for r in rows if r["dockable"] == "True"]
        skip = [r for r in rows if r["dockable"] == "False"]
        assert len(dockable) == 3  # PKT01, PKT02 (raw) + PKT01 (cl38_48)
        assert len(skip) == 1  # PKT03

    def test_budget_defaults_applied(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        run_pocket_reference_ingestion(tmp_phase3["phase2"], tmp_phase3["phase3"])
        rows = _load_csv(tmp_phase3["phase3"] / "phase3_candidate_reference_normalized.csv")
        primary = [r for r in rows if r["docking_priority"] == "primary"]
        secondary = [r for r in rows if r["docking_priority"] == "secondary"]
        for r in primary:
            assert r["recommended_exhaustiveness"] == "32"
            assert r["recommended_seeds"] == "3"
        for r in secondary:
            assert r["recommended_exhaustiveness"] == "16"
            assert r["recommended_seeds"] == "1"

    def test_validation_report_content(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        _, val = run_pocket_reference_ingestion(tmp_phase3["phase2"], tmp_phase3["phase3"])
        text = val.read_text()
        assert "PASS" in text
        assert "V1_NON_EMPTY" in text

    def test_missing_phase2_raises(self, tmp_path):
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        with pytest.raises(FileNotFoundError):
            run_pocket_reference_ingestion(tmp_path / "nonexistent", tmp_path / "out")


# ===========================================================================
# TG 3.1: Job Construction
# ===========================================================================

class TestJobConstruction:
    """Tests for TG 3.1."""

    def test_job_count(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        # 2 primary pockets (3 seeds each) + 1 secondary (1 seed) = 7 pocket-seeds
        # x 2 ligands = 14 jobs
        assert len(jobs) == 14

    def test_no_cross_receptor_mixing(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        for j in jobs:
            assert j["pocket_id"].startswith(j["receptor_id"])

    def test_skip_pocket_excluded(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        pocket_ids = set(j["pocket_id"] for j in jobs)
        assert "3GT8_raw_PKT03" not in pocket_ids

    def test_job_id_deterministic(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        ids = [j["job_id"] for j in jobs]
        assert len(ids) == len(set(ids)), "Job IDs must be unique"

    def test_box_table_includes_all_pockets(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        boxes = _load_csv(out / "phase3_job_box_table.csv")
        assert len(boxes) == 4  # All pockets including skip

    def test_docking_mode_focused(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        assert all(j["docking_mode"] == "focused" for j in jobs)


# ===========================================================================
# TG 3.2: Budget Policy — Unit Tests
# ===========================================================================

class TestBudgetPolicyUnit:
    """Unit tests for saturation and budget logic."""

    def test_saturation_triggers(self):
        from egfr_pipeline.phase3.budget_policy import evaluate_saturation, DEFAULT_BUDGET_POLICY
        state = {
            "acceptable_poses": 12,
            "total_poses_generated": 20,
            "seeds_completed": 2,
            "seeds_allocated": 3,
        }
        status, reason = evaluate_saturation(state, DEFAULT_BUDGET_POLICY)
        assert status == "saturated"
        assert "12" in reason

    def test_no_premature_saturation(self):
        from egfr_pipeline.phase3.budget_policy import evaluate_saturation, DEFAULT_BUDGET_POLICY
        state = {
            "acceptable_poses": 10,
            "total_poses_generated": 3,  # Below min_poses_before_saturation
            "seeds_completed": 1,
            "seeds_allocated": 3,
        }
        status, _ = evaluate_saturation(state, DEFAULT_BUDGET_POLICY)
        assert status == "open"

    def test_exhausted_when_all_seeds_done(self):
        from egfr_pipeline.phase3.budget_policy import evaluate_saturation, DEFAULT_BUDGET_POLICY
        state = {
            "acceptable_poses": 5,
            "total_poses_generated": 8,
            "seeds_completed": 3,
            "seeds_allocated": 3,
        }
        status, reason = evaluate_saturation(state, DEFAULT_BUDGET_POLICY)
        assert status == "exhausted"
        assert "seeds" in reason

    def test_pose_acceptable_within_window(self):
        from egfr_pipeline.phase3.budget_policy import is_pose_acceptable, DEFAULT_BUDGET_POLICY
        assert is_pose_acceptable(-7.5, -8.0, DEFAULT_BUDGET_POLICY) is True
        assert is_pose_acceptable(-5.5, -8.0, DEFAULT_BUDGET_POLICY) is False

    def test_reallocation_from_saturated(self):
        from egfr_pipeline.phase3.budget_policy import compute_reallocation, DEFAULT_BUDGET_POLICY
        statuses = [
            {"pocket_id": "P1", "status": "saturated", "seeds_remaining": 2,
             "docking_priority": "primary", "max_seeds": 3, "seeds_allocated": 3},
            {"pocket_id": "P2", "status": "open", "seeds_remaining": 1,
             "docking_priority": "primary", "max_seeds": 5, "seeds_allocated": 1},
        ]
        actions = compute_reallocation(statuses, DEFAULT_BUDGET_POLICY)
        assert len(actions) >= 1
        assert actions[0]["to_pocket"] == "P2"

    def test_no_reallocation_when_disabled(self):
        from egfr_pipeline.phase3.budget_policy import compute_reallocation
        policy = {"enable_reallocation": False}
        statuses = [
            {"pocket_id": "P1", "status": "saturated", "seeds_remaining": 2,
             "docking_priority": "primary", "max_seeds": 3, "seeds_allocated": 3},
        ]
        assert compute_reallocation(statuses, policy) == []


class TestBudgetPolicyIntegration:
    """Integration tests for TG 3.2."""

    def test_pocket_statuses_created(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        from egfr_pipeline.phase3.budget_policy import run_budget_policy
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])
        run_budget_policy(out)

        statuses = _load_csv(out / "pocket_search_status.csv")
        assert len(statuses) == 4
        open_count = sum(1 for s in statuses if s["status"] == "open")
        skip_count = sum(1 for s in statuses if s["status"] == "skipped")
        assert open_count == 3
        assert skip_count == 1

    def test_policy_document_created(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        from egfr_pipeline.phase3.budget_policy import run_budget_policy
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])
        run_budget_policy(out)

        policy_md = out / "phase3_budget_policy.md"
        assert policy_md.exists()
        text = policy_md.read_text()
        assert "Saturation Rules" in text
        assert "Reallocation" in text


# ===========================================================================
# TG 3.3: Execution Layer
# ===========================================================================

class TestExecutionLayer:
    """Tests for TG 3.3 (setup mode only — no Vina available)."""

    def test_setup_creates_script(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        from egfr_pipeline.phase3.budget_policy import run_budget_policy
        from egfr_pipeline.phase3.run_diverse_docking import run_setup
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])
        run_budget_policy(out)
        run_setup(out)

        assert (out / "run_phase3_docking.sh").exists()
        assert (out / "phase3_run_metadata.json").exists()

    def test_metadata_json_valid(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        from egfr_pipeline.phase3.budget_policy import run_budget_policy
        from egfr_pipeline.phase3.run_diverse_docking import run_setup
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])
        run_budget_policy(out)
        run_setup(out)

        meta = json.loads((out / "phase3_run_metadata.json").read_text())
        assert meta["pipeline"] == "phase3_diverse_docking"
        assert meta["total_jobs"] == 14
        assert "3GT8_raw" in meta["receptors"]

    def test_dry_run_validates(self, tmp_phase3):
        out = tmp_phase3["phase3"]
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        from egfr_pipeline.phase3.job_construction import run_job_construction
        from egfr_pipeline.phase3.budget_policy import run_budget_policy
        from egfr_pipeline.phase3.run_diverse_docking import dry_run
        from egfr_pipeline.phase3.budget_policy import DEFAULT_BUDGET_POLICY
        run_pocket_reference_ingestion(tmp_phase3["phase2"], out)
        run_job_construction(out, tmp_phase3["ligands"], tmp_phase3["receptors"])
        run_budget_policy(out)

        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        statuses = _load_csv(out / "pocket_search_status.csv")
        messages = dry_run(jobs, statuses, dict(DEFAULT_BUDGET_POLICY))
        ok_msgs = [m for m in messages if m.startswith("[OK]")]
        assert len(ok_msgs) >= 3


# ===========================================================================
# TG 3.4: Pose Attribution
# ===========================================================================

class TestPoseAttribution:
    """Tests for TG 3.4."""

    def test_pose_table_created(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        assert (out / "vina_pose_table.csv").exists()

    def test_provenance_fields_present(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        cols = _csv_columns(out / "vina_pose_table.csv")
        for field in ["candidate_pocket_id", "round_id", "seed",
                      "docking_mode", "job_id"]:
            assert field in cols, f"Missing provenance field: {field}"

    def test_backward_compatible_schema(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        cols = _csv_columns(out / "vina_pose_table.csv")
        legacy = ["receptor_id", "ligand_id", "pose_rank", "affinity",
                   "rmsd_lb", "rmsd_ub", "centroid_x", "centroid_y",
                   "centroid_z", "raw_pose_file", "pocket_id",
                   "contact_residues", "n_contact_residues"]
        for field in legacy:
            assert field in cols, f"Missing legacy field: {field}"

    def test_all_rows_have_pocket_attribution(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        rows = _load_csv(out / "vina_pose_table.csv")
        for r in rows:
            assert r["candidate_pocket_id"], f"Row missing pocket: {r['job_id']}"

    def test_provenance_note_created(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        note = out / "phase3_pose_provenance_note.md"
        assert note.exists()
        assert "candidate_pocket_id" in note.read_text()


# ===========================================================================
# TG 3.5: Diversity Validation — Unit Tests
# ===========================================================================

class TestDiversityUnit:
    """Unit tests for diversity metrics."""

    def test_shannon_entropy_uniform(self):
        """Uniform distribution should give max entropy."""
        from egfr_pipeline.phase3.diversity_validation import compute_diversity_metrics
        occ = [
            {"receptor_id": "R", "ligand_id": "L", "candidate_pocket_id": "P1",
             "n_poses": "10", "docking_priority": "primary", "pocket_status": "open"},
            {"receptor_id": "R", "ligand_id": "L", "candidate_pocket_id": "P2",
             "n_poses": "10", "docking_priority": "primary", "pocket_status": "open"},
        ]
        norm = [
            {"pocket_id": "P1", "receptor_id": "R", "dockable": "True"},
            {"pocket_id": "P2", "receptor_id": "R", "dockable": "True"},
        ]
        metrics = compute_diversity_metrics(occ, norm)
        assert len(metrics) == 1
        assert float(metrics[0]["normalized_entropy"]) == 1.0

    def test_shannon_entropy_concentrated(self):
        """All poses in one pocket should give zero entropy."""
        from egfr_pipeline.phase3.diversity_validation import compute_diversity_metrics
        occ = [
            {"receptor_id": "R", "ligand_id": "L", "candidate_pocket_id": "P1",
             "n_poses": "20", "docking_priority": "primary", "pocket_status": "open"},
            {"receptor_id": "R", "ligand_id": "L", "candidate_pocket_id": "P2",
             "n_poses": "0", "docking_priority": "primary", "pocket_status": "open"},
        ]
        norm = [
            {"pocket_id": "P1", "receptor_id": "R", "dockable": "True"},
            {"pocket_id": "P2", "receptor_id": "R", "dockable": "True"},
        ]
        metrics = compute_diversity_metrics(occ, norm)
        assert float(metrics[0]["normalized_entropy"]) == 0.0
        assert metrics[0]["diversity_verdict"] == "over_concentrated"

    def test_verdict_single_pocket(self):
        from egfr_pipeline.phase3.diversity_validation import compute_diversity_metrics
        occ = [
            {"receptor_id": "R", "ligand_id": "L", "candidate_pocket_id": "P1",
             "n_poses": "10", "docking_priority": "primary", "pocket_status": "open"},
        ]
        norm = [
            {"pocket_id": "P1", "receptor_id": "R", "dockable": "True"},
        ]
        metrics = compute_diversity_metrics(occ, norm)
        assert metrics[0]["diversity_verdict"] == "single_pocket"


class TestDiversityIntegration:
    """Integration tests for TG 3.5."""

    def test_occupancy_covers_all_pocket_ligand_combos(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        occ = _load_csv(out / "phase3_pocket_occupancy_summary.csv")
        # 4 pockets x 2 ligands = 8 entries (only for receptors with jobs)
        pocket_ligand_pairs = set(
            (r["candidate_pocket_id"], r["ligand_id"]) for r in occ
        )
        assert len(pocket_ligand_pairs) >= 6  # At least dockable combos

    def test_diversity_metrics_per_receptor_ligand(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        div = _load_csv(out / "phase3_diversity_metrics.csv")
        # Should have metrics for each (receptor, ligand) pair
        assert len(div) >= 2


# ===========================================================================
# TG 3.6: Phase 4 Export — Unit Tests
# ===========================================================================

class TestPhase4ExportUnit:
    """Unit tests for ligand support classification."""

    def test_strong_support(self):
        from egfr_pipeline.phase3.phase4_export import _classify_support
        assert _classify_support(10, "-7.5") == "strong"

    def test_moderate_support(self):
        from egfr_pipeline.phase3.phase4_export import _classify_support
        assert _classify_support(3, "-5.5") == "moderate"

    def test_weak_support(self):
        from egfr_pipeline.phase3.phase4_export import _classify_support
        assert _classify_support(1, "-3.0") == "weak"

    def test_no_support(self):
        from egfr_pipeline.phase3.phase4_export import _classify_support
        assert _classify_support(0, "") == "none"

    def test_pending_support(self):
        from egfr_pipeline.phase3.phase4_export import _classify_support
        assert _classify_support(5, "") == "pending_strong"
        assert _classify_support(2, "") == "pending_moderate"
        assert _classify_support(1, "") == "pending_weak"


class TestPhase4ExportIntegration:
    """Integration tests for TG 3.6."""

    def test_evidence_created(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        evidence = _load_csv(out / "phase4_docking_evidence_reference.csv")
        assert len(evidence) > 0

    def test_evidence_excludes_skip(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        evidence = _load_csv(out / "phase4_docking_evidence_reference.csv")
        pockets = set(r["candidate_pocket_id"] for r in evidence)
        assert "3GT8_raw_PKT03" not in pockets

    def test_evidence_has_perturbation_fields(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        evidence = _load_csv(out / "phase4_docking_evidence_reference.csv")
        for field in ["relationship_class", "overall_druggability_tier",
                      "state_class", "ligand_support_strength"]:
            for r in evidence:
                assert r.get(field), f"Missing {field} in evidence row"

    def test_handoff_note_created(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        assert (out / "phase3_to_phase4_handoff_note.md").exists()


# ===========================================================================
# TG 3.7: Review Report
# ===========================================================================

class TestReviewReport:
    """Tests for TG 3.7."""

    def test_report_created(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        report = out / "phase3_diverse_docking_report.md"
        assert report.exists()

    def test_report_sections(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        text = (out / "phase3_diverse_docking_report.md").read_text()
        for section in ["Executive Summary", "Pocket Catalog",
                        "Search Budget Summary", "Ligand Support Summary",
                        "Diversity Outcome", "Phase 4 Handoff",
                        "Caveats", "File Inventory"]:
            assert section in text, f"Missing section: {section}"

    def test_report_lists_all_pockets(self, tmp_phase3):
        out = _run_phase3_pipeline(tmp_phase3)
        text = (out / "phase3_diverse_docking_report.md").read_text()
        assert "3GT8_raw_PKT01" in text
        assert "3GT8_raw_PKT02" in text
        assert "3GT8_cl38_48_PKT01" in text


# ===========================================================================
# Schema Consistency
# ===========================================================================

class TestSchemaConsistency:
    """Verify declared schemas match actual CSV output."""

    def test_normalized_reference_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.pocket_reference_ingestion import NORMALIZED_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase3_candidate_reference_normalized.csv")
        assert actual == NORMALIZED_COLUMNS

    def test_job_table_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.job_construction import JOB_TABLE_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase3_docking_job_table.csv")
        assert actual == JOB_TABLE_COLUMNS

    def test_box_table_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.job_construction import BOX_TABLE_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase3_job_box_table.csv")
        assert actual == BOX_TABLE_COLUMNS

    def test_search_status_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.budget_policy import SEARCH_STATUS_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "pocket_search_status.csv")
        assert actual == SEARCH_STATUS_COLUMNS

    def test_pose_table_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.pose_attribution import POSE_TABLE_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "vina_pose_table.csv")
        assert actual == POSE_TABLE_COLUMNS

    def test_occupancy_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.diversity_validation import OCCUPANCY_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase3_pocket_occupancy_summary.csv")
        assert actual == OCCUPANCY_COLUMNS

    def test_diversity_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.diversity_validation import DIVERSITY_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase3_diversity_metrics.csv")
        assert actual == DIVERSITY_COLUMNS

    def test_evidence_schema(self, tmp_phase3):
        from egfr_pipeline.phase3.phase4_export import EVIDENCE_COLUMNS
        out = _run_phase3_pipeline(tmp_phase3)
        actual = _csv_columns(out / "phase4_docking_evidence_reference.csv")
        assert actual == EVIDENCE_COLUMNS


# ===========================================================================
# End-to-End
# ===========================================================================

class TestEndToEnd:
    """Full pipeline integration tests."""

    def test_multi_state_multi_pocket(self, tmp_phase3):
        """Full pipeline with 2 states, 4 pockets, 2 ligands."""
        out = _run_phase3_pipeline(tmp_phase3)
        evidence = _load_csv(out / "phase4_docking_evidence_reference.csv")
        assert len(evidence) > 0
        receptors = set(r["receptor_id"] for r in evidence)
        assert "3GT8_raw" in receptors

    def test_single_pocket(self, tmp_phase3_single_pocket):
        """Pipeline with only 1 dockable pocket."""
        out = _run_phase3_pipeline(tmp_phase3_single_pocket)
        jobs = _load_csv(out / "phase3_docking_job_table.csv")
        assert len(jobs) == 3  # 1 pocket x 1 ligand x 3 seeds
        div = _load_csv(out / "phase3_diversity_metrics.csv")
        assert div[0]["diversity_verdict"] == "single_pocket"

    def test_idempotent(self, tmp_phase3):
        """Running pipeline twice produces same output."""
        out1 = _run_phase3_pipeline(tmp_phase3)
        jobs1 = _load_csv(out1 / "phase3_docking_job_table.csv")
        evidence1 = _load_csv(out1 / "phase4_docking_evidence_reference.csv")

        # Run again
        out2 = _run_phase3_pipeline(tmp_phase3)
        jobs2 = _load_csv(out2 / "phase3_docking_job_table.csv")
        evidence2 = _load_csv(out2 / "phase4_docking_evidence_reference.csv")

        assert len(jobs1) == len(jobs2)
        assert len(evidence1) == len(evidence2)
        for j1, j2 in zip(jobs1, jobs2):
            assert j1["job_id"] == j2["job_id"]

    def test_no_dockable_pockets(self, tmp_phase3_no_dockable):
        """Pipeline where all pockets are skip."""
        from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
        env = tmp_phase3_no_dockable
        _, val = run_pocket_reference_ingestion(env["phase2"], env["phase3"])
        text = val.read_text()
        # Should have validation error about no dockable pockets
        assert "V9_NO_DOCKABLE" in text or "FAIL" in text
