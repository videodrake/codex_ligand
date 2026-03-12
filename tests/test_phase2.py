"""Comprehensive tests for Phase 2 pipeline (TG 2.0–2.7).

Tests cover:
  - Unit tests for classification/scoring logic
  - Schema validation for all CSV outputs
  - Cross-module data flow (field name consistency)
  - Edge cases (empty inputs, single state, missing data)
  - End-to-end pipeline integration
"""

import csv
import math
import shutil
import textwrap
from pathlib import Path
from typing import Dict, List

import pytest

# ---------------------------------------------------------------------------
# Fixtures: synthetic test data
# ---------------------------------------------------------------------------

RECEPTOR_STATES = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]


@pytest.fixture
def tmp_phase2(tmp_path):
    """Create a temporary Phase 2 output directory with synthetic data."""
    out = tmp_path / "phase2"
    out.mkdir()

    # Phase 1 patch reference (TG 2.0 output)
    _write(out / "phase2_patch_reference_normalized.csv", [
        "patch_residue_id,residue_num,residue_name,chain,lobe,construct_type,"
        "orientation_validated,robustness,n_states,states,max_occupancy,"
        "is_hotspot,pyrosetta_support,lightdock_support,method_agreement,"
        "confidence,phase1_evidence_source,notes",
        "LEU838,838,LEU,A,C_lobe,full_kinase_domain,True,robust,3,"
        "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100,0.85,True,True,True,both,"
        "high,pyrosetta+lightdock,",
        "ASP855,855,ASP,A,C_lobe,full_kinase_domain,True,robust,3,"
        "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100,0.78,True,True,False,pyrosetta_only,"
        "medium,pyrosetta,",
        "ILE857,857,ILE,A,C_lobe,full_kinase_domain,True,moderate,2,"
        "3GT8_raw;3GT8_cl38_48,0.62,True,True,False,pyrosetta_only,"
        "medium,pyrosetta,",
        "GLU866,866,GLU,A,C_lobe,full_kinase_domain,True,robust,3,"
        "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100,0.71,True,True,True,both,"
        "high,pyrosetta+lightdock,",
        "LEU819,819,LEU,A,N_lobe,full_kinase_domain,True,moderate,2,"
        "3GT8_raw;3GT8_cl85_100,0.55,True,True,False,pyrosetta_only,"
        "medium,pyrosetta,",
        "ALA822,822,ALA,A,N_lobe,full_kinase_domain,True,state_specific,1,"
        "3GT8_raw,0.51,True,True,False,pyrosetta_only,low,pyrosetta,",
        "VAL834,834,VAL,A,N_lobe,full_kinase_domain,True,state_specific,1,"
        "3GT8_raw,0.50,False,True,False,pyrosetta_only,low,pyrosetta,",
    ])

    # Raw pockets (TG 2.1 output) — multi-state + multi-tool
    _write(out / "candidate_pockets_raw.csv", [
        "receptor_id,candidate_pocket_id,proposal_source,pocket_rank,proposal_score,"
        "druggability_score,centroid_x,centroid_y,centroid_z,volume_A3,n_alpha_spheres,"
        "residue_ids,residue_nums,n_residues,box_center_x,box_center_y,box_center_z,"
        "box_size_x,box_size_y,box_size_z",
        # 3GT8_raw: fpocket
        "3GT8_raw,3GT8_raw_fpocket_P01,fpocket,1,0.82,0.71,44.5,32.3,18.6,890.0,42,"
        "LEU838;ASP855;ILE857;GLU866,838;855;857;866,4,44.5,32.3,18.6,14.0,13.5,13.0",
        "3GT8_raw,3GT8_raw_fpocket_P02,fpocket,2,0.56,0.44,28.5,45.5,12.6,570.0,28,"
        "LEU819;ALA822;VAL834,819;822;834,3,28.5,45.5,12.6,11.0,11.5,12.0",
        # 3GT8_raw: p2rank (overlaps with fpocket P01)
        "3GT8_raw,3GT8_raw_p2rank_P01,p2rank,1,0.78,0.68,44.8,32.0,18.2,910.0,0,"
        "LEU838;ASP855;ILE857;GLU866;MET769,838;855;857;866;769,5,44.8,32.0,18.2,15.0,14.0,13.5",
        # 3GT8_cl38_48: fpocket (similar location to raw PKT01)
        "3GT8_cl38_48,3GT8_cl38_48_fpocket_P01,fpocket,1,0.75,0.62,45.0,33.0,19.0,850.0,38,"
        "LEU838;ASP855;GLU866,838;855;866,3,45.0,33.0,19.0,13.0,12.5,12.0",
        # 3GT8_cl38_48: fpocket (distant pocket)
        "3GT8_cl38_48,3GT8_cl38_48_fpocket_P02,fpocket,2,0.35,0.20,80.0,60.0,50.0,300.0,14,"
        "ARG748;GLY750,748;750,2,80.0,60.0,50.0,10.0,10.0,10.0",
    ])

    return out


@pytest.fixture
def tmp_phase2_single_state(tmp_path):
    """Phase 2 data with only one receptor state (edge case)."""
    out = tmp_path / "phase2_single"
    out.mkdir()

    _write(out / "phase2_patch_reference_normalized.csv", [
        "patch_residue_id,residue_num,residue_name,chain,lobe,construct_type,"
        "orientation_validated,robustness,n_states,states,max_occupancy,"
        "is_hotspot,pyrosetta_support,lightdock_support,method_agreement,"
        "confidence,phase1_evidence_source,notes",
        "LEU838,838,LEU,A,C_lobe,full_kinase_domain,True,robust,3,"
        "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100,0.85,True,True,True,both,"
        "high,pyrosetta+lightdock,",
    ])

    _write(out / "candidate_pockets_raw.csv", [
        "receptor_id,candidate_pocket_id,proposal_source,pocket_rank,proposal_score,"
        "druggability_score,centroid_x,centroid_y,centroid_z,volume_A3,n_alpha_spheres,"
        "residue_ids,residue_nums,n_residues,box_center_x,box_center_y,box_center_z,"
        "box_size_x,box_size_y,box_size_z",
        "3GT8_raw,3GT8_raw_fpocket_P01,fpocket,1,0.80,0.65,44.5,32.3,18.6,890.0,42,"
        "LEU838,838,1,44.5,32.3,18.6,14.0,13.5,13.0",
    ])

    return out


@pytest.fixture
def tmp_phase2_empty(tmp_path):
    """Phase 2 dir with no pocket data (edge case)."""
    out = tmp_path / "phase2_empty"
    out.mkdir()
    _write(out / "phase2_patch_reference_normalized.csv", [
        "patch_residue_id,residue_num,residue_name,chain,lobe,construct_type,"
        "orientation_validated,robustness,n_states,states,max_occupancy,"
        "is_hotspot,pyrosetta_support,lightdock_support,method_agreement,"
        "confidence,phase1_evidence_source,notes",
    ])
    _write(out / "candidate_pockets_raw.csv", [
        "receptor_id,candidate_pocket_id,proposal_source,pocket_rank,proposal_score,"
        "druggability_score,centroid_x,centroid_y,centroid_z,volume_A3,n_alpha_spheres,"
        "residue_ids,residue_nums,n_residues,box_center_x,box_center_y,box_center_z,"
        "box_size_x,box_size_y,box_size_z",
    ])
    return out


def _write(path: Path, lines: List[str]):
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _load(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ===================================================================
# TG 2.2: Pocket Merge
# ===================================================================

class TestPocketMerge:
    """Tests for egfr_pipeline.phase2.pocket_merge."""

    def test_merge_overlapping_pockets(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)

        merged = _load(tmp_phase2 / "candidate_pockets.csv")
        # 3GT8_raw: fpocket_P01 + p2rank_P01 should merge (close centroids)
        # 3GT8_raw: fpocket_P02 stays separate
        # 3GT8_cl38_48: 2 pockets, distant, stay separate
        raw_states = {"3GT8_raw", "3GT8_cl38_48"}
        merged_states = {p["receptor_id"] for p in merged}
        assert merged_states == raw_states

        # fpocket_P01 and p2rank_P01 in 3GT8_raw should merge
        raw_pkt01 = [p for p in merged
                     if p["receptor_id"] == "3GT8_raw"
                     and "fpocket" in p.get("sources", "")
                     and "p2rank" in p.get("sources", "")]
        assert len(raw_pkt01) >= 1, "fpocket+p2rank should merge into one pocket"

    def test_never_merge_across_states(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)

        prov = _load(tmp_phase2 / "candidate_pocket_provenance.csv")
        for row in prov:
            # merged_pocket_id starts with a state, raw_pocket_id starts with same state
            merged_state = row["merged_pocket_id"].split("_PKT")[0]
            raw_state = "_".join(row["raw_pocket_id"].split("_")[:-2])
            # Handle potential naming issues by checking receptor_id
            assert row["receptor_id"] in merged_state

    def test_stable_pocket_ids(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)
        merged = _load(tmp_phase2 / "candidate_pockets.csv")

        ids = [p["pocket_id"] for p in merged]
        # All IDs unique
        assert len(ids) == len(set(ids))
        # All follow {state}_PKT{nn} format
        for pid in ids:
            assert "_PKT" in pid

    def test_provenance_completeness(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)

        prov = _load(tmp_phase2 / "candidate_pocket_provenance.csv")
        raw = _load(tmp_phase2 / "candidate_pockets_raw.csv")

        # Every raw pocket appears in provenance
        raw_ids = {p["candidate_pocket_id"] for p in raw}
        prov_raw_ids = {p["raw_pocket_id"] for p in prov}
        assert raw_ids == prov_raw_ids

    def test_merge_table_written(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)

        mt = _load(tmp_phase2 / "candidate_pocket_merge_table.csv")
        assert len(mt) > 0
        decisions = {r["merge_decision"] for r in mt}
        assert decisions <= {"merged", "kept_separate"}

    def test_empty_input_raises(self, tmp_phase2_empty):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        # Empty raw file has header only → 0 rows → should still work (empty merge)
        run_pocket_merge(tmp_phase2_empty)
        merged = _load(tmp_phase2_empty / "candidate_pockets.csv")
        assert len(merged) == 0


# ===================================================================
# TG 2.3: Patch Relationship Classification
# ===================================================================

class TestPatchRelationship:
    """Tests for egfr_pipeline.phase2.patch_relationship."""

    def _setup(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        run_pocket_merge(tmp_phase2)

    def test_all_pockets_classified(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_patch_relationship(tmp_phase2)

        merged = _load(tmp_phase2 / "candidate_pockets.csv")
        rels = _load(tmp_phase2 / "pocket_patch_relationship.csv")
        assert len(rels) == len(merged)

    def test_valid_relationship_classes(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_patch_relationship(tmp_phase2)

        rels = _load(tmp_phase2 / "pocket_patch_relationship.csv")
        valid = {"orthosteric_candidate", "rim_candidate",
                 "allosteric_candidate", "low_relevance_candidate"}
        for r in rels:
            assert r["relationship_class"] in valid, f"Invalid class: {r['relationship_class']}"

    def test_orthosteric_requires_overlap(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_patch_relationship(tmp_phase2)

        rels = _load(tmp_phase2 / "pocket_patch_relationship.csv")
        for r in rels:
            if r["relationship_class"] == "orthosteric_candidate":
                assert int(r["hotspot_overlap_count"]) >= 2
                assert float(r["hotspot_overlap_fraction"]) >= 0.25

    def test_raw_metrics_preserved(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_patch_relationship(tmp_phase2)

        metrics = _load(tmp_phase2 / "pocket_patch_relationship_metrics.csv")
        assert len(metrics) > 0
        for m in metrics:
            assert "pocket_residue_jaccard" in m
            assert "centroid_distance_A" in m

    def test_classification_basis_not_empty(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_patch_relationship(tmp_phase2)

        rels = _load(tmp_phase2 / "pocket_patch_relationship.csv")
        for r in rels:
            assert r["classification_basis"], f"Empty basis for {r['pocket_id']}"


# ===================================================================
# TG 2.3 Unit: classify_pocket_patch_relationship
# ===================================================================

class TestClassifyUnit:
    """Unit tests for classification logic."""

    def test_orthosteric(self):
        from egfr_pipeline.phase2.patch_relationship import classify_pocket_patch_relationship
        pocket = {
            "pocket_id": "test_PKT01", "receptor_id": "3GT8_raw",
            "centroid_x": "10.0", "centroid_y": "10.0", "centroid_z": "10.0",
            "residue_ids": "LEU838;ASP855;ILE857",
        }
        hotspots = {"LEU838", "ASP855", "ILE857", "GLU866"}
        centroid = (10.0, 10.0, 10.0)
        result = classify_pocket_patch_relationship(pocket, hotspots, centroid)
        assert result["relationship_class"] == "orthosteric_candidate"
        assert int(result["hotspot_overlap_count"]) == 3

    def test_low_relevance_no_overlap_far(self):
        from egfr_pipeline.phase2.patch_relationship import classify_pocket_patch_relationship
        pocket = {
            "pocket_id": "test_PKT02", "receptor_id": "3GT8_raw",
            "centroid_x": "100.0", "centroid_y": "100.0", "centroid_z": "100.0",
            "residue_ids": "ARG748;GLY750",
        }
        hotspots = {"LEU838", "ASP855"}
        centroid = (10.0, 10.0, 10.0)
        result = classify_pocket_patch_relationship(pocket, hotspots, centroid)
        assert result["relationship_class"] == "low_relevance_candidate"

    def test_allosteric_near_no_overlap(self):
        from egfr_pipeline.phase2.patch_relationship import classify_pocket_patch_relationship
        pocket = {
            "pocket_id": "test_PKT03", "receptor_id": "3GT8_raw",
            "centroid_x": "15.0", "centroid_y": "10.0", "centroid_z": "10.0",
            "residue_ids": "ARG748;GLY750",
        }
        hotspots = {"LEU838", "ASP855"}
        centroid = (10.0, 10.0, 10.0)
        result = classify_pocket_patch_relationship(pocket, hotspots, centroid)
        assert result["relationship_class"] == "allosteric_candidate"

    def test_rim_partial_overlap(self):
        from egfr_pipeline.phase2.patch_relationship import classify_pocket_patch_relationship
        pocket = {
            "pocket_id": "test_PKT04", "receptor_id": "3GT8_raw",
            "centroid_x": "50.0", "centroid_y": "50.0", "centroid_z": "50.0",
            "residue_ids": "LEU838;ARG748",
        }
        hotspots = {"LEU838", "ASP855", "ILE857", "GLU866"}
        centroid = (10.0, 10.0, 10.0)
        result = classify_pocket_patch_relationship(pocket, hotspots, centroid)
        assert result["relationship_class"] == "rim_candidate"

    def test_no_centroid_available(self):
        from egfr_pipeline.phase2.patch_relationship import classify_pocket_patch_relationship
        pocket = {
            "pocket_id": "test_PKT05", "receptor_id": "3GT8_raw",
            "centroid_x": "", "centroid_y": "", "centroid_z": "",
            "residue_ids": "ARG748",
        }
        hotspots = {"LEU838", "ASP855"}
        result = classify_pocket_patch_relationship(pocket, hotspots, None)
        assert result["relationship_class"] == "low_relevance_candidate"


# ===================================================================
# TG 2.4: Druggability Confidence
# ===================================================================

class TestDruggabilityConfidence:
    """Tests for egfr_pipeline.phase2.druggability_confidence."""

    def _setup(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)

    def test_all_pockets_scored(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        merged = _load(tmp_phase2 / "candidate_pockets.csv")
        summary = _load(tmp_phase2 / "druggability_proposal_summary.csv")
        assert len(summary) == len(merged)

    def test_confidence_levels_valid(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        summary = _load(tmp_phase2 / "druggability_proposal_summary.csv")
        for r in summary:
            assert r["druggability_confidence"] in {"high", "medium", "low"}

    def test_tier_levels_valid(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        summary = _load(tmp_phase2 / "druggability_proposal_summary.csv")
        for r in summary:
            assert r["overall_druggability_tier"] in {"tier_1", "tier_2", "tier_3"}

    def test_multi_source_detected(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        flags = _load(tmp_phase2 / "candidate_pocket_support_flags.csv")
        # At least one pocket should have both fpocket and p2rank
        has_consensus = any(
            r["fpocket_support"] == "true" and r["p2rank_support"] == "true"
            for r in flags
        )
        assert has_consensus, "Merged fpocket+p2rank pocket should show consensus"

    def test_normalized_scores_in_range(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        summary = _load(tmp_phase2 / "druggability_proposal_summary.csv")
        for r in summary:
            if r["normalized_proposal_score"]:
                val = float(r["normalized_proposal_score"])
                assert 0.0 <= val <= 1.0, f"Normalized score out of range: {val}"

    def test_ftmap_columns_present_but_empty(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_druggability_confidence(tmp_phase2)

        summary = _load(tmp_phase2 / "druggability_proposal_summary.csv")
        for r in summary:
            assert "ftmap_hotspot_count" in r
            assert r["hotspot_support_available"] == "false"


# ===================================================================
# TG 2.4 Unit: scoring logic
# ===================================================================

class TestDruggabilityUnit:
    """Unit tests for druggability scoring functions."""

    def test_assign_confidence_high(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_druggability_confidence
        assert assign_druggability_confidence(0.65, None) == "high"

    def test_assign_confidence_medium(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_druggability_confidence
        assert assign_druggability_confidence(0.35, None) == "medium"

    def test_assign_confidence_low(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_druggability_confidence
        assert assign_druggability_confidence(0.10, None) == "low"

    def test_assign_confidence_none(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_druggability_confidence
        assert assign_druggability_confidence(None, None) == "low"

    def test_assign_tier_1(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_overall_tier
        assert assign_overall_tier("high", True, "orthosteric_candidate") == "tier_1"
        assert assign_overall_tier("high", False, "rim_candidate") == "tier_1"

    def test_assign_tier_2(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_overall_tier
        assert assign_overall_tier("high", False, "low_relevance_candidate") == "tier_2"
        assert assign_overall_tier("medium", True, "orthosteric_candidate") == "tier_2"
        assert assign_overall_tier("medium", False, "rim_candidate") == "tier_2"

    def test_assign_tier_3(self):
        from egfr_pipeline.phase2.druggability_confidence import assign_overall_tier
        assert assign_overall_tier("low", False, "low_relevance_candidate") == "tier_3"
        assert assign_overall_tier("medium", False, "low_relevance_candidate") == "tier_3"

    def test_normalize_scores_single_value(self):
        from egfr_pipeline.phase2.druggability_confidence import normalize_scores_within_state
        pockets = [
            {"pocket_id": "A_PKT01", "receptor_id": "A", "score": "0.5"},
        ]
        result = normalize_scores_within_state(pockets, "score")
        assert result["A_PKT01"] == 0.5  # single value → 0.5 neutral

    def test_normalize_scores_range(self):
        from egfr_pipeline.phase2.druggability_confidence import normalize_scores_within_state
        pockets = [
            {"pocket_id": "A_PKT01", "receptor_id": "A", "score": "0.2"},
            {"pocket_id": "A_PKT02", "receptor_id": "A", "score": "0.8"},
        ]
        result = normalize_scores_within_state(pockets, "score")
        assert abs(result["A_PKT01"] - 0.0) < 1e-6
        assert abs(result["A_PKT02"] - 1.0) < 1e-6


# ===================================================================
# TG 2.5: Cross-State Alignment
# ===================================================================

class TestCrossStateAlignment:
    """Tests for egfr_pipeline.phase2.cross_state_alignment."""

    def _setup(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)

    def test_multi_state_comparison(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_cross_state_alignment(tmp_phase2)

        comp = _load(tmp_phase2 / "candidate_pocket_cross_state_comparison.csv")
        # Should have comparisons between 3GT8_raw and 3GT8_cl38_48
        assert len(comp) > 0
        states_in_comp = set()
        for c in comp:
            states_in_comp.add(c["receptor_id_a"])
            states_in_comp.add(c["receptor_id_b"])
        assert len(states_in_comp) == 2

    def test_match_types_valid(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_cross_state_alignment(tmp_phase2)

        comp = _load(tmp_phase2 / "candidate_pocket_cross_state_comparison.csv")
        valid = {"match", "shifted", "no_match"}
        for c in comp:
            assert c["match_type"] in valid

    def test_state_classes_valid(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_cross_state_alignment(tmp_phase2)

        classes = _load(tmp_phase2 / "candidate_pocket_state_classes.csv")
        valid = {"state_robust_pocket", "state_shifted_pocket",
                 "state_specific_pocket", "uncertain_alignment"}
        for c in classes:
            assert c["state_class"] in valid

    def test_similar_pockets_match(self, tmp_phase2):
        """3GT8_raw PKT01 and 3GT8_cl38_48 PKT01 have similar centroids → should match."""
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_cross_state_alignment(tmp_phase2)

        comp = _load(tmp_phase2 / "candidate_pocket_cross_state_comparison.csv")
        # Find the comparison between the two PKT01s
        matches = [c for c in comp if c["match_type"] == "match"]
        assert len(matches) >= 1, "Similar pockets across states should match"

    def test_single_state_all_specific(self, tmp_phase2_single_state):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment

        run_pocket_merge(tmp_phase2_single_state)
        run_patch_relationship(tmp_phase2_single_state)
        run_druggability_confidence(tmp_phase2_single_state)
        run_cross_state_alignment(tmp_phase2_single_state)

        classes = _load(tmp_phase2_single_state / "candidate_pocket_state_classes.csv")
        for c in classes:
            assert c["state_class"] == "state_specific_pocket"

    def test_druggability_enrichment(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_cross_state_alignment(tmp_phase2)

        classes = _load(tmp_phase2 / "candidate_pocket_state_classes.csv")
        # Should have druggability columns carried from TG 2.4
        for c in classes:
            assert "druggability_confidence" in c
            assert "overall_druggability_tier" in c


# ===================================================================
# TG 2.5 Unit: match classification
# ===================================================================

class TestCrossStateUnit:
    """Unit tests for cross-state alignment logic."""

    def test_close_centroids_match(self):
        from egfr_pipeline.phase2.cross_state_alignment import _classify_match
        mt, _ = _classify_match(dist=5.0, jaccard=0.0, n_shared=0)
        assert mt == "match"

    def test_high_jaccard_match(self):
        from egfr_pipeline.phase2.cross_state_alignment import _classify_match
        mt, _ = _classify_match(dist=20.0, jaccard=0.30, n_shared=3)
        assert mt == "match"

    def test_shifted_centroid(self):
        from egfr_pipeline.phase2.cross_state_alignment import _classify_match
        mt, _ = _classify_match(dist=12.0, jaccard=0.05, n_shared=1)
        assert mt == "shifted"

    def test_no_match_far(self):
        from egfr_pipeline.phase2.cross_state_alignment import _classify_match
        mt, _ = _classify_match(dist=50.0, jaccard=0.0, n_shared=0)
        assert mt == "no_match"


# ===================================================================
# TG 2.6: Phase 3 Export
# ===================================================================

class TestPhase3Export:
    """Tests for egfr_pipeline.phase2.phase3_export."""

    def _setup(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)
        run_cross_state_alignment(tmp_phase2)

    def test_export_all_pockets(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_phase3_export(tmp_phase2)

        merged = _load(tmp_phase2 / "candidate_pockets.csv")
        export = _load(tmp_phase2 / "phase3_candidate_pocket_reference.csv")
        assert len(export) == len(merged)

    def test_docking_priority_valid(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_phase3_export(tmp_phase2)

        export = _load(tmp_phase2 / "phase3_candidate_pocket_reference.csv")
        valid = {"primary", "secondary", "exploratory", "skip"}
        for r in export:
            assert r["docking_priority"] in valid, f"Invalid priority: {r['docking_priority']}"

    def test_required_fields_present(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_phase3_export(tmp_phase2)

        export = _load(tmp_phase2 / "phase3_candidate_pocket_reference.csv")
        required = ["receptor_id", "pocket_id", "centroid_x", "centroid_y",
                     "centroid_z", "relationship_class", "docking_priority"]
        for r in export:
            for field in required:
                assert r.get(field, "").strip(), f"Missing {field} in {r['pocket_id']}"

    def test_pocket_ids_unique(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_phase3_export(tmp_phase2)

        export = _load(tmp_phase2 / "phase3_candidate_pocket_reference.csv")
        ids = [r["pocket_id"] for r in export]
        assert len(ids) == len(set(ids))

    def test_handoff_note_written(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_phase3_export(tmp_phase2)

        note = tmp_phase2 / "phase2_to_phase3_handoff_note.md"
        assert note.exists()
        content = note.read_text(encoding="utf-8")
        assert "Phase 3" in content
        assert "primary" in content


# ===================================================================
# TG 2.6 Unit: docking priority
# ===================================================================

class TestDockingPriorityUnit:
    """Unit tests for docking priority assignment."""

    def test_tier1_always_primary(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_1", "orthosteric_candidate", "")
        assert pri == "primary"
        pri, _ = assign_docking_priority("tier_1", "low_relevance_candidate", "")
        assert pri == "primary"

    def test_tier2_orthosteric_primary(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_2", "orthosteric_candidate", "")
        assert pri == "primary"

    def test_tier2_allosteric_secondary(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_2", "allosteric_candidate", "")
        assert pri == "secondary"

    def test_tier3_orthosteric_secondary(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_3", "orthosteric_candidate", "")
        assert pri == "secondary"

    def test_tier3_low_relevance_skip(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_3", "low_relevance_candidate", "")
        assert pri == "skip"

    def test_tier3_allosteric_exploratory(self):
        from egfr_pipeline.phase2.phase3_export import assign_docking_priority
        pri, _ = assign_docking_priority("tier_3", "allosteric_candidate", "")
        assert pri == "exploratory"


# ===================================================================
# TG 2.7: Review Report
# ===================================================================

@pytest.mark.reporting
class TestReviewReport:
    """Tests for egfr_pipeline.phase2.review_report."""

    def _setup(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)
        run_cross_state_alignment(tmp_phase2)
        run_phase3_export(tmp_phase2)

    def test_report_generated(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.review_report import run_review_report
        path = run_review_report(tmp_phase2)
        assert path.exists()
        content = path.read_text(encoding="utf-8")
        assert len(content) > 500

    def test_report_has_all_sections(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.review_report import run_review_report
        path = run_review_report(tmp_phase2)
        content = path.read_text(encoding="utf-8")

        sections = [
            "Executive Summary",
            "Pocket Proposal Sources",
            "Merged Pocket Catalog",
            "PPI Patch Relationship",
            "Druggability Confidence",
            "Cross-State Pocket Alignment",
            "Phase 3 Docking Recommendation",
            "Limitations and Caveats",
            "Phase 2 Output File Inventory",
        ]
        for section in sections:
            assert section in content, f"Missing section: {section}"

    def test_report_mentions_all_pockets(self, tmp_phase2):
        self._setup(tmp_phase2)
        from egfr_pipeline.phase2.review_report import run_review_report
        path = run_review_report(tmp_phase2)
        content = path.read_text(encoding="utf-8")

        merged = _load(tmp_phase2 / "candidate_pockets.csv")
        for p in merged:
            assert p["pocket_id"] in content


# ===================================================================
# End-to-End Integration: Full Pipeline
# ===================================================================

class TestEndToEnd:
    """Full pipeline integration test (TG 2.2 → 2.7)."""

    def test_full_pipeline_multi_state(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        from egfr_pipeline.phase2.review_report import run_review_report

        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)
        run_cross_state_alignment(tmp_phase2)
        run_phase3_export(tmp_phase2)
        run_review_report(tmp_phase2)

        # All expected files exist
        expected_files = [
            "candidate_pockets.csv",
            "candidate_pocket_merge_table.csv",
            "candidate_pocket_provenance.csv",
            "pocket_patch_relationship.csv",
            "pocket_patch_relationship_metrics.csv",
            "druggability_proposal_summary.csv",
            "candidate_pocket_support_flags.csv",
            "candidate_pocket_cross_state_comparison.csv",
            "candidate_pocket_state_classes.csv",
            "phase3_candidate_pocket_reference.csv",
            "phase2_to_phase3_handoff_note.md",
            "phase2_candidate_pocket_report.md",
        ]
        for fname in expected_files:
            assert (tmp_phase2 / fname).exists(), f"Missing: {fname}"

    def test_full_pipeline_single_state(self, tmp_phase2_single_state):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        from egfr_pipeline.phase2.phase3_export import run_phase3_export
        from egfr_pipeline.phase2.review_report import run_review_report

        d = tmp_phase2_single_state
        run_pocket_merge(d)
        run_patch_relationship(d)
        run_druggability_confidence(d)
        run_cross_state_alignment(d)
        run_phase3_export(d)
        run_review_report(d)

        report = (d / "phase2_candidate_pocket_report.md").read_text(encoding="utf-8")
        assert "state_specific_pocket" in report

    def test_idempotency(self, tmp_phase2):
        """Running pipeline twice produces identical output."""
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence

        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)

        first = (tmp_phase2 / "druggability_proposal_summary.csv").read_text(encoding="utf-8")

        # Re-run
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)

        second = (tmp_phase2 / "druggability_proposal_summary.csv").read_text(encoding="utf-8")
        assert first == second, "Pipeline should be deterministic"


# ===================================================================
# Schema cross-check
# ===================================================================

class TestSchemaConsistency:
    """Verify declared schemas match actual output."""

    def test_merged_pocket_schema(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge, MERGED_POCKET_COLUMNS
        run_pocket_merge(tmp_phase2)
        with open(tmp_phase2 / "candidate_pockets.csv") as f:
            actual = csv.DictReader(f).fieldnames
        assert list(actual) == MERGED_POCKET_COLUMNS

    def test_relationship_schema(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import (
            run_patch_relationship, RELATIONSHIP_COLUMNS, METRICS_COLUMNS
        )
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)

        with open(tmp_phase2 / "pocket_patch_relationship.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == RELATIONSHIP_COLUMNS
        with open(tmp_phase2 / "pocket_patch_relationship_metrics.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == METRICS_COLUMNS

    def test_druggability_schema(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import (
            run_druggability_confidence, SUMMARY_COLUMNS, SUPPORT_FLAGS_COLUMNS
        )
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)

        with open(tmp_phase2 / "druggability_proposal_summary.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == SUMMARY_COLUMNS
        with open(tmp_phase2 / "candidate_pocket_support_flags.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == SUPPORT_FLAGS_COLUMNS

    def test_cross_state_schema(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import (
            run_cross_state_alignment, COMPARISON_COLUMNS, STATE_CLASS_COLUMNS
        )
        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)
        run_cross_state_alignment(tmp_phase2)

        with open(tmp_phase2 / "candidate_pocket_cross_state_comparison.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == COMPARISON_COLUMNS
        with open(tmp_phase2 / "candidate_pocket_state_classes.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == STATE_CLASS_COLUMNS

    def test_export_schema(self, tmp_phase2):
        from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
        from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
        from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
        from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
        from egfr_pipeline.phase2.phase3_export import run_phase3_export, EXPORT_COLUMNS

        run_pocket_merge(tmp_phase2)
        run_patch_relationship(tmp_phase2)
        run_druggability_confidence(tmp_phase2)
        run_cross_state_alignment(tmp_phase2)
        run_phase3_export(tmp_phase2)

        with open(tmp_phase2 / "phase3_candidate_pocket_reference.csv") as f:
            assert list(csv.DictReader(f).fieldnames) == EXPORT_COLUMNS
