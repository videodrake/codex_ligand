"""Group 1 tests: experimental data integration (F-1).

Covers:
1. Ligand diversity CSV format and values (AC-1.1)
2. is_atp_site determination logic at 50% boundary (AC-1.2)
3. Ko et al. consistency check: PASS/WARNING/FAIL (AC-1.3)
4. Sheet contact columns integrity (AC-1.5)
5. Backward compatibility (new columns don't break existing parsers)
"""

import csv
import io
import logging
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ── 1. Ligand diversity CSV (AC-1.1) ─────────────────────────────────────

def test_ligand_diversity_csv_exists():
    csv_path = PROJECT_ROOT / "output" / "ligand_diversity_assessment.csv"
    assert csv_path.exists(), "ligand_diversity_assessment.csv not found"


def test_ligand_diversity_csv_format():
    csv_path = PROJECT_ROOT / "output" / "ligand_diversity_assessment.csv"
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    descriptor_rows = [r for r in rows if r["row_type"] == "descriptor"]
    pairwise_rows = [r for r in rows if r["row_type"] == "pairwise"]

    assert len(descriptor_rows) == 3
    assert len(pairwise_rows) == 3

    for pr in pairwise_rows:
        tanimoto = float(pr["tanimoto_similarity"])
        assert 0.0 <= tanimoto <= 1.0
        assert pr["diversity_judgment"] in ("diverse", "moderate", "similar")

    for dr in descriptor_rows:
        assert float(dr["MW"]) > 0
        assert dr["ligand_a"] in ("173940", "97806", "VAX-C12_0")


# ── 2. is_atp_site determination (AC-1.2) ────────────────────────────────

def test_is_atp_site_above_50_percent():
    """Pocket with >50% ATP site residues → is_atp_site=True."""
    from egfr_pipeline.vina.pocket_summary import summarize_pose_rows

    # Build fake pose rows: all residues are in ATP site (e.g., hinge region)
    poses = [
        {
            "receptor_id": "test",
            "pocket_id": "P001",
            "ligand_id": "lig1",
            "centroid_x": "0", "centroid_y": "0", "centroid_z": "0",
            "affinity": "-8.0",
            "pose_rank": "1",
            "contact_residues": "LEU788;THR790;GLN791;MET793;PRO794",
        },
    ]
    pocket_rows, _, _ = summarize_pose_rows(poses)
    assert len(pocket_rows) == 1
    # All 5 residues are in ATP site (788, 790, 791, 793, 794) → 100%
    assert pocket_rows[0]["is_atp_site"] is True


def test_is_atp_site_below_50_percent():
    """Pocket with <50% ATP site residues → is_atp_site=False."""
    from egfr_pipeline.vina.pocket_summary import summarize_pose_rows

    # 1 ATP residue (790) + 4 non-ATP residues (860, 870, 880, 890)
    poses = [
        {
            "receptor_id": "test",
            "pocket_id": "P002",
            "ligand_id": "lig1",
            "centroid_x": "0", "centroid_y": "0", "centroid_z": "0",
            "affinity": "-6.0",
            "pose_rank": "1",
            "contact_residues": "THR790;ALA860;ALA870;ALA880;ALA890",
        },
    ]
    pocket_rows, _, _ = summarize_pose_rows(poses)
    assert pocket_rows[0]["is_atp_site"] is False


def test_is_atp_site_borderline():
    """Pocket with 40-60% overlap → is_atp_site_borderline=True."""
    from egfr_pipeline.vina.pocket_summary import summarize_pose_rows

    # 2 ATP (790, 791) + 3 non-ATP (860, 870, 880) → 40%
    poses = [
        {
            "receptor_id": "test",
            "pocket_id": "P003",
            "ligand_id": "lig1",
            "centroid_x": "0", "centroid_y": "0", "centroid_z": "0",
            "affinity": "-6.0",
            "pose_rank": "1",
            "contact_residues": "THR790;GLN791;ALA860;ALA870;ALA880",
        },
    ]
    pocket_rows, _, _ = summarize_pose_rows(poses)
    assert pocket_rows[0]["is_atp_site_borderline"] is True


# ── 3. Ko et al. consistency check (AC-1.3) ──────────────────────────────

def _write_patch_csv(path, rows):
    """Helper to write a mock patch_reference CSV."""
    from egfr_pipeline.phase1.review_report import PATCH_REFERENCE_COLUMNS
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PATCH_REFERENCE_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _make_hotspot_row(resnum, is_hotspot=True):
    return {
        "residue_id": f"RES{resnum}",
        "residue_num": str(resnum),
        "residue_name": "ALA",
        "chain": "B",
        "lobe_label": "partner",
        "evidence_source": "primary",
        "construct_type": "extended_beta_meander",
        "orientation_validation_status": "orientation_validated",
        "robustness_class": "robust",
        "n_states_present": "3",
        "states_present": "3GT8_raw;EGFR_160-185;EGFR_170-200",
        "global_max_occupancy": "0.8",
        "is_hotspot_any_state": str(is_hotspot),
        "pyrosetta_evidence": "True",
        "lightdock_evidence": "False",
        "method_agreement": "pyrosetta_only",
        "confidence": "medium",
        "notes": "",
    }


def test_ko_check_pass(tmp_path):
    """5 active face residues in hotspot → PASS (no exception)."""
    from run_production import _validate_ko_consistency

    csv_path = tmp_path / "patch.csv"
    rows = [_make_hotspot_row(r) for r in [961, 962, 963, 968, 971]]
    _write_patch_csv(csv_path, rows)
    _validate_ko_consistency(csv_path)  # Should not raise


def test_ko_check_fail(tmp_path):
    """Only 2 active face residues → FAIL (ValueError)."""
    from run_production import _validate_ko_consistency

    csv_path = tmp_path / "patch.csv"
    rows = [_make_hotspot_row(r) for r in [961, 962, 977, 985]]
    _write_patch_csv(csv_path, rows)
    with pytest.raises(ValueError, match="Ko et al. FAIL"):
        _validate_ko_consistency(csv_path)


def test_ko_check_warning_sheet10(tmp_path, caplog):
    """Sheet 10/11 residues in hotspot → WARNING logged."""
    from run_production import _validate_ko_consistency

    csv_path = tmp_path / "patch.csv"
    # 5 active face + 1 sheet 10
    rows = [_make_hotspot_row(r) for r in [961, 962, 963, 968, 971, 977]]
    _write_patch_csv(csv_path, rows)
    with caplog.at_level(logging.WARNING):
        _validate_ko_consistency(csv_path)
    assert "sheet 10/11" in caplog.text.lower() or "WARNING" in caplog.text


# ── 4. Sheet contact columns (AC-1.5) ────────────────────────────────────

def test_sheet_contact_columns():
    """Pocket summary includes sheet contact columns with integer values."""
    from egfr_pipeline.vina.pocket_summary import summarize_pose_rows

    # Pose touching sheet 8 residue 962 (MYO1D side)
    poses = [
        {
            "receptor_id": "test",
            "pocket_id": "P010",
            "ligand_id": "lig1",
            "centroid_x": "0", "centroid_y": "0", "centroid_z": "0",
            "affinity": "-7.0",
            "pose_rank": "1",
            "contact_residues": "VAL962;SER964;ALA870;ALA880",
        },
    ]
    pocket_rows, _, _ = summarize_pose_rows(poses)
    p = pocket_rows[0]
    assert "contacts_sheet_8_9" in p
    assert "contacts_sheet_10_11" in p
    assert "contacts_sheet_12" in p
    assert isinstance(p["contacts_sheet_8_9"], int)
    assert p["contacts_sheet_8_9"] == 2  # residues 962 and 964 are in sheet 8
    assert p["contacts_sheet_10_11"] == 0
    assert p["contacts_sheet_12"] == 0


# ── 5. Backward compatibility ─────────────────────────────────────────────

def test_pocket_summary_backward_compat():
    """New columns don't break existing pocket_rows structure."""
    from egfr_pipeline.vina.pocket_summary import summarize_pose_rows

    poses = [
        {
            "receptor_id": "test",
            "pocket_id": "P099",
            "ligand_id": "lig1",
            "centroid_x": "1", "centroid_y": "2", "centroid_z": "3",
            "affinity": "-5.5",
            "pose_rank": "1",
            "contact_residues": "ALA860;ALA870",
        },
    ]
    pocket_rows, drug_map_rows, occ_rows = summarize_pose_rows(poses)
    p = pocket_rows[0]

    # Original columns still present
    assert "receptor_id" in p
    assert "pocket_id" in p
    assert "n_pose" in p
    assert "best_affinity" in p
    assert "union_contact_residues" in p

    # New columns present
    assert "is_atp_site" in p
    assert "contacts_sheet_8_9" in p
    assert "exclusion_reason" not in p  # exclusion_reason is in verdict, not pocket_summary
