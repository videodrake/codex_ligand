"""Group 2 tests: Vina blind docking bias quantification (F-2).

Covers:
1. Pose distribution CSV format and region fraction sum (AC-2.1)
2. Residue mapping comparison mismatch detection (AC-2.2)
3. Bootstrap-Verdict backward compatibility (AC-2.3)
4. Existing valid_sites.csv — no unintended verdict flips (AC-2.3)
5. Region classification consistency (same residue → same region) (AC-2.1)
6. Affinity distribution analysis output (AC-2.4)
7. Vina scoring bias documentation (AC-2.5)
"""

import csv
from pathlib import Path
from typing import Dict, List

import pytest


PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ── 1. Pose distribution CSV format & fraction sum (AC-2.1) ───────────


def test_pose_region_classifier_fraction_sums_to_one():
    """Region fractions per (receptor_id, ligand_id) should sum to ~1.0."""
    from egfr_pipeline.vina.pose_region_classifier import classify_poses_by_region

    poses = [
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "LEU788;THR790;GLN791"},
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "ALA860;ALA870;ALA880"},
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "ALA699;ALA700;ALA701"},
    ]
    dist_rows, warnings = classify_poses_by_region(poses)
    assert len(dist_rows) > 0

    total_fraction = sum(float(r["fraction"]) for r in dist_rows)
    assert abs(total_fraction - 1.0) < 0.01, f"Fractions sum to {total_fraction}, expected ~1.0"


def test_pose_region_classifier_csv_fields():
    """Output rows have the required DISTRIBUTION_FIELDS."""
    from egfr_pipeline.vina.pose_region_classifier import (
        DISTRIBUTION_FIELDS,
        classify_poses_by_region,
    )

    poses = [
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "ALA860;ALA870;ALA880"},
    ]
    dist_rows, _ = classify_poses_by_region(poses)
    assert len(dist_rows) >= 1
    for field in DISTRIBUTION_FIELDS:
        assert field in dist_rows[0], f"Missing field: {field}"


def test_pose_region_classifier_c_lobe_surface_warning():
    """C-lobe surface < 10% triggers WARNING (AC-2.1)."""
    from egfr_pipeline.vina.pose_region_classifier import classify_poses_by_region

    # All 10 poses go to ATP site; 0 go to c_lobe_surface → < 10%
    poses = [
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "LEU788;THR790;GLN791;MET793;PRO794"}
        for _ in range(10)
    ]
    _, warnings = classify_poses_by_region(poses)
    assert any("c_lobe_surface" in w.lower() or "c-lobe" in w.lower() for w in warnings), (
        f"Expected C-lobe surface warning, got: {warnings}"
    )


def test_pose_region_classifier_multiple_receptor_ligand():
    """Fractions are computed per (receptor_id, ligand_id), not globally."""
    from egfr_pipeline.vina.pose_region_classifier import classify_poses_by_region

    poses = [
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "ALA860;ALA870;ALA880"},
        {"receptor_id": "R2", "ligand_id": "L1", "contact_residues": "LEU788;THR790;GLN791"},
    ]
    dist_rows, _ = classify_poses_by_region(poses)

    # Each (receptor, ligand) group should have fractions summing to 1.0
    from collections import defaultdict
    frac_sums: Dict[tuple, float] = defaultdict(float)
    for row in dist_rows:
        key = (row["receptor_id"], row["ligand_id"])
        frac_sums[key] += float(row["fraction"])

    for key, total in frac_sums.items():
        assert abs(total - 1.0) < 0.01, f"{key} fractions sum to {total}"


# ── 2. Residue mapping comparison mismatch detection (AC-2.2) ─────────


def test_parse_pdb_residue_identity(tmp_path):
    """parse_pdb_residue_identity extracts (resnum, resname) from PDB."""
    from egfr_pipeline.validate import parse_pdb_residue_identity

    pdb_content = (
        "ATOM      1  CA  ALA A 699      10.000 10.000 10.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A 700      12.000 10.000 10.000  1.00 20.00           C\n"
        "ATOM      3  CA  VAL A 701      14.000 10.000 10.000  1.00 20.00           C\n"
        "END\n"
    )
    pdb_path = tmp_path / "test.pdb"
    pdb_path.write_text(pdb_content)

    result = parse_pdb_residue_identity(pdb_path)
    assert "A" in result
    assert result["A"][699] == "ALA"
    assert result["A"][700] == "GLY"
    assert result["A"][701] == "VAL"


def test_residue_numbering_mismatch_detected(tmp_path):
    """Deliberately changed amino acid at a residue → WARNING or FAIL."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    # Two PDBs: same residue numbers but one has ALA→GLY at residue 860
    pdb_a = tmp_path / "state_a.pdb"
    pdb_b = tmp_path / "state_b.pdb"

    lines_a = []
    lines_b = []
    for idx, resnum in enumerate([860, 870, 880, 890], 1):
        lines_a.append(
            f"ATOM  {idx:5d}  CA  ALA A{resnum:4d}    "
            f"{10.0 + idx * 2:8.3f}{10.000:8.3f}{10.000:8.3f}  1.00 20.00           C"
        )
        resname = "GLY" if resnum == 860 else "ALA"
        lines_b.append(
            f"ATOM  {idx:5d}  CA  {resname} A{resnum:4d}    "
            f"{10.0 + idx * 2:8.3f}{10.000:8.3f}{10.000:8.3f}  1.00 20.00           C"
        )

    pdb_a.write_text("\n".join(lines_a) + "\nEND\n")
    pdb_b.write_text("\n".join(lines_b) + "\nEND\n")

    config = {
        "receptors": [
            {"id": "state_a", "pdb": str(pdb_a)},
            {"id": "state_b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    alignment_rows: List[dict] = []
    check_residue_numbering(config, result, alignment_rows=alignment_rows)

    # Should detect at least one mismatch at residue 860
    assert len(alignment_rows) >= 1 or result.warnings or result.failures, (
        "Expected mismatch detection for residue 860 ALA→GLY"
    )


def test_residue_numbering_no_mismatch(tmp_path):
    """Identical PDBs → no mismatches, exit_code=0."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    pdb_content = (
        "ATOM      1  CA  ALA A 860      10.000 10.000 10.000  1.00 20.00           C\n"
        "ATOM      2  CA  ALA A 870      12.000 10.000 10.000  1.00 20.00           C\n"
        "END\n"
    )
    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    pdb_a.write_text(pdb_content)
    pdb_b.write_text(pdb_content)

    config = {
        "receptors": [
            {"id": "a", "pdb": str(pdb_a)},
            {"id": "b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    check_residue_numbering(config, result)
    assert result.exit_code == 0


# ── 3. Bootstrap-Verdict backward compatibility (AC-2.3) ──────────────


def test_bootstrap_confidence_categories():
    """_bootstrap_confidence maps frac to correct categories."""
    from egfr_pipeline.verdict import _bootstrap_confidence

    assert _bootstrap_confidence(0.90) == "high"
    assert _bootstrap_confidence(0.80) == "high"
    assert _bootstrap_confidence(0.79) == "medium"
    assert _bootstrap_confidence(0.50) == "medium"
    assert _bootstrap_confidence(0.49) == "low"
    assert _bootstrap_confidence(0.0) == "low"
    assert _bootstrap_confidence("") == "not_assessed"
    assert _bootstrap_confidence(None) == "not_assessed"
    assert _bootstrap_confidence("invalid") == "not_assessed"


def test_bootstrap_not_run_backward_compat():
    """Without bootstrap data, score_pocket uses existing logic (no crash)."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-7.5",
        "n_pose": "10",
        "n_ligand": "2",
        # No pocket_exists_frac — simulates bootstrap not run
        "dominant_ligand_fraction": "0.60",
        "ligand_pose_entropy": "0.90",
    }

    total, verdict, reasons, vina_score, ppi_score, cross_score, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert verdict in {"STRONG", "MODERATE", "WEAK"}
    assert total > 0
    assert vina_score > 0


def test_bootstrap_low_frac_zeros_stability_pts():
    """pocket_exists_frac < stability_min → stability_pts = 0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-8.0",
        "n_pose": "15",
        "n_ligand": "3",
        "pocket_exists_frac": "0.20",
        "dominant_ligand_fraction": "0.55",
        "ligand_pose_entropy": "1.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 0.0


def test_bootstrap_high_frac_full_stability_pts():
    """pocket_exists_frac >= stability_great → stability_pts = 10.0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-8.0",
        "n_pose": "15",
        "n_ligand": "3",
        "pocket_exists_frac": "0.85",
        "dominant_ligand_fraction": "0.55",
        "ligand_pose_entropy": "1.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 10.0


def test_bootstrap_medium_frac_dampens_convergence():
    """pocket_exists_frac < stability_good → convergence_pts *= 0.5."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-8.5",
        "n_pose": "12",
        "n_ligand": "2",
        "pocket_exists_frac": "0.45",  # < 0.60 (stability_good)
        "dominant_ligand_fraction": "0.95",
        "ligand_pose_entropy": "0.20",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    # Convergence should be halved (damped)
    assert raw["vina_convergence_pts"] == 7.5  # 15 * 0.5


# ── 4. No unintended verdict flips (AC-2.3) ──────────────────────────
# This is validated by test_generate_verdict_prefers_stable_recurrent_site
# in test_verdict.py. Here we add a complementary check: generating a
# verdict WITH bootstrap data should not crash, and pocket ordering
# should respect bootstrap stability.

def test_verdict_stable_pocket_outranks_unstable(tmp_path):
    """Pocket with high bootstrap frac should score higher than low frac."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    base = {
        "best_affinity": "-7.5",
        "n_pose": "10",
        "n_ligand": "2",
        "dominant_ligand_fraction": "0.60",
        "ligand_pose_entropy": "0.90",
    }

    stable = dict(base, pocket_exists_frac="0.90")
    unstable = dict(base, pocket_exists_frac="0.20")

    def _score(pocket):
        total, *_ = score_pocket(
            pocket,
            ppi_agreement=None,
            cross_receptor_matches=[],
            cross_support=None,
            thresholds=dict(DEFAULT_THRESHOLDS),
            has_ppi_data=False,
        )
        return total

    assert _score(stable) > _score(unstable)


# ── 5. Region classification consistency (AC-2.1) ────────────────────


def test_region_classification_consistency():
    """Same residue always maps to the same region via get_region()."""
    from egfr_pipeline.region_definitions import get_region

    # Check a set of known residues
    known_mappings = {
        788: "atp_site",   # hinge
        790: "atp_site",   # hinge
        860: "c_lobe_surface",
        700: "n_lobe",
        835: "atp_site",   # catalytic
    }
    for resnum, expected in known_mappings.items():
        result = get_region(resnum)
        assert result == expected, f"get_region({resnum}) = {result}, expected {expected}"
        # Calling again should give the same result (consistency)
        assert get_region(resnum) == result


def test_region_classification_outside_kinase_returns_none():
    """Residues outside kinase domain (699-1007) → None."""
    from egfr_pipeline.region_definitions import get_region

    assert get_region(500) is None
    assert get_region(1100) is None
    assert get_region(0) is None


def test_pose_classifier_consistent_with_region_definitions():
    """_classify_single_pose uses get_region() and returns consistent results."""
    from egfr_pipeline.vina.pose_region_classifier import _classify_single_pose

    # All ATP site residues → "atp_site"
    assert _classify_single_pose("LEU788;THR790;GLN791;MET793") == "atp_site"

    # All c_lobe_surface residues → "c_lobe_surface"
    assert _classify_single_pose("ALA860;ALA870;ALA880") == "c_lobe_surface"

    # Mixed: no majority → "mixed"
    assert _classify_single_pose("LEU788;ALA860") == "mixed"

    # No classifiable residues → "unknown"
    assert _classify_single_pose("ALA1200;ALA1300") == "unknown"

    # Empty → "unknown"
    assert _classify_single_pose("") == "unknown"


# ── 6. Affinity distribution analysis (AC-2.4) ───────────────────────


def test_affinity_distribution_analysis_output():
    """Analyze affinity distribution returns required fields."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    pocket_rows = [
        {"best_affinity": "-9.0", "is_atp_site": "False"},
        {"best_affinity": "-7.5", "is_atp_site": "False"},
        {"best_affinity": "-6.0", "is_atp_site": "False"},
        {"best_affinity": "-5.5", "is_atp_site": "False"},
        {"best_affinity": "-8.5", "is_atp_site": "False"},
    ]

    result = analyze_affinity_distribution(pocket_rows)

    assert result["n_pockets"] == 5
    assert "p25" in result
    assert "p50_median" in result
    assert "p75" in result
    assert "p90" in result
    assert result["discriminative_power"] in ("sufficient", "insufficient")
    assert result["recommendation"] is not None


def test_affinity_distribution_excludes_atp_site():
    """ATP site pockets are excluded from the analysis."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    pocket_rows = [
        {"best_affinity": "-11.0", "is_atp_site": "True"},
        {"best_affinity": "-10.5", "is_atp_site": "True"},
        {"best_affinity": "-6.5", "is_atp_site": "False"},
    ]

    result = analyze_affinity_distribution(pocket_rows)
    assert result["n_pockets"] == 1  # Only the non-ATP pocket


def test_affinity_distribution_empty_input():
    """No valid pockets → no_data."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    result = analyze_affinity_distribution([])
    assert result["n_pockets"] == 0
    assert result["discriminative_power"] == "no_data"


def test_affinity_distribution_discriminative_sufficient():
    """Spread distribution → sufficient discriminative power."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    # Create pockets with balanced distribution across brackets
    pocket_rows = (
        [{"best_affinity": str(-9.0 + i * 0.1), "is_atp_site": "False"} for i in range(10)]
        + [{"best_affinity": str(-7.0 + i * 0.1), "is_atp_site": "False"} for i in range(10)]
        + [{"best_affinity": str(-5.5 + i * 0.1), "is_atp_site": "False"} for i in range(10)]
    )

    result = analyze_affinity_distribution(pocket_rows)
    assert result["discriminative_power"] == "sufficient"


def test_affinity_distribution_discriminative_insufficient():
    """Concentrated distribution → insufficient discriminative power."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    # 80% of pockets in weak bracket (> -6.5)
    pocket_rows = [
        {"best_affinity": str(-5.0 + i * 0.05), "is_atp_site": "False"} for i in range(20)
    ] + [
        {"best_affinity": "-9.0", "is_atp_site": "False"},
        {"best_affinity": "-7.0", "is_atp_site": "False"},
    ]

    result = analyze_affinity_distribution(pocket_rows)
    assert result["discriminative_power"] == "insufficient"


# ── 7. Vina scoring bias documentation (AC-2.5) ──────────────────────


def test_vina_bias_documentation_exists():
    """결과 해석 가이드에 소수성 과대평가/hydrophobic bias 문구 존재."""
    guide_path = PROJECT_ROOT / "docs" / "manual_vina.md"
    assert guide_path.exists(), "manual_vina.md not found"

    content = guide_path.read_text(encoding="utf-8")
    assert "소수성 과대평가" in content, "Missing '소수성 과대평가' in interpretation guide"
    assert "C-lobe" in content and ("surface" in content or "Surface" in content), (
        "Missing C-lobe surface affinity guidance"
    )


def test_vina_bias_doc_clobe_affinity_guidance():
    """C-lobe surface 포켓의 -5~-7 kcal/mol 해석 지침 존재."""
    guide_path = PROJECT_ROOT / "docs" / "manual_vina.md"
    content = guide_path.read_text(encoding="utf-8")
    assert "-5" in content and "-7" in content, (
        "Missing C-lobe surface affinity range guidance (-5 ~ -7)"
    )
    assert "의미 있는 결합" in content or "meaningful binding" in content.lower(), (
        "Missing explanation that -5~-7 is meaningful for surface pockets"
    )


# ══════════════════════════════════════════════════════════════════════
# Edge Case tests — PRD EC-2.1, EC-2.2, EC-2.3 + boundary values
# ══════════════════════════════════════════════════════════════════════


# ── EC-2.1: C-lobe surface 0 poses → specific WARNING text ───────────


def test_ec21_c_lobe_surface_zero_poses_warning():
    """EC-2.1: 0 C-lobe surface poses → WARNING mentioning Workflow B."""
    from egfr_pipeline.vina.pose_region_classifier import classify_poses_by_region

    # All poses go to n_lobe — zero in c_lobe_surface
    poses = [
        {"receptor_id": "R1", "ligand_id": "L1", "contact_residues": "ALA699;ALA700;ALA701"}
        for _ in range(5)
    ]
    _, warnings = classify_poses_by_region(poses)
    assert len(warnings) >= 1, "Expected at least one warning for 0 C-lobe surface poses"
    zero_warning = [w for w in warnings if "no poses" in w.lower() or "0%" in w.lower() or "미도달" in w.lower()]
    assert len(zero_warning) >= 1, f"Expected zero-pose specific warning, got: {warnings}"
    # Should mention Workflow B fallback
    assert any("workflow b" in w.lower() or "focused docking" in w.lower() for w in warnings), (
        f"Warning should mention Workflow B/focused docking fallback, got: {warnings}"
    )


# ── EC-2.2: >= 10 unknown mismatches → FAIL ──────────────────────────


def _make_pdb_with_residues(path, residues, chain="A"):
    """Helper: create minimal PDB with specified (resnum, resname) pairs."""
    lines = []
    for idx, (resnum, resname) in enumerate(residues, 1):
        lines.append(
            f"ATOM  {idx:5d}  CA  {resname:3s} {chain}{resnum:4d}    "
            f"{10.0 + idx * 2:8.3f}{10.000:8.3f}{10.000:8.3f}  1.00 20.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def test_ec22_ten_or_more_mismatches_fail(tmp_path):
    """EC-2.2: >= 10 unknown identity mismatches → FAIL (exit_code=2)."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    # 15 residues with 12 mismatches (use c_lobe region residues to avoid known mutations)
    resnums = list(range(860, 875))  # 15 residues
    residues_a = [(r, "ALA") for r in resnums]
    # Change 12 of them to GLY
    residues_b = [(r, "GLY" if i < 12 else "ALA") for i, r in enumerate(resnums)]

    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    _make_pdb_with_residues(pdb_a, residues_a)
    _make_pdb_with_residues(pdb_b, residues_b)

    config = {
        "receptors": [
            {"id": "a", "pdb": str(pdb_a)},
            {"id": "b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    alignment_rows: List[dict] = []
    check_residue_numbering(config, result, alignment_rows=alignment_rows)

    assert result.exit_code == 2, (
        f"Expected FAIL (exit_code=2) for >= 10 mismatches, got exit_code={result.exit_code}"
    )
    assert len([r for r in alignment_rows if r["status"] == "unexpected"]) >= 10


def test_ec22_fewer_than_ten_mismatches_warning(tmp_path):
    """EC-2.2: < 10 unknown mismatches → WARNING (exit_code=1), not FAIL."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    resnums = list(range(860, 870))  # 10 residues
    residues_a = [(r, "ALA") for r in resnums]
    # Change 5 to GLY
    residues_b = [(r, "GLY" if i < 5 else "ALA") for i, r in enumerate(resnums)]

    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    _make_pdb_with_residues(pdb_a, residues_a)
    _make_pdb_with_residues(pdb_b, residues_b)

    config = {
        "receptors": [
            {"id": "a", "pdb": str(pdb_a)},
            {"id": "b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    check_residue_numbering(config, result)

    assert result.exit_code == 1, (
        f"Expected WARNING (exit_code=1) for < 10 mismatches, got exit_code={result.exit_code}"
    )


def test_ec22_alignment_rows_have_region_info(tmp_path):
    """AC-2.2: alignment_rows include region classification from get_region."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    residues_a = [(860, "ALA"), (870, "ALA")]
    residues_b = [(860, "GLY"), (870, "ALA")]

    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    _make_pdb_with_residues(pdb_a, residues_a)
    _make_pdb_with_residues(pdb_b, residues_b)

    config = {
        "receptors": [
            {"id": "a", "pdb": str(pdb_a)},
            {"id": "b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    alignment_rows: List[dict] = []
    check_residue_numbering(config, result, alignment_rows=alignment_rows)

    assert len(alignment_rows) >= 1
    row = alignment_rows[0]
    assert "region" in row
    assert row["region"] in ("n_lobe", "atp_site", "c_lobe_surface", "c_lobe_core", "outside_kinase")
    assert "status" in row
    assert row["status"] in ("known_mutation", "unexpected")


# ── EC-2.3: All bootstrap > 0.8 → no existing result changes ─────────


def test_ec23_all_high_bootstrap_no_convergence_damping():
    """EC-2.3: When all pockets have frac > 0.8, convergence is NOT damped."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-8.5",
        "n_pose": "12",
        "n_ligand": "2",
        "pocket_exists_frac": "0.90",  # > 0.80 → stability_great
        "dominant_ligand_fraction": "0.95",
        "ligand_pose_entropy": "0.20",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    # With frac >= stability_good (0.60), convergence is NOT halved by stability
    # base_convergence = 15.0 (n_pose=12 >= n_pose_great=10)
    # n_ligand=2, so single-ligand dominance damping (n_ligand <= 1) does NOT apply
    # stability factor = 1.0 (frac >= stability_good)
    assert raw["vina_convergence_pts"] == 15.0
    assert raw["vina_stability_pts"] == 10.0


# ── Stability boundary values ────────────────────────────────────────


def test_stability_boundary_exactly_at_min():
    """pocket_exists_frac == stability_min (0.40) → stability_pts = 3.0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-7.0",
        "n_pose": "5",
        "n_ligand": "1",
        "pocket_exists_frac": "0.40",  # exactly stability_min
        "dominant_ligand_fraction": "1.0",
        "ligand_pose_entropy": "0.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 3.0


def test_stability_boundary_exactly_at_good():
    """pocket_exists_frac == stability_good (0.60) → stability_pts = 6.0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-7.0",
        "n_pose": "5",
        "n_ligand": "1",
        "pocket_exists_frac": "0.60",  # exactly stability_good
        "dominant_ligand_fraction": "1.0",
        "ligand_pose_entropy": "0.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 6.0


def test_stability_boundary_exactly_at_great():
    """pocket_exists_frac == stability_great (0.80) → stability_pts = 10.0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-7.0",
        "n_pose": "5",
        "n_ligand": "1",
        "pocket_exists_frac": "0.80",  # exactly stability_great
        "dominant_ligand_fraction": "1.0",
        "ligand_pose_entropy": "0.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 10.0


def test_stability_boundary_just_below_min():
    """pocket_exists_frac = 0.39 (< stability_min 0.40) → stability_pts = 0."""
    from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, score_pocket

    pocket = {
        "best_affinity": "-7.0",
        "n_pose": "5",
        "n_ligand": "1",
        "pocket_exists_frac": "0.39",
        "dominant_ligand_fraction": "1.0",
        "ligand_pose_entropy": "0.0",
    }

    _, _, _, _, _, _, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert raw["vina_stability_pts"] == 0.0


# ── Malformed / edge inputs ──────────────────────────────────────────


def test_classify_poses_empty_input():
    """Empty pose list → no distribution rows, no crash."""
    from egfr_pipeline.vina.pose_region_classifier import classify_poses_by_region

    dist_rows, warnings = classify_poses_by_region([])
    assert dist_rows == []
    assert warnings == []


def test_classify_poses_malformed_contact_residues():
    """Malformed contact_residues → gracefully classified as unknown."""
    from egfr_pipeline.vina.pose_region_classifier import _classify_single_pose

    # Various malformed inputs
    assert _classify_single_pose(";;;") == "unknown"
    assert _classify_single_pose("GARBAGE_DATA") == "unknown"
    # Note: "X:???999" extracts resnum 999 which IS in kinase domain → valid classification
    # Use truly out-of-range residue numbers for "unknown"
    assert _classify_single_pose("X:???9999") == "unknown"


def test_bootstrap_confidence_string_float_input():
    """_bootstrap_confidence handles string-encoded floats (from CSV reading)."""
    from egfr_pipeline.verdict import _bootstrap_confidence

    assert _bootstrap_confidence("0.85") == "high"
    assert _bootstrap_confidence("0.50") == "medium"
    assert _bootstrap_confidence("0.30") == "low"
    assert _bootstrap_confidence("0.0") == "low"


def test_affinity_distribution_all_atp_site():
    """All pockets are ATP site → no_data (all filtered out)."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    pocket_rows = [
        {"best_affinity": "-11.0", "is_atp_site": "True"},
        {"best_affinity": "-10.0", "is_atp_site": "1"},
        {"best_affinity": "-9.5", "is_atp_site": "yes"},
    ]

    result = analyze_affinity_distribution(pocket_rows)
    assert result["n_pockets"] == 0
    assert result["discriminative_power"] == "no_data"


def test_affinity_distribution_missing_affinity_field():
    """Pockets with missing or empty best_affinity are skipped."""
    from scripts.analyze_affinity_distribution import analyze_affinity_distribution

    pocket_rows = [
        {"best_affinity": "", "is_atp_site": "False"},
        {"best_affinity": None, "is_atp_site": "False"},
        {"is_atp_site": "False"},  # no best_affinity key at all
        {"best_affinity": "-7.0", "is_atp_site": "False"},
    ]

    result = analyze_affinity_distribution(pocket_rows)
    assert result["n_pockets"] == 1  # Only the last one is valid


def test_single_receptor_residue_check_warns(tmp_path):
    """Only 1 receptor PDB → WARNING (need at least 2 for comparison)."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    pdb_path = tmp_path / "only.pdb"
    _make_pdb_with_residues(pdb_path, [(860, "ALA")])

    config = {"receptors": [{"id": "only", "pdb": str(pdb_path)}]}
    result = ValidationResult()
    check_residue_numbering(config, result)

    assert any("at least 2" in w for w in result.warnings)


def test_charmm_his_variants_not_mismatched(tmp_path):
    """CHARMM HIS variants (HSD/HSE/HSP) normalized to HIS → no mismatch."""
    from egfr_pipeline.validate import ValidationResult, check_residue_numbering

    residues_a = [(860, "HIS"), (870, "ALA")]
    residues_b = [(860, "HSD"), (870, "ALA")]

    pdb_a = tmp_path / "a.pdb"
    pdb_b = tmp_path / "b.pdb"
    _make_pdb_with_residues(pdb_a, residues_a)
    _make_pdb_with_residues(pdb_b, residues_b)

    config = {
        "receptors": [
            {"id": "a", "pdb": str(pdb_a)},
            {"id": "b", "pdb": str(pdb_b)},
        ]
    }
    result = ValidationResult()
    alignment_rows: List[dict] = []
    check_residue_numbering(config, result, alignment_rows=alignment_rows)

    # HIS vs HSD should NOT be counted as a mismatch
    unexpected = [r for r in alignment_rows if r["status"] == "unexpected"]
    assert len(unexpected) == 0, f"HIS/HSD mismatch should be normalized, got: {unexpected}"
