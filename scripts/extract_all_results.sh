#!/bin/bash
# =============================================================================
# EGFR-MYO1D 파이프라인 — HPC 결과 전체 추출 스크립트
#
# 사용법:
#   cd /work4/hwang/onepack/my_second_project/codex_ligand2
#   bash scripts/extract_all_results.sh
#
# 결과:
#   results_export/ 디렉토리에 모든 CSV/MD 복사
#   results_export.tar.gz 생성 → 이 환경으로 전달
# =============================================================================

set -euo pipefail

EXPORT_DIR="results_export"
rm -rf "$EXPORT_DIR"
mkdir -p "$EXPORT_DIR"

echo "================================================"
echo "  EGFR-MYO1D 전체 결과 추출"
echo "================================================"

# ---------------------------------------------------------------------------
# 1. Workflow A — PPI 결합 잔기 (핵심)
# ---------------------------------------------------------------------------
WA="output/workflow_a"

echo ""
echo "[1/6] Workflow A: PPI postprocess"
mkdir -p "$EXPORT_DIR/wf_a/phase3_ppi_postprocess"
for f in \
    ppi_pyrosetta_residues.csv \
    ppi_pyrosetta_summary.csv \
    ppi_interface_patch_table.csv \
; do
    src="$WA/phase3_ppi_postprocess/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_a/phase3_ppi_postprocess/" && echo "  OK: $f" || echo "  SKIP: $f (not found)"
done

echo ""
echo "[2/6] Workflow A: Verdict + Report"
mkdir -p "$EXPORT_DIR/wf_a/phase5_verdict"
mkdir -p "$EXPORT_DIR/wf_a/phase4_vina_postprocess"
mkdir -p "$EXPORT_DIR/wf_a/phase6_report"
for f in valid_sites.csv cross_method_agreement.csv; do
    src="$WA/phase5_verdict/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_a/phase5_verdict/" && echo "  OK: $f" || echo "  SKIP: $f"
done
for f in vina_pocket_table.csv vina_pocket_comparison.csv vina_drug_pocket_map.csv; do
    src="$WA/phase4_vina_postprocess/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_a/phase4_vina_postprocess/" && echo "  OK: $f" || echo "  SKIP: $f"
done
# Report
for f in project_report.txt; do
    src="$WA/phase6_report/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_a/phase6_report/" && echo "  OK: $f" || echo "  SKIP: $f"
done

# ---------------------------------------------------------------------------
# 2. Workflow B Phase 1 — PPI 패치 분석 (상태별)
# ---------------------------------------------------------------------------
WB="output/workflow_b"

echo ""
echo "[3/6] Workflow B Phase 1: PPI patch analysis (per-state)"
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
    mkdir -p "$EXPORT_DIR/wf_b/phase1_ppi_analysis/$state"
    for f in \
        orientation_filter_log.csv \
        ppi_cluster_summary.csv \
        ppi_hotspot_residues.csv \
        ppi_interface_patch_table.csv \
        pyrosetta_interface_models.csv \
        pyrosetta_interface_residue_table.csv \
    ; do
        src="$WB/phase1_ppi_analysis/$state/$f"
        [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_b/phase1_ppi_analysis/$state/" && echo "  OK: $state/$f" || echo "  SKIP: $state/$f"
    done
done

# Phase 1 cross-state 핸드오프
for f in \
    phase1_downstream_patch_reference.csv \
    phase1_cross_state_comparison.csv \
    phase1_review_report.md \
    lightdock_cross_validation_summary.csv \
; do
    src="$WB/phase1_ppi_analysis/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_b/phase1_ppi_analysis/" && echo "  OK: $f" || echo "  SKIP: $f"
done

# ---------------------------------------------------------------------------
# 3. Workflow B Phase 2 — Pocket analysis (전체)
# ---------------------------------------------------------------------------
echo ""
echo "[4/6] Workflow B Phase 2: Pocket analysis"
mkdir -p "$EXPORT_DIR/wf_b/phase2_pocket_analysis"
for f in \
    candidate_pockets.csv \
    candidate_pockets_raw.csv \
    candidate_pocket_provenance.csv \
    pocket_patch_relationship.csv \
    pocket_patch_relationship_metrics.csv \
    druggability_proposal_summary.csv \
    candidate_pocket_support_flags.csv \
    candidate_pocket_state_classes.csv \
    candidate_pocket_cross_state_comparison.csv \
    phase3_candidate_pocket_reference.csv \
    phase2_patch_reference_normalized.csv \
    phase2_candidate_pocket_report.md \
    phase2_relationship_classification_note.md \
    phase2_to_phase3_handoff_note.md \
    candidate_pocket_source_summary.csv \
    candidate_pocket_merge_table.csv \
; do
    src="$WB/phase2_pocket_analysis/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_b/phase2_pocket_analysis/" && echo "  OK: $f" || echo "  SKIP: $f"
done

# ---------------------------------------------------------------------------
# 4. Workflow B Phase 3 — Focused Vina docking
# ---------------------------------------------------------------------------
echo ""
echo "[5/6] Workflow B Phase 3: Focused Vina docking"
mkdir -p "$EXPORT_DIR/wf_b/phase3_focused_docking"
for f in \
    phase3_round_log.csv \
    phase4_docking_evidence_reference.csv \
    phase3_docking_summary.md \
    phase3_review_report.md \
; do
    src="$WB/phase3_focused_docking/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_b/phase3_focused_docking/" && echo "  OK: $f" || echo "  SKIP: $f"
done

# ---------------------------------------------------------------------------
# 5. Workflow B Phase 4 — Perturbation scoring (최종)
# ---------------------------------------------------------------------------
echo ""
echo "[6/6] Workflow B Phase 4: Perturbation scoring (final)"
mkdir -p "$EXPORT_DIR/wf_b/phase4_scoring"
for f in \
    phase4_final_review_table.csv \
    phase4_expanded_evidence_table.csv \
    phase4_axis_definition_table.csv \
    phase4_axis_scores.csv \
    phase4_state_interpretation.csv \
    phase4_evidence_normalized.csv \
    phase4_evidence_validation.md \
    phase4_score_framework.md \
    final_candidate_classes.csv \
    phase4_mechanistic_classification_note.md \
    perturbation_axis_scores.csv \
    perturbation_candidate_table.csv \
    phase4_ranking_method_note.md \
    phase4_accessibility_note.md \
    phase4_review_expanded.csv \
    phase4_final_report.md \
    integrated_phase4_report.md \
; do
    src="$WB/phase4_scoring/$f"
    [ -f "$src" ] && cp "$src" "$EXPORT_DIR/wf_b/phase4_scoring/" && echo "  OK: $f" || echo "  SKIP: $f"
done

# ---------------------------------------------------------------------------
# 6. 파일 목록 + 통계
# ---------------------------------------------------------------------------
echo ""
echo "================================================"
echo "  추출 완료 — 파일 목록"
echo "================================================"
find "$EXPORT_DIR" -type f | sort | while read -r fpath; do
    size=$(du -h "$fpath" | cut -f1)
    echo "  $size  $fpath"
done

TOTAL_FILES=$(find "$EXPORT_DIR" -type f | wc -l)
TOTAL_SIZE=$(du -sh "$EXPORT_DIR" | cut -f1)
echo ""
echo "  총 파일: ${TOTAL_FILES}개"
echo "  총 크기: ${TOTAL_SIZE}"

# tar.gz 생성
tar czf results_export.tar.gz "$EXPORT_DIR"
TAR_SIZE=$(du -sh results_export.tar.gz | cut -f1)
echo ""
echo "  아카이브: results_export.tar.gz ($TAR_SIZE)"
echo ""
echo "================================================"
echo "  이 파일을 Claude Code 환경으로 전달하세요:"
echo "  results_export.tar.gz"
echo "================================================"
