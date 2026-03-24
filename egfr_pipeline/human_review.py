"""Production 완료 후 context-aware human review checklist 생성 모듈."""

import csv
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _read_csv_rows(path: Path) -> List[dict]:
    """CSV를 DictReader로 읽어서 행 리스트 반환. 파일 없으면 빈 리스트."""
    if not path.exists():
        return []
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def _read_text(path: Path) -> str:
    """텍스트 파일 읽기. 없으면 빈 문자열."""
    if not path.exists():
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def _parse_contact_residue_nums(union_contact_residues: str) -> set:
    """'A:LEU844;X:ARG841;...' → {844, 841, ...}"""
    nums = set()
    if not union_contact_residues:
        return nums
    for token in re.split(r"[;,\s]+", union_contact_residues.strip()):
        token = token.strip()
        if not token:
            continue
        # "A:LEU844" → 844
        m = re.search(r"(\d+)\s*$", token)
        if m:
            nums.add(int(m.group(1)))
    return nums


def _get_atp_site_residues() -> set:
    """ATP_SITE_RESIDUES를 region_definitions에서 가져오거나 fallback."""
    try:
        from egfr_pipeline.region_definitions import ATP_SITE_RESIDUES
        return set(ATP_SITE_RESIDUES)
    except ImportError:
        return {
            719, 720, 721, 722, 723, 743, 744, 745, 761, 762, 763, 764, 765,
            766, 767, 768, 788, 789, 790, 791, 792, 793, 794, 795, 796,
            831, 832, 833, 834, 835, 836, 837, 854, 855, 856, 857, 858,
        }


def _safe_float(val, default=0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _safe_int(val, default=0) -> int:
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Checklist item data class
# ---------------------------------------------------------------------------

class _Item:
    __slots__ = ("id", "priority", "title", "context", "action", "pymol_hint")

    def __init__(self, id: str, priority: str, title: str,
                 context: str, action: str, pymol_hint: str = ""):
        self.id = id
        self.priority = priority
        self.title = title
        self.context = context
        self.action = action
        self.pymol_hint = pymol_hint


# ---------------------------------------------------------------------------
# Rule implementations
# ---------------------------------------------------------------------------

def _rule1_strong_visual(
    valid_sites: List[dict],
    pocket_table: List[dict],
) -> List[_Item]:
    """규칙 1: 모든 STRONG 포켓에 대해 시각적 검증 (MUST)."""
    # pocket_table에서 receptor+pocket별 top residue 추출
    pocket_residues: Dict[Tuple[str, str], str] = {}
    for row in pocket_table:
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        union = row.get("union_contact_residues", "").strip()
        if rid and pid and union:
            # top 5 residue numbers
            nums = sorted(_parse_contact_residue_nums(union))
            pocket_residues[(rid, pid)] = "+".join(str(n) for n in nums[:5])

    items = []
    for row in valid_sites:
        if row.get("verdict", "").strip() != "STRONG":
            continue
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        score = row.get("confidence_score", "")
        vina_q = row.get("vina_quality_score", "")
        ppi_p = row.get("ppi_proximity_score", "")
        cross_r = row.get("cross_receptor_score", "")
        evidence = row.get("evidence_profile", "")

        resi_sel = pocket_residues.get((rid, pid), "")
        pymol = ""
        if resi_sel:
            pymol = f"select pocket_{pid}, receptor and resi {resi_sel}"

        items.append(_Item(
            id=f"H-5.1-{rid}_{pid}",
            priority="MUST",
            title=f"{rid} {pid} 시각적 검증",
            context=(
                f"Score {score}/100, Vina {vina_q}, PPI {ppi_p}, Cross {cross_r}\n"
                f"         Evidence: {evidence}"
            ),
            action="PyMOL에서 이 포켓의 pose 분포가 접근 가능한 표면 cavity에 있는지 확인",
            pymol_hint=pymol,
        ))
    return items


def _rule2_atp_overlap(
    valid_sites: List[dict],
    pocket_table: List[dict],
) -> List[_Item]:
    """규칙 2: ATP site 겹침 경고 (MUST)."""
    atp_residues = _get_atp_site_residues()

    # pocket_table → (rid, pid) → contact residue nums
    pocket_contacts: Dict[Tuple[str, str], set] = {}
    for row in pocket_table:
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        union = row.get("union_contact_residues", "").strip()
        if rid and pid:
            pocket_contacts[(rid, pid)] = _parse_contact_residue_nums(union)

    items = []
    for row in valid_sites:
        if row.get("verdict", "").strip() != "STRONG":
            continue
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        contact_nums = pocket_contacts.get((rid, pid), set())
        overlap = contact_nums & atp_residues
        if overlap:
            overlap_sorted = sorted(overlap)
            items.append(_Item(
                id=f"H-4.1-{rid}_{pid}",
                priority="MUST",
                title=f"ATP site 겹침 확인: {len(overlap)}개 잔기",
                context=f"겹치는 잔기: {overlap_sorted}",
                action=(
                    "PyMOL에서 이 잔기들이 포켓 중심인지 가장자리인지 확인. "
                    "포켓 중심이 ATP site이면 false positive 가능."
                ),
            ))
    return items


def _rule3_ppi_moderate(
    valid_sites: List[dict],
    agreement: List[dict],
) -> List[_Item]:
    """규칙 3: PPI proximity가 moderate인 STRONG 포켓 (SHOULD)."""
    strong_set = set()
    for row in valid_sites:
        if row.get("verdict", "").strip() == "STRONG":
            strong_set.add((
                row.get("receptor_id", "").strip(),
                row.get("pocket_id", "").strip(),
            ))

    items = []
    for row in agreement:
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        if (rid, pid) not in strong_set:
            continue
        if row.get("spatial_proximity", "").strip() != "moderate":
            continue
        dist = row.get("spatial_dist_A", "")
        n_shared = row.get("n_shared_residues", "")
        level = row.get("agreement_level", "")
        items.append(_Item(
            id=f"H-5.2-{rid}_{pid}",
            priority="SHOULD",
            title=f"PPI proximity {dist}A 구조적 해석",
            context=f"Shared residues: {n_shared}개, Agreement: {level}",
            action=(
                "포켓과 PPI surface 사이에 연결된 cavity가 있는지, "
                "protein core가 가로막는지 PyMOL에서 확인"
            ),
        ))
    return items


def _rule4_cross_state_missing(
    valid_sites: List[dict],
) -> List[_Item]:
    """규칙 4: Cross-state 매칭 없는 STRONG 포켓 (SHOULD)."""
    items = []
    for row in valid_sites:
        if row.get("verdict", "").strip() != "STRONG":
            continue
        if _safe_float(row.get("cross_receptor_score", "1")) != 0:
            continue
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        score = row.get("confidence_score", "")
        vina_q = row.get("vina_quality_score", "")

        is_3gt8 = rid == "3GT8_raw"
        extra_ctx = ""
        if is_3gt8:
            extra_ctx = "\n         3GT8의 V924R(PDB 948) crystallization mutation 영향 가능"
        action = (
            "PyMOL에서 residue 948과 이 포켓의 거리 측정. 10A 이내면 artifact 가능성"
            if is_3gt8
            else "다른 state에서 이 영역의 conformation 차이 확인"
        )
        items.append(_Item(
            id=f"H-5.3-{rid}_{pid}",
            priority="SHOULD",
            title=f"{rid} 단독 STRONG 포켓 — cross-state 매칭 없음",
            context=(
                f"Score {score}, Vina {vina_q}. 다른 receptor state에서 매칭 없음."
                + extra_ctx
            ),
            action=action,
        ))
    return items


def _rule5_bootstrap_low(
    valid_sites: List[dict],
) -> List[_Item]:
    """규칙 5: Bootstrap stability 낮은 포켓 (NOTE)."""
    items = []
    for row in valid_sites:
        stability = _safe_float(row.get("pocket_stability", "1.0"), 1.0)
        if stability >= 0.6:
            continue
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        n_pose = row.get("n_pose", "")
        items.append(_Item(
            id=f"H-4.2-{rid}_{pid}",
            priority="NOTE",
            title=f"Bootstrap stability {stability:.2f} (낮음)",
            context=f"0.6 미만은 리샘플링에서 포켓이 불안정. n_pose={n_pose}",
            action=(
                "이 포켓이 소수 pose에 의존하는지 확인. "
                "stability가 낮아도 cross-state에서 재현되면 유의미할 수 있음."
            ),
        ))
    return items


def _rule6_validation_warnings(
    validation_text: str,
) -> List[_Item]:
    """규칙 6: Validation warning 기반 (SHOULD/NOTE)."""
    if not validation_text:
        return []

    # "Warning summary:" 이후 텍스트에서 warning 추출
    idx = validation_text.find("Warning summary:")
    if idx < 0:
        # 대체: "WARNING" 포함 행 수집
        idx = validation_text.find("WARNING")
        if idx < 0:
            return []

    warning_block = validation_text[idx:]
    lines = [ln.strip() for ln in warning_block.splitlines() if ln.strip()]

    items = []
    seen_types = set()
    for line in lines:
        lower = line.lower()
        if "known mutation" in lower and "mutation" not in seen_types:
            seen_types.add("mutation")
            items.append(_Item(
                id="H-6.1-mutation",
                priority="SHOULD",
                title="Known mutation warning 검토",
                context=line,
                action=(
                    "해당 mutation이 STRONG 포켓의 핵심 잔기에 해당하는지 확인"
                ),
            ))
        elif "chain" in lower and "chain" not in seen_types:
            seen_types.add("chain")
            items.append(_Item(
                id="H-6.1-chain",
                priority="NOTE",
                title="Chain ID 차이 확인됨",
                context=line,
                action="Chain ID 차이가 결과 해석에 영향을 주는지 확인",
            ))
        elif "optional" in lower and "optional" not in seen_types:
            seen_types.add("optional")
            items.append(_Item(
                id="H-6.2-optional",
                priority="NOTE",
                title="Optional 파일 부재",
                context=line,
                action="누락이 의도적인지 (예: AFM 미실행) 확인",
            ))
    return items


def _rule7_ppi_low_reproducibility(
    ppi_residues: List[dict],
    agreement: List[dict],
    valid_sites: List[dict],
) -> List[_Item]:
    """규칙 7: PPI seed 재현성 낮은 잔기가 STRONG 포켓의 shared residue에 포함 (NOTE)."""
    # STRONG 포켓 set
    strong_set = set()
    for row in valid_sites:
        if row.get("verdict", "").strip() == "STRONG":
            strong_set.add((
                row.get("receptor_id", "").strip(),
                row.get("pocket_id", "").strip(),
            ))

    # agreement에서 STRONG 포켓의 shared_residue_list 수집
    strong_shared_nums: set = set()
    for row in agreement:
        rid = row.get("receptor_id", "").strip()
        pid = row.get("pocket_id", "").strip()
        if (rid, pid) not in strong_set:
            continue
        shared_raw = row.get("shared_residue_list", "").strip()
        if shared_raw:
            for token in re.split(r"[;,\s]+", shared_raw):
                m = re.search(r"(\d+)", token)
                if m:
                    strong_shared_nums.add(int(m.group(1)))

    if not strong_shared_nums:
        return []

    # PPI residues에서 C-lobe + low frac + occupancy >= 0.05
    low_robust = []
    for row in ppi_residues:
        if row.get("lobe_label", "").strip() != "C-lobe":
            continue
        frac = _safe_float(row.get("frac_runs_supporting", "1.0"), 1.0)
        occ = _safe_float(row.get("occupancy", "0"), 0.0)
        if frac >= 0.6 or occ < 0.05:
            continue
        resnum = _safe_int(row.get("residue_num", "0"))
        if resnum in strong_shared_nums:
            low_robust.append(f"{row.get('residue_name', '')} {resnum} (frac={frac:.2f})")

    if not low_robust:
        return []

    return [_Item(
        id="H-3.2-low_robustness",
        priority="NOTE",
        title=f"{len(low_robust)}개 PPI 잔기의 seed 재현성 낮음 (frac_runs < 0.6)",
        context=f"{', '.join(low_robust)}이 STRONG 포켓의 shared residue에 포함됨",
        action="이 잔기가 verdict 점수에 기여한 정도를 ppi_overlap_pts에서 확인",
    )]


# ---------------------------------------------------------------------------
# Format + write
# ---------------------------------------------------------------------------

def _format_item(item: _Item) -> str:
    """단일 checklist item을 마크다운 블록으로 포맷."""
    lines = [f"### [ ] {item.id}: {item.title}"]
    lines.append(f"- **Context**: {item.context}")
    lines.append(f"- **Action**: {item.action}")
    if item.pymol_hint:
        lines.append(f"- **PyMOL**: `{item.pymol_hint}`")
    if item.priority in ("MUST", "SHOULD"):
        lines.append("- **판단 결과**: _(여기에 기록)_")
    lines.append("")
    return "\n".join(lines)


def generate_human_review_checklist(project_root: Path) -> Path:
    """Production 완료 후 context-aware human review checklist를 생성한다.

    읽는 파일:
      - valid_sites.csv
      - cross_method_agreement.csv
      - vina_pocket_table.csv
      - ppi_pyrosetta_residues.csv
      - validation_summary.txt (또는 step7_validate/validation_summary.txt)
      - vina_consensus_sites.csv
      - combined_residue_evidence.csv

    생성하는 파일:
      - {project_root}/human_review_checklist.md

    Returns:
      생성된 체크리스트 파일 경로
    """
    from egfr_pipeline.paths import (
        wa_phase3_ppi_postprocess,
        wa_phase4_vina_postprocess,
        wa_phase5_verdict,
        wa_phase7_validation,
    )

    verdict_dir = wa_phase5_verdict()
    vina_post_dir = wa_phase4_vina_postprocess()
    ppi_post_dir = wa_phase3_ppi_postprocess()
    validation_dir = wa_phase7_validation()

    # Load data
    valid_sites = _read_csv_rows(verdict_dir / "valid_sites.csv")
    agreement = _read_csv_rows(verdict_dir / "cross_method_agreement.csv")
    pocket_table = _read_csv_rows(vina_post_dir / "vina_pocket_table.csv")
    ppi_residues = _read_csv_rows(ppi_post_dir / "ppi_pyrosetta_residues.csv")

    # validation_summary.txt: try multiple locations
    validation_text = _read_text(validation_dir / "validation_summary.txt")
    if not validation_text:
        validation_text = _read_text(project_root / "validation_summary.txt")

    # Collect items via all 7 rules
    all_items: List[_Item] = []

    if valid_sites:
        all_items.extend(_rule1_strong_visual(valid_sites, pocket_table))
        all_items.extend(_rule2_atp_overlap(valid_sites, pocket_table))
        all_items.extend(_rule3_ppi_moderate(valid_sites, agreement))
        all_items.extend(_rule4_cross_state_missing(valid_sites))
        all_items.extend(_rule5_bootstrap_low(valid_sites))
    else:
        all_items.append(_Item(
            id="H-0.1-no_data",
            priority="NOTE",
            title="valid_sites.csv 데이터 없음",
            context="verdict 결과 파일이 비어있거나 존재하지 않음",
            action="Phase 5 재실행 필요 여부 확인",
        ))

    all_items.extend(_rule6_validation_warnings(validation_text))

    if ppi_residues and agreement and valid_sites:
        all_items.extend(_rule7_ppi_low_reproducibility(
            ppi_residues, agreement, valid_sites,
        ))

    # Classify by priority
    must_items = [i for i in all_items if i.priority == "MUST"]
    should_items = [i for i in all_items if i.priority == "SHOULD"]
    note_items = [i for i in all_items if i.priority == "NOTE"]

    # Count verdicts
    n_strong = sum(1 for r in valid_sites if r.get("verdict", "").strip() == "STRONG")
    n_moderate = sum(1 for r in valid_sites if r.get("verdict", "").strip() == "MODERATE")
    n_weak = sum(1 for r in valid_sites if r.get("verdict", "").strip() == "WEAK")

    project_name = project_root.resolve().name
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Build markdown
    lines = [
        "# Human Review Checklist",
        "",
        f"Generated: {timestamp} after production run",
        f"Project: {project_name}",
        f"Based on: valid_sites.csv ({n_strong} STRONG, {n_moderate} MODERATE, {n_weak} WEAK)",
        "",
        "---",
        "",
    ]

    # MUST
    lines.append("## MUST — 반드시 확인 (진행 전 필수)")
    lines.append("")
    if must_items:
        for item in must_items:
            lines.append(_format_item(item))
    else:
        lines.append("_(해당 항목 없음)_")
        lines.append("")

    lines.append("---")
    lines.append("")

    # SHOULD
    lines.append("## SHOULD — 권장 확인")
    lines.append("")
    if should_items:
        for item in should_items:
            lines.append(_format_item(item))
    else:
        lines.append("_(해당 항목 없음)_")
        lines.append("")

    lines.append("---")
    lines.append("")

    # NOTE
    lines.append("## NOTE — 참고 사항")
    lines.append("")
    if note_items:
        for item in note_items:
            lines.append(_format_item(item))
    else:
        lines.append("_(해당 항목 없음)_")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("체크 완료 후 이 파일에 판단 결과를 기록하고 보관하세요.")
    lines.append("각 항목: ✓ (확인, 문제없음) / ✗ (문제 발견, 조치 필요) / — (해당 없음)")
    lines.append("")

    # Write
    output_path = project_root / "human_review_checklist.md"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")

    return output_path


def checklist_summary(project_root: Path) -> Tuple[int, int, int]:
    """체크리스트의 MUST/SHOULD/NOTE 항목 수를 반환. 파일 없으면 (0, 0, 0)."""
    checklist = project_root / "human_review_checklist.md"
    if not checklist.exists():
        return 0, 0, 0
    content = checklist.read_text(encoding="utf-8")
    n_must = content.count("Priority: MUST") if "Priority: MUST" in content else 0
    n_should = content.count("Priority: SHOULD") if "Priority: SHOULD" in content else 0
    n_note = content.count("Priority: NOTE") if "Priority: NOTE" in content else 0
    # Fallback: section headers 기반
    if n_must == 0 and n_should == 0 and n_note == 0:
        # H-x.x ID 패턴 기반 카운트
        must_section = content.split("## SHOULD")[0] if "## SHOULD" in content else ""
        should_section = content.split("## SHOULD")[1].split("## NOTE")[0] if "## SHOULD" in content and "## NOTE" in content else ""
        note_section = content.split("## NOTE")[1] if "## NOTE" in content else ""
        n_must = must_section.count("### [ ]")
        n_should = should_section.count("### [ ]")
        n_note = note_section.count("### [ ]")
    return n_must, n_should, n_note
