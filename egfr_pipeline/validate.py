"""Output validation, regression checks, and handoff readiness.

Runs after a pipeline execution to verify:
1. Core output files exist
2. Receptor/ligand IDs are consistent across all output tables
3. CSV field structures match expected schemas (regression guard)
4. Residue numbering is consistent across receptor PDBs
5. Repository contains required handoff documents

Exit code 0 = all checks pass, 1 = warnings only, 2 = failures found.
"""
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.config import load_config, project_root_from_config


# ---------------------------------------------------------------------------
# Result tracking
# ---------------------------------------------------------------------------

class ValidationResult:
    def __init__(self):
        self.passes: List[str] = []
        self.warnings: List[str] = []
        self.failures: List[str] = []

    def ok(self, msg: str):
        self.passes.append(msg)

    def warn(self, msg: str):
        self.warnings.append(msg)

    def fail(self, msg: str):
        self.failures.append(msg)

    @property
    def exit_code(self) -> int:
        if self.failures:
            return 2
        if self.warnings:
            return 1
        return 0

    def summary(self) -> str:
        lines = []
        lines.append(f"\n{'='*60}")
        lines.append("VALIDATION SUMMARY")
        lines.append(f"{'='*60}")
        lines.append(f"  PASS:    {len(self.passes)}")
        lines.append(f"  WARN:    {len(self.warnings)}")
        lines.append(f"  FAIL:    {len(self.failures)}")
        lines.append("")
        if self.failures:
            lines.append("FAILURES:")
            for f in self.failures:
                lines.append(f"  [FAIL] {f}")
            lines.append("")
        if self.warnings:
            lines.append("WARNINGS:")
            for w in self.warnings:
                lines.append(f"  [WARN] {w}")
            lines.append("")
        if not self.failures and not self.warnings:
            lines.append("All checks passed.")
        lines.append(f"Exit code: {self.exit_code}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_csv_safe(path: Path) -> Tuple[List[dict], List[str]]:
    """Load CSV, return (rows, fieldnames). Empty if file missing."""
    if not path.exists():
        return [], []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return rows, list(reader.fieldnames) if reader.fieldnames else []


# ---------------------------------------------------------------------------
# 8.1 Output existence and consistency
# ---------------------------------------------------------------------------

CORE_OUTPUTS = [
    "vina_pose_table.csv",
    "vina_pocket_table.csv",
    "vina_drug_pocket_map.csv",
]

OPTIONAL_OUTPUTS = [
    "vina_pocket_comparison.csv",
    "ppi_pyrosetta_residues.csv",
    "ppi_pyrosetta_summary.csv",
    "ppi_afm_residues.csv",
    "project_report.txt",
    "combined_residue_evidence.csv",
    "cross_method_agreement.csv",
    "valid_sites.csv",
]


def check_output_existence(project_root: Path, result: ValidationResult):
    """8.1: Verify core outputs exist."""
    for name in CORE_OUTPUTS:
        path = project_root / name
        if path.exists():
            rows, _ = load_csv_safe(path)
            if rows:
                result.ok(f"{name} exists ({len(rows)} rows)")
            else:
                result.warn(f"{name} exists but is empty")
        else:
            result.fail(f"{name} missing")

    for name in OPTIONAL_OUTPUTS:
        path = project_root / name
        if path.exists():
            result.ok(f"{name} exists (optional)")
        else:
            result.warn(f"{name} not found (optional)")


def check_id_consistency(project_root: Path, config: dict, result: ValidationResult):
    """8.1: Verify receptor/ligand IDs are consistent across outputs."""
    expected_receptors = {r["id"] for r in config.get("receptors", [])}
    expected_ligands = {l["id"] for l in config.get("ligands", [])}

    if not expected_receptors:
        result.warn("No receptors defined in config")
        return

    files_to_check = {
        "vina_pose_table.csv": ("receptor_id", "ligand_id"),
        "vina_pocket_table.csv": ("receptor_id", None),
        "vina_drug_pocket_map.csv": ("receptor_id", "ligand_id"),
    }

    for filename, (rec_col, lig_col) in files_to_check.items():
        rows, _ = load_csv_safe(project_root / filename)
        if not rows:
            continue

        found_receptors = {r[rec_col] for r in rows if rec_col in r}
        unexpected_recs = found_receptors - expected_receptors
        missing_recs = expected_receptors - found_receptors

        if unexpected_recs:
            result.fail(f"{filename}: unexpected receptor IDs: {unexpected_recs}")
        elif missing_recs:
            result.warn(f"{filename}: missing receptor IDs: {missing_recs}")
        else:
            result.ok(f"{filename}: receptor IDs consistent")

        if lig_col:
            found_ligands = {r[lig_col] for r in rows if lig_col in r}
            unexpected_ligs = found_ligands - expected_ligands
            if unexpected_ligs:
                result.fail(f"{filename}: unexpected ligand IDs: {unexpected_ligs}")
            else:
                result.ok(f"{filename}: ligand IDs consistent")


# ---------------------------------------------------------------------------
# 8.2 Regression checks (CSV schema stability)
# ---------------------------------------------------------------------------

EXPECTED_SCHEMAS = {
    "vina_pose_table.csv": [
        "receptor_id", "ligand_id", "pose_rank", "affinity",
        "rmsd_lb", "rmsd_ub", "centroid_x", "centroid_y", "centroid_z",
        "raw_pose_file", "pocket_id", "contact_residues", "n_contact_residues",
    ],
    "vina_pocket_table.csv": [
        "receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z",
        "n_pose", "n_ligand", "best_affinity", "mean_affinity",
        "union_contact_residues", "top_residues",
    ],
    "vina_drug_pocket_map.csv": [
        "receptor_id", "ligand_id", "dominant_pocket_id",
        "dominant_pocket_pose_count", "dominant_pocket_fraction",
        "best_affinity", "best_pose_rank", "top_pose_residues",
        "alternative_pockets", "is_multimodal_binding",
    ],
    "vina_pocket_comparison.csv": [
        "receptor_a", "pocket_a", "receptor_b", "pocket_b",
        "centroid_dist", "residue_jaccard", "residue_overlap_coeff",
        "shared_residues", "n_shared_residues",
        "residues_only_a", "residues_only_b",
        "n_residues_a", "n_residues_b",
        "shared_ligands", "n_shared_ligands", "n_ligands_a", "n_ligands_b",
        "affinity_a", "affinity_b", "n_pose_a", "n_pose_b",
        "same_patch_candidate",
    ],
    "ppi_pyrosetta_residues.csv": [
        "receptor_id", "source", "residue_id", "residue_num",
        "frequency_final_ranking", "frequency_cluster_summary",
        "n_models_final_ranking", "occupancy",
        "mean_interface_delta_e", "best_interface_delta_e",
    ],
    "combined_residue_evidence.csv": [
        "receptor_id", "residue_id", "vina_pockets", "n_vina_pockets",
        "ppi_occupancy", "ppi_frequency", "ppi_delta_e", "evidence_sources",
    ],
    "cross_method_agreement.csv": [
        "receptor_id", "pocket_id", "n_vina_residues", "n_ppi_residues",
        "n_shared_residues", "jaccard", "overlap_coeff", "shared_residue_list",
        "ppi_mean_occupancy_of_shared", "spatial_dist_A", "spatial_proximity",
        "vina_best_affinity_kcal", "ppi_best_dg_REU", "agreement_level",
    ],
    "valid_sites.csv": [
        "receptor_id", "pocket_id", "verdict", "confidence_score",
        "vina_quality_score", "ppi_proximity_score", "cross_receptor_score",
        "ppi_data_available", "best_affinity", "n_pose", "n_ligand",
        "spatial_dist_to_ppi", "n_shared_with_ppi",
        "cross_receptor_matches", "reasons",
    ],
}


def check_csv_schemas(project_root: Path, result: ValidationResult):
    """8.2: Verify CSV field structures haven't changed."""
    for filename, expected_fields in EXPECTED_SCHEMAS.items():
        path = project_root / filename
        if not path.exists():
            continue
        _, actual_fields = load_csv_safe(path)
        if not actual_fields:
            result.warn(f"{filename}: could not read header")
            continue

        missing = set(expected_fields) - set(actual_fields)
        extra = set(actual_fields) - set(expected_fields)

        if missing:
            result.fail(f"{filename}: missing expected columns: {sorted(missing)}")
        elif extra:
            result.warn(f"{filename}: extra columns (non-breaking): {sorted(extra)}")
        else:
            result.ok(f"{filename}: schema matches expected ({len(expected_fields)} columns)")


# ---------------------------------------------------------------------------
# 8.3 Residue numbering consistency
# ---------------------------------------------------------------------------

def parse_pdb_residue_set(pdb_path: Path) -> Dict[str, Set[int]]:
    """Extract {chain: set of residue numbers} from PDB."""
    chain_residues: Dict[str, Set[int]] = defaultdict(set)
    if not pdb_path.exists():
        return {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip() or "?"
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            chain_residues[chain].add(resnum)
    return dict(chain_residues)


def check_residue_numbering(config: dict, result: ValidationResult):
    """8.3: Compare residue numbering across receptor PDBs."""
    receptors = config.get("receptors", [])
    pdb_data: Dict[str, Dict[str, Set[int]]] = {}

    for rec in receptors:
        pdb_path = Path(rec.get("pdb", ""))
        if not pdb_path.exists():
            result.warn(f"Receptor PDB not found for {rec['id']}: {pdb_path}")
            continue
        pdb_data[rec["id"]] = parse_pdb_residue_set(pdb_path)

    if len(pdb_data) < 2:
        result.warn("Need at least 2 receptor PDBs for numbering comparison")
        return

    # Compare all pairs
    receptor_ids = sorted(pdb_data.keys())
    for i, rec_a in enumerate(receptor_ids):
        for rec_b in receptor_ids[i+1:]:
            chains_a = pdb_data[rec_a]
            chains_b = pdb_data[rec_b]

            # Find the primary chain in each (largest residue count)
            main_chain_a = max(chains_a, key=lambda c: len(chains_a[c])) if chains_a else None
            main_chain_b = max(chains_b, key=lambda c: len(chains_b[c])) if chains_b else None

            if not main_chain_a or not main_chain_b:
                result.warn(f"{rec_a} vs {rec_b}: could not identify main chains")
                continue

            res_a = chains_a[main_chain_a]
            res_b = chains_b[main_chain_b]
            overlap = res_a & res_b
            only_a = res_a - res_b
            only_b = res_b - res_a

            if main_chain_a != main_chain_b:
                result.warn(
                    f"{rec_a}(chain {main_chain_a}) vs {rec_b}(chain {main_chain_b}): "
                    f"different chain IDs -- cross-receptor residue comparison requires chain-stripping"
                )

            if overlap:
                overlap_pct = len(overlap) / max(len(res_a), len(res_b)) * 100
                result.ok(
                    f"{rec_a} vs {rec_b}: {len(overlap)} shared residue numbers "
                    f"({overlap_pct:.0f}% of larger set), "
                    f"{len(only_a)} only in {rec_a}, {len(only_b)} only in {rec_b}"
                )
                if overlap_pct < 50:
                    result.warn(
                        f"{rec_a} vs {rec_b}: low overlap ({overlap_pct:.0f}%) -- "
                        f"residue-level comparison may be unreliable"
                    )
            else:
                result.fail(
                    f"{rec_a} vs {rec_b}: NO overlapping residue numbers -- "
                    f"cross-receptor comparison is unsafe"
                )


# ---------------------------------------------------------------------------
# 8.4 Handoff readiness
# ---------------------------------------------------------------------------

HANDOFF_DOCS = [
    ("README.md", True),
    ("CLAUDE.md", True),
    ("CODEX_CONTINUATION_2026-03-09.md", True),
    ("codex_handoff_egfr_myo_1_d_pipeline_spec.md", True),
    ("tasks_egfr_myo_1_d_pipeline.md", True),
    ("prd_egfr_myo_1_d_pipeline.md", False),
    ("brief_egfr_myo_1_d_pipeline.md", False),
    ("MANUAL.md", False),
]

HANDOFF_MODULES = [
    "run_docking.py",
    "parse_vina_results.py",
    "extract_contacts.py",
    "cluster_pockets.py",
    "summarize_pockets.py",
    "compare_pockets.py",
    "extract_ppi_residues.py",
    "generate_report.py",
    "validate_outputs.py",
]


def check_handoff_readiness(repo_root: Path, result: ValidationResult):
    """8.4: Verify repository has required docs and modules."""
    for name, required in HANDOFF_DOCS:
        path = repo_root / name
        if path.exists():
            size = path.stat().st_size
            if size > 100:
                result.ok(f"Doc {name} exists ({size} bytes)")
            else:
                result.warn(f"Doc {name} exists but very small ({size} bytes)")
        elif required:
            result.fail(f"Required doc missing: {name}")
        else:
            result.warn(f"Optional doc missing: {name}")

    for name in HANDOFF_MODULES:
        path = repo_root / name
        if path.exists():
            result.ok(f"Module {name} exists")
        else:
            result.fail(f"Module missing: {name}")

    # Check config examples
    example_config = repo_root / "config" / "example-project.yaml"
    if example_config.exists():
        result.ok("Example project config exists")
    else:
        result.warn("config/example-project.yaml missing")


# ---------------------------------------------------------------------------
# Traceability check
# ---------------------------------------------------------------------------

def check_traceability(project_root: Path, result: ValidationResult):
    """8.1: Verify parsed outputs remain connected to raw sources."""
    pose_rows, _ = load_csv_safe(project_root / "vina_pose_table.csv")
    if not pose_rows:
        return

    missing_raw = 0
    checked = 0
    for row in pose_rows:
        raw_file = row.get("raw_pose_file", "")
        if raw_file:
            checked += 1
            if not Path(raw_file).exists():
                missing_raw += 1

    if checked == 0:
        result.warn("No raw_pose_file references in pose table")
    elif missing_raw == 0:
        result.ok(f"All {checked} raw pose files are accessible")
    elif missing_raw == checked:
        result.warn(f"All {checked} raw pose files are inaccessible (expected if moved from original run location)")
    else:
        result.warn(f"{missing_raw}/{checked} raw pose files are inaccessible")


# ---------------------------------------------------------------------------
# Main validation runner
# ---------------------------------------------------------------------------

def run_validation(
    config_path: str,
    repo_root: Optional[str] = None,
) -> ValidationResult:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    repo = Path(repo_root) if repo_root else Path(".")

    result = ValidationResult()

    print("Running validation checks...")
    print(f"  Config: {config_path}")
    print(f"  Project root: {project_root}")
    print(f"  Repo root: {repo}")
    print()

    # 8.1
    print("[8.1] Output existence and consistency...")
    check_output_existence(project_root, result)
    check_id_consistency(project_root, config, result)
    check_traceability(project_root, result)

    # 8.2
    print("[8.2] CSV schema regression checks...")
    check_csv_schemas(project_root, result)

    # 8.3
    print("[8.3] Residue numbering consistency...")
    check_residue_numbering(config, result)

    # 8.4
    print("[8.4] Handoff readiness...")
    check_handoff_readiness(repo, result)

    return result
