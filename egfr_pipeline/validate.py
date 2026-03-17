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
    "vina_postprocess_coverage.csv",
    "ppi_pyrosetta_residues.csv",
    "ppi_pyrosetta_summary.csv",
    "ppi_afm_residues.csv",
    "project_report.txt",
    "combined_residue_evidence.csv",
    "cross_method_agreement.csv",
    "valid_sites.csv",
    "vina_consensus_sites.csv",
    "vina_pocket_bootstrap.csv",
    "vina_contact_distances.csv",
    "vina_pocket_residue_occupancy.csv",
    "vina_clustering_merge_log.csv",
    "vina_clustering_parameters.json",
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
        "vina_postprocess_coverage.csv": ("receptor_id", "ligand_id"),
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
        "raw_pose_file", "docking_mode", "exhaustiveness", "requested_n_poses",
        "energy_range", "cpu_per_job", "vina_seed", "pose_source_status",
        "pocket_id", "contact_residues", "n_contact_residues", "contact_distances",
    ],
    "vina_pocket_table.csv": [
        "receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z",
        "n_pose", "n_ligand", "best_affinity", "mean_affinity",
        "union_contact_residues", "top_residues",
        "centroid_spread_A", "affinity_std", "affinity_iqr",
        "dominant_ligand_fraction", "ligand_pose_entropy",
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
        "centroid_dist_bootstrap_ci",
    ],
    "ppi_pyrosetta_residues.csv": [
        "receptor_id", "partner_id", "source", "chain", "residue_id", "residue_num",
        "residue_name", "lobe_label", "construct_type", "orientation_validation_status",
        "n_runs_total", "n_runs_supporting", "frac_runs_supporting", "supporting_seed_indices",
        "frequency_final_ranking", "frequency_cluster_summary",
        "n_models_final_ranking", "occupancy",
        "mean_interface_delta_e", "best_interface_delta_e",
    ],
    "ppi_pyrosetta_summary.csv": [
        "receptor_id", "partner_id", "source", "construct_type", "orientation_validation_status",
        "n_runs_total", "n_runs_completed", "seed_indices",
        "n_final_models", "n_clusters", "n_interface_residues",
        "n_nlobe_interface_residues", "n_clobe_interface_residues",
        "top_residues", "best_dg", "mean_dg", "best_dsasa",
    ],
    "ppi_pyrosetta_residue_long.csv": [
        "model_id", "receptor_id", "partner_id", "source", "construct_type",
        "orientation_validation_status", "seed_index", "run_label", "run_dir",
        "rank", "cluster_id", "chain", "residue_id", "residue_num",
        "residue_name", "lobe_label", "delta_e_total", "delta_e_fa_atr",
        "delta_e_fa_rep", "delta_e_fa_sol", "delta_e_fa_elec", "source_file",
    ],
    "ppi_pyrosetta_model_table.csv": [
        "model_id", "receptor_id", "partner_id", "source", "construct_type",
        "orientation_validation_status", "seed_index", "run_label", "run_dir",
        "rank", "cluster_id", "dG_separated", "dSASA", "sc_value", "packstat",
        "nres_int", "n_receptor_interface_residues", "n_partner_interface_residues",
        "n_nlobe_interface_residues", "n_clobe_interface_residues",
        "receptor_interface_residues", "partner_interface_residues", "source_file",
    ],
    "ppi_afm_residues.csv": [
        "receptor_id", "source", "residue_id", "residue_num",
        "min_ca_distance",
    ],
    "combined_residue_evidence.csv": [
        "receptor_id", "residue_id", "vina_pockets", "n_vina_pockets",
        "vina_best_affinity", "vina_best_pocket_stability",
        "vina_best_confidence_score", "vina_cross_receptor_support",
        "vina_evidence_profile",
        "ppi_occupancy", "ppi_frequency", "ppi_delta_e",
        "ppi_frac_runs_supporting", "ppi_partners", "evidence_sources",
    ],
    "cross_method_agreement.csv": [
        "receptor_id", "pocket_id", "n_vina_residues", "n_ppi_residues",
        "n_shared_residues", "jaccard", "overlap_coeff", "shared_residue_list",
        "ppi_mean_occupancy_of_shared", "spatial_dist_A", "spatial_proximity",
        "closest_ppi_partner", "n_ppi_partners_near",
        "ppi_frac_runs_supporting", "ppi_best_frac_runs_supporting",
        "ppi_best_interface_delta_e", "ppi_shared_partner_count",
        "vina_best_affinity_kcal", "ppi_best_dg_REU", "agreement_level",
    ],
    "valid_sites.csv": [
        "receptor_id", "pocket_id", "verdict", "confidence_score",
        "vina_quality_score", "ppi_proximity_score", "cross_receptor_score",
        "vina_affinity_pts", "vina_convergence_pts", "vina_stability_pts",
        "vina_diversity_pts", "vina_consensus_pts",
        "ppi_spatial_pts", "ppi_overlap_pts", "ppi_reproducibility_pts",
        "cross_receptor_pts", "cross_receptor_support_pts",
        "score_denominator",
        "ppi_data_available", "best_affinity", "n_pose", "n_ligand",
        "dominant_ligand_fraction", "ligand_pose_entropy",
        "spatial_dist_to_ppi", "closest_ppi_partner", "n_ppi_partners_near",
        "n_shared_with_ppi", "ppi_frac_runs_supporting", "ppi_best_interface_delta_e",
        "cross_receptor_support",
        "cross_receptor_matches", "consensus_site_id",
        "exp_sensitivity", "exp_specificity", "exp_enrichment", "exp_rank_impact",
        "pocket_stability", "evidence_profile", "reason_tags", "decision_trace",
        "reasons",
    ],
    "vina_consensus_sites.csv": [
        "consensus_site_id", "n_receptors", "receptor_list", "pocket_list",
        "centroid_x", "centroid_y", "centroid_z",
        "best_affinity", "total_n_ligand", "total_n_pose", "consensus_residues",
    ],
    "vina_pocket_bootstrap.csv": [
        "receptor_id", "pocket_id", "pocket_exists_frac", "centroid_std_A",
        "affinity_mean", "affinity_std", "affinity_iqr",
        "n_pose_mean", "n_pose_std", "n_replicates", "sample_fraction",
        "stability_scope",
    ],
    "vina_postprocess_coverage.csv": [
        "receptor_id", "ligand_id", "raw_pose_file", "status",
        "parsed_n_poses", "requested_n_poses", "coverage_fraction",
        "docking_mode", "exhaustiveness", "energy_range",
        "cpu_per_job", "base_seed", "vina_seed",
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


def parse_pdb_residue_identity(pdb_path: Path) -> Dict[str, Dict[int, str]]:
    """Extract {chain: {resnum: resname}} from PDB (CA atoms only for speed).

    Returns the 3-letter residue name for each residue number, enabling
    identity verification across PDBs (same number = same amino acid?).
    """
    chain_res: Dict[str, Dict[int, str]] = defaultdict(dict)
    if not pdb_path.exists():
        return {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21].strip() or "?"
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            chain_res[chain][resnum] = resname
    return dict(chain_res)


def _detect_offset(
    seq_a: Dict[int, str], seq_b: Dict[int, str],
) -> Optional[int]:
    """Detect a systematic residue number offset between two PDB sequences.

    Tries offsets from -500 to +500. Returns the offset that maximizes
    the number of matching (resnum, resname) pairs, or None if no
    significant offset is found.
    """
    if not seq_a or not seq_b:
        return None

    # First check offset=0 (ideal case)
    direct_match = sum(1 for r in seq_a if r in seq_b and seq_a[r] == seq_b[r])
    min_len = min(len(seq_a), len(seq_b))
    if min_len == 0:
        return None
    if direct_match / min_len > 0.8:
        return 0

    # Try candidate offsets based on start-residue difference
    sorted_a = sorted(seq_a.keys())
    sorted_b = sorted(seq_b.keys())
    candidate_offsets = set()
    # Try difference of first/last residues
    for sa in sorted_a[:3] + sorted_a[-3:]:
        for sb in sorted_b[:3] + sorted_b[-3:]:
            candidate_offsets.add(sb - sa)
    # Also try common small offsets
    for o in range(-500, 501, 1):
        candidate_offsets.add(o)

    best_offset = 0
    best_matches = direct_match

    for offset in candidate_offsets:
        if offset == 0:
            continue
        matches = sum(
            1 for r in seq_a
            if (r + offset) in seq_b and seq_a[r] == seq_b[r + offset]
        )
        if matches > best_matches:
            best_matches = matches
            best_offset = offset

    if best_matches / min_len > 0.5 and best_offset != 0:
        return best_offset
    return 0 if direct_match / min_len > 0.3 else None


def check_residue_numbering(config: dict, result: ValidationResult):
    """8.3: Compare residue numbering and identity across receptor PDBs.

    Three-level check:
      1. Number overlap (do the same residue numbers appear?)
      2. Identity match (same number = same amino acid?)
      3. Offset detection (systematic shift in numbering?)
    """
    receptors = config.get("receptors", [])
    pdb_nums: Dict[str, Dict[str, Set[int]]] = {}
    pdb_ids: Dict[str, Dict[str, Dict[int, str]]] = {}

    for rec in receptors:
        pdb_path = Path(rec.get("pdb", ""))
        if not pdb_path.exists():
            result.warn(f"Receptor PDB not found for {rec['id']}: {pdb_path}")
            continue
        pdb_nums[rec["id"]] = parse_pdb_residue_set(pdb_path)
        pdb_ids[rec["id"]] = parse_pdb_residue_identity(pdb_path)

    if len(pdb_nums) < 2:
        result.warn("Need at least 2 receptor PDBs for numbering comparison")
        return

    receptor_ids = sorted(pdb_nums.keys())
    for i, rec_a in enumerate(receptor_ids):
        for rec_b in receptor_ids[i+1:]:
            chains_a = pdb_nums[rec_a]
            chains_b = pdb_nums[rec_b]

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
                    f"different chain IDs — chain-stripping applied for comparison"
                )

            # Level 1: Number overlap
            if overlap:
                overlap_pct = len(overlap) / max(len(res_a), len(res_b)) * 100
                result.ok(
                    f"{rec_a} vs {rec_b}: {len(overlap)} shared residue numbers "
                    f"({overlap_pct:.0f}% of larger set), "
                    f"{len(only_a)} only in {rec_a}, {len(only_b)} only in {rec_b}"
                )
                if overlap_pct < 50:
                    result.warn(
                        f"{rec_a} vs {rec_b}: low overlap ({overlap_pct:.0f}%) — "
                        f"residue-level comparison may be unreliable"
                    )
            else:
                result.fail(
                    f"{rec_a} vs {rec_b}: NO overlapping residue numbers — "
                    f"cross-receptor comparison is unsafe"
                )

            # Level 2: Identity match (same resnum = same resname?)
            id_a = pdb_ids.get(rec_a, {}).get(main_chain_a, {})
            id_b = pdb_ids.get(rec_b, {}).get(main_chain_b, {})

            if id_a and id_b and overlap:
                mismatches = []
                matched = 0
                for resnum in sorted(overlap):
                    name_a = id_a.get(resnum)
                    name_b = id_b.get(resnum)
                    if name_a and name_b:
                        if name_a == name_b:
                            matched += 1
                        else:
                            # Allow CHARMM name variants (HSD/HSE/HSP → HIS)
                            norm_a = "HIS" if name_a in ("HSD", "HSE", "HSP") else name_a
                            norm_b = "HIS" if name_b in ("HSD", "HSE", "HSP") else name_b
                            if norm_a == norm_b:
                                matched += 1
                            else:
                                mismatches.append((resnum, name_a, name_b))

                if mismatches:
                    # Separate known mutations from unexpected mismatches
                    known = [(r, a, b) for r, a, b in mismatches if r in KNOWN_MUTATIONS]
                    unknown = [(r, a, b) for r, a, b in mismatches if r not in KNOWN_MUTATIONS]

                    if known:
                        kn_str = ", ".join(f"{r}:{a}≠{b}(known mutation)" for r, a, b in known)
                        result.warn(
                            f"{rec_a} vs {rec_b}: {len(known)} known mutation difference(s): {kn_str}"
                        )
                    if unknown:
                        n_mm = len(unknown)
                        examples = unknown[:5]
                        ex_str = ", ".join(f"{r}:{a}≠{b}" for r, a, b in examples)
                        if n_mm > 5:
                            ex_str += f" ...+{n_mm - 5} more"
                        result.fail(
                            f"{rec_a} vs {rec_b}: {n_mm} UNEXPECTED identity mismatches! "
                            f"Same number, different amino acid: {ex_str}"
                        )
                elif matched > 0:
                    result.ok(
                        f"{rec_a} vs {rec_b}: {matched}/{len(overlap)} "
                        f"residue identities verified (same number = same amino acid)"
                    )

            # Level 3: Offset detection
            if id_a and id_b:
                offset = _detect_offset(id_a, id_b)
                if offset is None:
                    result.warn(
                        f"{rec_a} vs {rec_b}: could not determine numbering relationship"
                    )
                elif offset != 0:
                    # Verify the offset
                    offset_matches = sum(
                        1 for r in id_a
                        if (r + offset) in id_b and id_a[r] == id_b.get(r + offset, "")
                    )
                    result.fail(
                        f"{rec_a} vs {rec_b}: NUMBERING OFFSET DETECTED! "
                        f"offset={offset:+d} ({offset_matches} residues match with shift). "
                        f"e.g., {rec_a}:res{sorted(id_a.keys())[0]} = "
                        f"{rec_b}:res{sorted(id_a.keys())[0] + offset}. "
                        f"Cross-receptor residue comparison is INVALID without correction!"
                    )


# ---------------------------------------------------------------------------
# 8.3b Known mutations check
# ---------------------------------------------------------------------------

# Known crystallization/engineering mutations that differ from wild-type.
# Format: {resnum: (wt_resname, mutant_resname, receptor_id, description)}
KNOWN_MUTATIONS = {
    948: ("VAL", "ARG", "3GT8_raw",
          "V924R (UniProt V924, PDB 948) — Jura et al. 2009 Cell. "
          "Asymmetric dimer interface mutation for crystallization. "
          "αI helix region; may affect PPI interface energies near C02 site."),
}


def check_known_mutations(config: dict, result: ValidationResult):
    """8.3b: Flag known engineering mutations in receptor PDBs.

    Checks if crystallization artifacts or engineering mutations are present
    in receptor structures used for docking. These may distort binding
    energies, especially for PPI docking near the mutation site.
    """
    for rec in config.get("receptors", []):
        pdb_path = Path(rec.get("pdb", ""))
        rec_id = rec["id"]
        if not pdb_path.exists():
            continue

        identity = parse_pdb_residue_identity(pdb_path)
        # Get primary chain
        if not identity:
            continue
        main_chain = max(identity, key=lambda c: len(identity[c]))
        seq = identity[main_chain]

        for resnum, (wt, mut, affected_rec, desc) in KNOWN_MUTATIONS.items():
            actual = seq.get(resnum)
            if actual is None:
                continue
            # Normalize CHARMM names
            actual_norm = "HIS" if actual in ("HSD", "HSE", "HSP") else actual

            if actual_norm == mut and rec_id == affected_rec:
                result.warn(
                    f"{rec_id}: KNOWN MUTATION at residue {resnum}: "
                    f"{actual}(mutant) present, wild-type={wt}. {desc}"
                )
            elif actual_norm == wt:
                result.ok(
                    f"{rec_id}: residue {resnum} is wild-type {wt} (correct)"
                )


# ---------------------------------------------------------------------------
# 8.4 Handoff readiness
# ---------------------------------------------------------------------------

HANDOFF_DOCS = [
    ("README.md", True),
    ("CLAUDE.md", True),
    ("docs/manual_execution.md", False),
]

HANDOFF_MODULES = [
    "egfr_pipeline/__init__.py",
    "egfr_pipeline/config.py",
    "egfr_pipeline/vina/vina_executor.py",
    "egfr_pipeline/vina/parse_poses.py",
    "egfr_pipeline/vina/pose_contacts.py",
    "egfr_pipeline/vina/pocket_cluster.py",
    "egfr_pipeline/vina/pocket_summary.py",
    "egfr_pipeline/vina/cross_receptor.py",
    "egfr_pipeline/vina/pocket_stability.py",
    "egfr_pipeline/verdict.py",
    "egfr_pipeline/report.py",
    "egfr_pipeline/validate.py",
    "main.py",
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


def check_vina_postprocess_coverage(project_root: Path, config: dict, result: ValidationResult):
    """8.1: Verify Step 4 coverage/provenance spans all expected docking pairs."""
    coverage_rows, _ = load_csv_safe(project_root / "vina_postprocess_coverage.csv")
    if not coverage_rows:
        result.warn("vina_postprocess_coverage.csv not found or empty")
        return

    expected_pairs = len(config.get("receptors", [])) * len(config.get("ligands", []))
    if len(coverage_rows) == expected_pairs:
        result.ok(f"vina_postprocess_coverage.csv covers all {expected_pairs} expected receptor/ligand pairs")
    else:
        result.warn(
            f"vina_postprocess_coverage.csv covers {len(coverage_rows)}/{expected_pairs} expected receptor/ligand pairs"
        )

    problematic = [
        row for row in coverage_rows
        if str(row.get("status", "")).strip().lower() != "parsed"
    ]
    if not problematic:
        result.ok("Vina postprocess coverage shows all expected raw pose files were parsed")
    else:
        labels = [f"{row.get('receptor_id', '?')}/{row.get('ligand_id', '?')}={row.get('status', 'unknown')}" for row in problematic]
        result.warn(f"Vina postprocess coverage has non-parsed inputs: {labels}")


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
    check_vina_postprocess_coverage(project_root, config, result)

    # 8.2
    print("[8.2] CSV schema regression checks...")
    check_csv_schemas(project_root, result)

    # 8.3
    print("[8.3] Residue numbering consistency...")
    check_residue_numbering(config, result)

    # 8.3b
    print("[8.3b] Known mutations check...")
    check_known_mutations(config, result)

    # 8.4
    print("[8.4] Handoff readiness...")
    check_handoff_readiness(repo, result)

    return result
