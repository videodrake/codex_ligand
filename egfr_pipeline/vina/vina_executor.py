#!/usr/bin/env python3
"""
AutoDock Vina Docking Pipeline (Blind / Focused)

Usage:
  python run_docking.py                    # 인터랙티브 모드 (메뉴 선택)
  python run_docking.py --mode blind       # CLI 모드
  python run_docking.py --mode focused --region clobe
  python run_docking.py --mode focused --center -38.2 0.2 -73.8 --box-size 30 30 30
  python run_docking.py --mode blind --ligands input/97806_ligand.pdbqt
  python run_docking.py --mode blind --n-pockets 5
  python run_docking.py --mode blind --n-pockets 5 --max-per-pocket 3 --cluster-radius 8
  python run_docking.py --mode blind --exclude-zone -38.2 0.2 -73.8 10
  python run_docking.py --mode blind --smiles "CCO" "CC(=O)Oc1ccccc1C(=O)O"
  python run_docking.py --config output/2026-02-26_clobe_3gt8_monomer_anp/config.yaml
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
import hashlib
import json
import logging
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger("pipeline.vina")

# standalone 실행 시 handler가 없으면 StreamHandler 자동 추가
if not logging.getLogger("pipeline").handlers and not logger.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    logger.addHandler(_handler)
    logger.setLevel(logging.INFO)
from egfr_pipeline import paths
from egfr_pipeline.config import resolve_resource_config
from egfr_pipeline.runtime import cap_worker_count, resolve_runtime_resources

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# ============================================================
# 스크립트 기준 디렉토리
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent.parent.parent
INPUT_DIR = SCRIPT_DIR / "input"
OUTPUT_DIR = SCRIPT_DIR / "output"
PROJECT_RECEPTOR_IDS = (
    "3GT8_raw",
    "EGFR_160-185",
    "EGFR_170-200",
)

# ============================================================
# Region 프리셋
# ============================================================
REGION_PRESETS = {
    "clobe": {
        "center": [-38.2, 0.2, -73.8],
        "box_size": [30, 30, 30],
        "description": "C-lobe pocket (catalytic loop 808-810 + A-loop PRO848 + AP-2 helix 987-991)",
    },
}

_VINA_CONSTRUCTOR_FALLBACK_WARNED = False
_VINA_WRITE_POSES_FALLBACK_WARNED = False
_VINA_ENERGIES_FALLBACK_WARNED = False


# ============================================================
# Config 저장/로드 (YAML 우선, JSON fallback)
# ============================================================
def save_config(config: dict, path: Path):
    path = Path(path)
    if HAS_YAML:
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    else:
        json_path = path.with_suffix(".json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        logger.info("pyyaml 미설치 → JSON으로 저장: %s", json_path)


def load_config(path: str) -> dict:
    path = Path(path)
    if path.suffix in (".yaml", ".yml"):
        if not HAS_YAML:
            logger.error("YAML config를 로드하려면 pyyaml 설치 필요: pip install pyyaml")
            sys.exit(1)
        with open(path, encoding="utf-8") as f:
            return yaml.safe_load(f)
    else:
        with open(path, encoding="utf-8") as f:
            return json.load(f)


def is_project_config(config: dict) -> bool:
    return isinstance(config.get("receptors"), list) and bool(config.get("receptors"))


def _candidate_project_source_paths(entry: dict, kind: str) -> List[Path]:
    """Return possible source files that can be prepared into the configured PDBQT."""
    pdbqt_path = Path(entry.get("pdbqt", "")) if entry.get("pdbqt") else None
    candidates: List[Path] = []

    if kind == "receptor":
        if entry.get("pdb"):
            candidates.append(Path(entry["pdb"]))
        if pdbqt_path is not None:
            stem = pdbqt_path.stem.replace("_receptor", "")
            candidates.append(pdbqt_path.with_name(f"{stem}.pdb"))
    else:
        for key in ("sdf", "mol2", "pdb"):
            if entry.get(key):
                candidates.append(Path(entry[key]))
        if pdbqt_path is not None:
            stem = pdbqt_path.stem.replace("_ligand", "")
            for suffix in (".sdf", ".mol2", ".pdb"):
                candidates.append(pdbqt_path.with_name(f"{stem}{suffix}"))

    deduped: List[Path] = []
    seen = set()
    for path in candidates:
        resolved = str(path)
        if resolved not in seen:
            deduped.append(path)
            seen.add(resolved)
    return deduped


def _resolve_project_source_path(entry: dict, kind: str) -> Optional[Path]:
    for candidate in _candidate_project_source_paths(entry, kind):
        if candidate.exists():
            return candidate
    return None


def validate_project_config(config: dict) -> None:
    receptors = config.get("receptors")
    ligands = config.get("ligands")
    if not isinstance(receptors, list) or len(receptors) != 3:
        raise ValueError("Project config must define exactly 3 receptors.")
    if not isinstance(ligands, list) or not ligands:
        raise ValueError("Project config must define at least 1 ligand.")

    receptor_ids = [item.get("id") for item in receptors]
    if set(receptor_ids) != set(PROJECT_RECEPTOR_IDS):
        expected = ", ".join(PROJECT_RECEPTOR_IDS)
        raise ValueError(f"Project config receptor ids must be exactly: {expected}")

    for receptor in receptors:
        if not receptor.get("pdbqt"):
            raise ValueError(f"Missing receptor pdbqt path for {receptor.get('id')}")
        pdbqt_path = Path(receptor["pdbqt"])
        pdb_path = receptor.get("pdb")
        if pdb_path and not Path(pdb_path).exists():
            raise FileNotFoundError(f"Receptor pdb file not found: {pdb_path}")
        if not pdbqt_path.exists() and _resolve_project_source_path(receptor, "receptor") is None:
            raise FileNotFoundError(
                f"Receptor pdbqt file not found and no receptor source file is available: {pdbqt_path}"
            )

    for ligand in ligands:
        ligand_id = ligand.get("id")
        ligand_path = ligand.get("pdbqt")
        if not ligand_id or not ligand_path:
            raise ValueError("Each ligand entry must include id and pdbqt.")
        if not Path(ligand_path).exists() and _resolve_project_source_path(ligand, "ligand") is None:
            raise FileNotFoundError(
                f"Ligand pdbqt file not found and no ligand source file is available: {ligand_path}"
            )


def ensure_project_config_inputs_ready(config: dict) -> None:
    """Prepare missing project-config PDBQT inputs from configured source files."""
    for receptor in config.get("receptors", []):
        pdbqt_path = Path(receptor["pdbqt"])
        source_path = _resolve_project_source_path(receptor, "receptor")

        if pdbqt_path.exists() and (source_path is None or not prepared_output_is_stale(source_path, pdbqt_path)):
            continue

        if source_path is None:
            raise FileNotFoundError(
                f"Receptor pdbqt file not found and no receptor source file is available: {pdbqt_path}"
            )

        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        action = "Rebuilding stale receptor PDBQT" if pdbqt_path.exists() else "Preparing receptor PDBQT"
        logger.info(
            "%s for %s: %s -> %s",
            action, receptor.get('id'), source_path, pdbqt_path,
        )
        if not prepare_receptor(source_path, pdbqt_path):
            raise RuntimeError(
                f"Failed to prepare receptor PDBQT for {receptor.get('id')}: {source_path}"
            )

    for ligand in config.get("ligands", []):
        pdbqt_path = Path(ligand["pdbqt"])
        source_path = _resolve_project_source_path(ligand, "ligand")

        if pdbqt_path.exists() and (source_path is None or not prepared_output_is_stale(source_path, pdbqt_path)):
            continue

        if source_path is None:
            raise FileNotFoundError(
                f"Ligand pdbqt file not found and no ligand source file is available: {pdbqt_path}"
            )

        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        action = "Rebuilding stale ligand PDBQT" if pdbqt_path.exists() else "Preparing ligand PDBQT"
        logger.info(
            "%s for %s: %s -> %s",
            action, ligand.get('id'), source_path, pdbqt_path,
        )
        if not prepare_ligand(source_path, pdbqt_path):
            raise RuntimeError(
                f"Failed to prepare ligand PDBQT for {ligand.get('id')}: {source_path}"
            )


def apply_config_to_args(args, config: dict):
    vina_cfg = config.get("vina", {}) if isinstance(config.get("vina"), dict) else {}
    postprocess_cfg = config.get("postprocess", {}) if isinstance(config.get("postprocess"), dict) else {}

    if is_project_config(config):
        validate_project_config(config)
        args.project_name = config.get("project_name")
        args.output_root = config.get("output_root")
        args.max_workers = int(config.get("max_workers", getattr(args, "max_workers", 16)))
        args.receptors = [entry["pdbqt"] for entry in config["receptors"]]
        args.receptor_ids = {entry["pdbqt"]: entry["id"] for entry in config["receptors"]}
        args.receptor_pdb_map = {
            entry["pdbqt"]: entry.get("pdb") for entry in config["receptors"] if entry.get("pdb")
        }
        args.receptor_source_types = {
            entry["pdbqt"]: entry.get("source_type") for entry in config["receptors"]
        }
        args.ligands = [entry["pdbqt"] for entry in config["ligands"]]
        args.ligand_ids = {entry["pdbqt"]: entry["id"] for entry in config["ligands"]}
    else:
        if args.mode is None:
            args.mode = config.get("mode")
        if args.receptor is None:
            args.receptor = config.get("receptor")
        if args.receptor_pdb is None:
            args.receptor_pdb = config.get("receptor_pdb")
        if args.ligands is None:
            args.ligands = config.get("ligands")
        if args.center is None:
            args.center = config.get("center")
        if args.box_size is None:
            args.box_size = config.get("box_size")
        if args.region is None:
            args.region = config.get("region")
        if args.exclude_zone is None:
            args.exclude_zone = config.get("exclude_zone")
        if getattr(args, "n_pockets", None) is None:
            args.n_pockets = config.get("n_pockets")

    if args.mode is None:
        args.mode = config.get("mode", vina_cfg.get("mode"))
    if args.center is None:
        args.center = config.get("center", vina_cfg.get("center"))
    if args.box_size is None:
        args.box_size = config.get("box_size", vina_cfg.get("box_size"))
    if args.region is None:
        args.region = config.get("region", vina_cfg.get("region"))
    if args.exclude_zone is None:
        args.exclude_zone = config.get("exclude_zone", config.get("exclude_zones"))  # singular preferred; plural is legacy alias
    if getattr(args, "n_pockets", None) is None:
        args.n_pockets = config.get("n_pockets", vina_cfg.get("n_pockets"))

    args.max_per_pocket = config.get(
        "max_per_pocket",
        vina_cfg.get("max_per_pocket", getattr(args, "max_per_pocket", 3)),
    )
    args.cluster_radius = config.get(
        "cluster_radius",
        vina_cfg.get("cluster_radius", getattr(args, "cluster_radius", 5.0)),
    )
    args.exhaustiveness = config.get(
        "exhaustiveness", vina_cfg.get("exhaustiveness", args.exhaustiveness)
    )
    args.n_poses = config.get("n_poses", vina_cfg.get("n_poses", args.n_poses))
    args.energy_range = float(
        config.get(
            "energy_range",
            vina_cfg.get("energy_range", getattr(args, "energy_range", 3.0)),
        )
    )
    args.padding = config.get("padding", vina_cfg.get("padding", args.padding))
    args.min_box = config.get("min_box", vina_cfg.get("min_box", args.min_box))
    cpu_value = config.get("cpu", vina_cfg.get("cpu", getattr(args, "cpu", 0)))
    args.cpu = 0 if cpu_value in (None, "") else int(cpu_value)
    args.seed = normalize_optional_int(
        config.get("seed", vina_cfg.get("seed", getattr(args, "seed", None)))
    )
    parallel_value = config.get(
        "parallel_receptors",
        vina_cfg.get(
            "parallel_receptors",
            getattr(args, "parallel_receptors", 1),
        ),
    )
    args.parallel_receptors = int(parallel_value or 1)
    args.max_workers = int(config.get("max_workers", getattr(args, "max_workers", 16)))
    if getattr(args, "smiles", None) is None:
        args.smiles = config.get("smiles")
    args.parse_results = bool(postprocess_cfg.get("parse_results", getattr(args, "parse_results", False)))
    args.extract_contacts = bool(postprocess_cfg.get("extract_contacts", getattr(args, "extract_contacts", False)))
    args.contact_cutoff = float(postprocess_cfg.get("contact_cutoff", getattr(args, "contact_cutoff", 4.0)))
    args.cluster_pockets = bool(postprocess_cfg.get("cluster_pockets", getattr(args, "cluster_pockets", False)))
    args.pocket_cutoff = float(postprocess_cfg.get("pocket_cutoff", getattr(args, "pocket_cutoff", 8.0)))
    args.summarize_pockets = bool(postprocess_cfg.get("summarize_pockets", getattr(args, "summarize_pockets", False)))
    args.compare_pockets = bool(postprocess_cfg.get("compare_pockets", getattr(args, "compare_pockets", False)))
    args.comparison_centroid_cutoff = float(postprocess_cfg.get("comparison_centroid_cutoff", getattr(args, "comparison_centroid_cutoff", 15.0)))
    args.extract_ppi_residues = bool(postprocess_cfg.get("extract_ppi_residues", getattr(args, "extract_ppi_residues", False)))
    args.generate_report = bool(postprocess_cfg.get("generate_report", getattr(args, "generate_report", False)))
    args.validate_outputs = bool(postprocess_cfg.get("validate_outputs", getattr(args, "validate_outputs", False)))
    return args


def resolve_receptor_id(args, receptor_path: Path) -> str:
    receptor_map = getattr(args, "receptor_ids", {}) or {}
    return receptor_map.get(str(receptor_path), receptor_path.stem.replace("_receptor", ""))


def resolve_ligand_id(args, ligand_path: Path) -> str:
    ligand_map = getattr(args, "ligand_ids", {}) or {}
    return ligand_map.get(str(ligand_path), ligand_path.stem.replace("_ligand", ""))


def resolve_receptor_pdb_override(args, receptor_path: Path) -> Optional[Path]:
    receptor_pdb_map = getattr(args, "receptor_pdb_map", {}) or {}
    pdb_path = receptor_pdb_map.get(str(receptor_path))
    return Path(pdb_path) if pdb_path else None


def normalize_optional_int(value: Any) -> Optional[int]:
    if value in (None, "", "auto"):
        return None
    return int(value)


def derive_docking_seed(
    base_seed: Optional[int],
    receptor_id: str,
    ligand_id: str,
) -> Optional[int]:
    if base_seed is None:
        return None

    digest = hashlib.sha256(
        f"{base_seed}:{receptor_id}:{ligand_id}".encode("utf-8")
    ).hexdigest()
    return (int(digest[:8], 16) % 2147483646) + 1


def resolve_ligand_worker_count(
    max_workers: int,
    n_jobs: int,
    *,
    cpu_per_job: int = 0,
    available_cores: Optional[int] = None,
) -> int:
    """Cap ligand-level worker count to avoid oversubscribing visible cores."""
    worker_count = min(int(max_workers), max(int(n_jobs), 1))
    if available_cores is None:
        visible_cores = resolve_runtime_resources().effective_cpus
    else:
        visible_cores = available_cores

    if cpu_per_job <= 0:
        if worker_count > 1:
            logger.warning(
                "cpu=0 allows each Vina ligand job to use all visible cores; "
                "forcing ligand-level parallelism to 1."
            )
        return 1

    capped = cap_worker_count(
        requested_workers=worker_count,
        n_jobs=n_jobs,
        cpu_per_job=cpu_per_job,
        available_cpus=max(1, int(visible_cores or 1)),
    )
    if capped < worker_count:
        logger.warning(
            "Requested ligand CPU budget exceeds visible cores: "
            "%d > %s. Capping ligand workers %d -> %d.",
            worker_count * cpu_per_job, visible_cores, worker_count, capped,
        )
        worker_count = capped

    return max(1, worker_count)


def prepared_output_is_stale(source_path: Path, output_path: Path) -> bool:
    if not output_path.exists():
        return True

    try:
        return source_path.exists() and source_path.stat().st_mtime > output_path.stat().st_mtime
    except OSError:
        return True


def build_receptor_output_dir(args, receptor_id: str) -> Optional[Path]:
    output_root = getattr(args, "output_root", None)
    if not output_root:
        return None
    config = {"output_root": output_root}
    return paths.wa_phase1_vina_receptor(config, receptor_id)


def append_job_status(out_dir: Path, record: Dict[str, Any]) -> None:
    status_path = out_dir / "job_status.tsv"
    file_exists = status_path.exists()
    with open(status_path, "a", encoding="utf-8") as handle:
        if not file_exists:
            handle.write(
                "receptor_id\tligand_id\tstatus\toutput_file\tseed\tcpu\tenergy_range\treturned_poses\terror\n"
            )
        handle.write(
            f"{record['receptor_id']}\t{record['ligand_id']}\t{record['status']}\t"
            f"{record['output_file']}\t{record.get('seed', '')}\t{record.get('cpu', '')}\t"
            f"{record.get('energy_range', '')}\t{record.get('returned_poses', '')}\t"
            f"{record.get('error', '')}\n"
        )


def execute_ligand_job(job: Dict[str, Any]) -> Dict[str, Any]:
    try:
        energies = run_docking(
            job["receptor_path"],
            job["ligand_path"],
            job["output_file"],
            job["center"],
            job["box_size"],
            job["exhaustiveness"],
            job["n_poses"],
            cpu=job.get("cpu"),
            seed=job.get("seed"),
            energy_range=job.get("energy_range", 3.0),
        )
        return {
            "status": "ok",
            "receptor_id": job["receptor_id"],
            "ligand_id": job["ligand_id"],
            "ligand_name": job["ligand_name"],
            "output_file": str(job["output_file"]),
            "energies": energies,
            "seed": job.get("seed"),
            "cpu": job.get("cpu"),
            "energy_range": job.get("energy_range"),
            "returned_poses": len(energies) if energies is not None else 0,
        }
    except Exception as exc:
        return {
            "status": "error",
            "receptor_id": job["receptor_id"],
            "ligand_id": job["ligand_id"],
            "ligand_name": job["ligand_name"],
            "output_file": str(job["output_file"]),
            "seed": job.get("seed"),
            "cpu": job.get("cpu"),
            "energy_range": job.get("energy_range"),
            "returned_poses": 0,
            "error": str(exc),
        }


# ============================================================
# 파일 분류 (receptor / ligand / smiles / unknown)
# ============================================================
def count_atoms(filepath: Path) -> int:
    """PDB/PDBQT/CIF 파일의 ATOM/HETATM 행 수를 센다."""
    count = 0
    try:
        with open(filepath, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    count += 1
    except Exception:
        pass
    return count


def classify_file(filepath: Path) -> str:
    """파일을 receptor/ligand/smiles/unknown 으로 분류.

    분류 우선순위:
      1. 네이밍 규칙 (_receptor, _ligand)
      2. 확장자 (SDF/MOL2 = ligand, SMI = smiles)
      3. 원자 수 (>500 = receptor, ≤500 = ligand)
    """
    ext = filepath.suffix.lower()
    stem = filepath.stem

    # 1) 명시적 네이밍 규칙
    if "_receptor" in stem:
        return "receptor"
    if "_ligand" in stem:
        return "ligand"

    # 2) 확장자 기반
    if ext == ".pml":
        return "style"
    if ext == ".smi":
        return "smiles"
    if ext in (".sdf", ".mol2"):
        return "ligand"  # 소분자 포맷은 항상 ligand

    # 3) PDB/CIF/PDBQT → 원자 수로 판단
    if ext in (".pdb", ".cif", ".pdbqt"):
        atom_count = count_atoms(filepath)
        if atom_count > 500:
            return "receptor"
        elif atom_count > 0:
            return "ligand"

    return "unknown"


# ============================================================
# 자동 감지
# ============================================================
def auto_detect_receptor(input_dir: Path):
    """input/에서 *_receptor.pdbqt 자동 감지. 복수일 경우 리스트 반환."""
    receptors = sorted(input_dir.glob("*_receptor.pdbqt"))
    if len(receptors) == 1:
        logger.info("Receptor: %s", receptors[0].name)
        return receptors
    elif len(receptors) > 1:
        logger.info("Receptor %d개 감지:", len(receptors))
        for r in receptors:
            logger.info("  > %s", r.name)
        return receptors
    else:
        logger.error("input/에서 *_receptor.pdbqt를 찾을 수 없습니다.")
        logger.error("  --receptor로 직접 지정하거나 --prepare-receptor로 변환하세요.")
        sys.exit(1)


def auto_detect_receptor_pdb(receptor_pdbqt: Path, input_dir: Path):
    # receptor PDBQT 이름에서 PDB 유추: xxx_receptor.pdbqt → xxx.pdb
    stem = receptor_pdbqt.name.replace("_receptor.pdbqt", "")
    pdb_candidate = input_dir / f"{stem}.pdb"
    if pdb_candidate.exists():
        logger.info("Receptor PDB: %s", pdb_candidate.name)
        return pdb_candidate

    # 단일 PDB 파일이면 사용
    pdbs = sorted(input_dir.glob("*.pdb"))
    if len(pdbs) == 1:
        logger.info("Receptor PDB: %s", pdbs[0].name)
        return pdbs[0]

    logger.warning("%s에 대응하는 PDB를 찾을 수 없습니다.", receptor_pdbqt.name)
    if pdbs:
        logger.warning("  후보 파일:")
        for p in pdbs:
            logger.warning("    %s", p)
    return None


def auto_detect_ligands(input_dir: Path):
    ligands = sorted(input_dir.glob("*_ligand.pdbqt"))
    if not ligands:
        logger.error("input/에서 *_ligand.pdbqt를 찾을 수 없습니다.")
        sys.exit(1)
    logger.info("Ligands (%d):", len(ligands))
    for lig in ligands:
        logger.info("  %s", lig.name)
    return ligands


# ============================================================
# Receptor PDB → PDBQT 변환
# ============================================================
def prepare_receptor(pdb_file: Path, pdbqt_file: Path):
    pdb_str = str(pdb_file)
    pdbqt_str = str(pdbqt_file)

    commands = [
        ("ADFR", f"prepare_receptor -r {pdb_str} -o {pdbqt_str} -A hydrogens"),
        ("MGLTools", f"pythonsh $(which prepare_receptor4.py) -r {pdb_str} -o {pdbqt_str} -A hydrogens"),
        ("OpenBabel", f"obabel {pdb_str} -O {pdbqt_str} -xr"),
    ]

    for label, cmd in commands:
        logger.debug("Receptor %s 시도: %s", label, cmd)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0 and pdbqt_file.exists():
            logger.info("Receptor %s로 변환 성공 → %s", label, pdbqt_file.name)
            return True
        else:
            logger.debug("Receptor %s 실패: %s", label, result.stderr.strip()[:200])

    logger.error("receptor PDBQT 변환 실패. prepare_receptor 또는 obabel 설치 확인.")
    return False


# ============================================================
# Ligand SDF/MOL2/PDB → PDBQT 변환
# ============================================================
def prepare_ligand(input_file: Path, pdbqt_file: Path):
    in_str = str(input_file)
    out_str = str(pdbqt_file)

    commands = [
        ("Meeko", f"mk_prepare_ligand.py -i {in_str} -o {out_str}"),
        ("ADFR", f"prepare_ligand -l {in_str} -o {out_str}"),
        ("MGLTools", f"pythonsh $(which prepare_ligand4.py) -l {in_str} -o {out_str}"),
        ("OpenBabel", f"obabel {in_str} -O {out_str} --gen3d"),
    ]

    for label, cmd in commands:
        logger.debug("Ligand %s 시도: %s", label, cmd)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0 and pdbqt_file.exists():
            logger.info("Ligand %s로 변환 성공 → %s", label, pdbqt_file.name)
            return True
        else:
            logger.debug("Ligand %s 실패: %s", label, result.stderr.strip()[:200])

    logger.error("ligand PDBQT 변환 실패: %s", input_file.name)
    logger.error("  mk_prepare_ligand (meeko), prepare_ligand, 또는 obabel 설치 확인.")
    return False


# ============================================================
# SMILES → 3D SDF 변환
# ============================================================
def smiles_to_sdf(smiles: str, sdf_file: Path, name: str = "ligand"):
    """SMILES 문자열을 3D 좌표가 포함된 SDF 파일로 변환.

    RDKit (설치 시) 우선 사용, 없으면 OpenBabel fallback.
    """
    # 방법 1: RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            # ETKDGv3 실패 시 랜덤 좌표로 재시도
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        mol.SetProp("_Name", name)
        writer = Chem.SDWriter(str(sdf_file))
        writer.write(mol)
        writer.close()
        logger.info("SMILES RDKit으로 3D 구조 생성 → %s", sdf_file.name)
        return True
    except ImportError:
        pass
    except Exception as e:
        logger.debug("SMILES RDKit 실패: %s", e)

    # 방법 2: OpenBabel
    cmd = f'obabel -:"{smiles}" -O "{sdf_file}" --gen3d --minimize --ff MMFF94'
    logger.debug("SMILES OpenBabel 시도: %s", cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0 and sdf_file.exists() and sdf_file.stat().st_size > 0:
        logger.info("SMILES OpenBabel로 3D 구조 생성 → %s", sdf_file.name)
        return True
    else:
        logger.debug("SMILES OpenBabel 실패: %s", result.stderr.strip()[:200])

    logger.error("SMILES→SDF 변환 실패. RDKit 또는 OpenBabel 설치 필요.")
    return False


# ============================================================
# SMILES 파일 (.smi) 자동 변환
# ============================================================
def auto_convert_smiles_file(input_dir: Path, smi_file: Path):
    """단일 .smi 파일의 모든 SMILES → SDF → PDBQT 변환.

    SMI 파일 포맷 (각 줄):
        SMILES [이름]          # 탭 또는 공백 구분
        # 주석 줄은 무시
        CC(=O)Oc1ccccc1C(=O)O  aspirin
        CCO                              # 이름 없으면 파일명_0, 파일명_1, ...
    """
    lines = smi_file.read_text(encoding="utf-8", errors="ignore").strip().splitlines()
    unnamed_idx = 0
    converted = 0
    failed = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        smiles = parts[0]
        if len(parts) > 1:
            name = parts[1]
        else:
            name = f"{smi_file.stem}_{unnamed_idx}"
            unnamed_idx += 1

        pdbqt_file = input_dir / f"{name}_ligand.pdbqt"
        sdf_file = input_dir / f"{name}_ligand.sdf"

        if pdbqt_file.exists():
            logger.debug("Skip %s → %s (이미 존재)", name, pdbqt_file.name)
            continue

        # Step 1: SMILES → SDF (3D 좌표 생성)
        if not sdf_file.exists():
            logger.debug("SMILES→SDF %s: %s", name, smiles)
            if not smiles_to_sdf(smiles, sdf_file, name):
                failed += 1
                continue

        # Step 2: SDF → PDBQT (기존 prepare_ligand 재사용)
        logger.debug("SDF→PDBQT %s → %s", sdf_file.name, pdbqt_file.name)
        if prepare_ligand(sdf_file, pdbqt_file):
            converted += 1
        else:
            failed += 1

    return converted, failed


def auto_convert_smiles(input_dir: Path):
    """input/ 내 모든 *.smi 파일을 자동 변환. (하위 호환)"""
    smi_files = sorted(input_dir.glob("*.smi"))
    if not smi_files:
        return
    logger.info("SMILES .smi 파일 감지, 자동 변환 시작...")
    total_c, total_f = 0, 0
    for smi_file in smi_files:
        c, f = auto_convert_smiles_file(input_dir, smi_file)
        total_c += c
        total_f += f
    if total_c:
        logger.info("SMILES %d개 변환 완료", total_c)
    if total_f:
        logger.warning("SMILES %d개 변환 실패", total_f)


# ============================================================
# 자동 변환: input/ 내 PDB/SDF → PDBQT
# ============================================================
def auto_convert_receptors(input_dir: Path):
    """input/ 내 *.pdb 중 receptor PDBQT가 없는 것을 자동 변환."""
    pdb_files = sorted(input_dir.glob("*.pdb"))
    if not pdb_files:
        return

    converted = 0
    for pdb in pdb_files:
        pdbqt = input_dir / f"{pdb.stem}_receptor.pdbqt"
        if pdbqt.exists():
            logger.debug("Skip %s → %s (이미 존재)", pdb.name, pdbqt.name)
        else:
            logger.debug("변환 %s → %s", pdb.name, pdbqt.name)
            if prepare_receptor(pdb, pdbqt):
                converted += 1
            else:
                logger.warning("receptor 변환 실패: %s", pdb.name)

    if converted:
        logger.info("Receptor %d개 변환 완료", converted)


def auto_convert_ligands(input_dir: Path):
    """input/ 내 SMILES/SDF/MOL2 → PDBQT 자동 변환."""
    # SMILES → SDF 먼저 처리 (생성된 SDF는 아래에서 PDBQT로 변환)
    auto_convert_smiles(input_dir)

    # 기존 로직: SDF/MOL2 → PDBQT
    source_exts = (".sdf", ".mol2")
    source_files = []
    for ext in source_exts:
        source_files.extend(input_dir.glob(f"*_ligand{ext}"))
    source_files = sorted(set(source_files))

    if not source_files:
        return

    converted = 0
    failed = 0
    for src in source_files:
        # 97806_ligand.sdf → 97806_ligand.pdbqt
        pdbqt = src.with_suffix(".pdbqt")
        if pdbqt.exists():
            logger.debug("Skip %s → %s (이미 존재)", src.name, pdbqt.name)
        else:
            logger.debug("변환 %s → %s", src.name, pdbqt.name)
            if prepare_ligand(src, pdbqt):
                converted += 1
            else:
                failed += 1

    if converted:
        logger.info("Ligand %d개 변환 완료", converted)
    if failed:
        logger.warning("Ligand %d개 변환 실패", failed)


# ============================================================
# 통합 전처리 — 스캔 → 분류 → 변환 → 요약
# ============================================================
def preprocess_inputs(input_dir: Path):
    """input/ 내 모든 분자 파일을 스캔, 분류, 변환하여 PDBQT 준비.

    Returns:
        (ready_receptors, ready_ligands): PDBQT Path 리스트
    """
    ALL_EXTS = {".pdb", ".pdbqt", ".sdf", ".mol2", ".smi", ".cif", ".pml"}
    all_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in ALL_EXTS
    )

    receptors_src = []    # 변환 필요한 receptor 소스
    ligands_src = []      # 변환 필요한 ligand 소스
    smiles_files = []     # SMILES 파일
    style_files = []      # PyMOL 색상 스크립트
    ready_receptors = []  # 이미 PDBQT인 receptor
    ready_ligands = []    # 이미 PDBQT인 ligand
    skipped = []          # 분류 불가

    # ── 1단계: 파일 분류 ──
    for f in all_files:
        cls = classify_file(f)
        ext = f.suffix.lower()

        if cls == "smiles":
            smiles_files.append(f)
        elif cls == "style":
            style_files.append(f)
        elif cls == "receptor":
            if ext == ".pdbqt":
                ready_receptors.append(f)
            else:
                receptors_src.append(f)
        elif cls == "ligand":
            if ext == ".pdbqt":
                ready_ligands.append(f)
            else:
                ligands_src.append(f)
        else:
            skipped.append(f)

    # ── 2단계: SMILES → SDF → PDBQT ──
    for smi_file in smiles_files:
        logger.debug("SMILES %s 변환 중...", smi_file.name)
        auto_convert_smiles_file(input_dir, smi_file)

    # ── 3단계: Receptor 변환 (PDB/MOL2/CIF → PDBQT) ──
    converted_r = 0
    for src in receptors_src:
        stem = src.stem.replace("_receptor", "")
        pdbqt = input_dir / f"{stem}_receptor.pdbqt"
        if pdbqt.exists():
            logger.debug("Skip %s → %s (이미 존재)", src.name, pdbqt.name)
            ready_receptors.append(pdbqt)
        else:
            logger.debug("변환 %s → %s", src.name, pdbqt.name)
            if prepare_receptor(src, pdbqt):
                ready_receptors.append(pdbqt)
                converted_r += 1
            else:
                logger.warning("receptor 변환 실패: %s", src.name)

    # ── 4단계: Ligand 변환 (SDF/MOL2/PDB → PDBQT) ──
    converted_l = 0
    failed_l = 0
    for src in ligands_src:
        stem = src.stem.replace("_ligand", "")
        pdbqt = input_dir / f"{stem}_ligand.pdbqt"
        if pdbqt.exists():
            logger.debug("Skip %s → %s (이미 존재)", src.name, pdbqt.name)
            ready_ligands.append(pdbqt)
        else:
            logger.debug("변환 %s → %s", src.name, pdbqt.name)
            if prepare_ligand(src, pdbqt):
                ready_ligands.append(pdbqt)
                converted_l += 1
            else:
                failed_l += 1

    # SMILES에서 생성된 PDBQT도 포함
    extra_ligands = sorted(input_dir.glob("*_ligand.pdbqt"))
    ready_ligands = sorted(set(ready_ligands + extra_ligands))
    ready_receptors = sorted(set(ready_receptors))

    # ── 5단계: 전처리 요약 ──
    summary_lines = [
        f"전처리 완료 — Receptor: {len(ready_receptors)}개, Ligand: {len(ready_ligands)}개",
    ]
    if converted_r or converted_l:
        summary_lines.append(f"  변환: receptor {converted_r}개, ligand {converted_l}개")
    if failed_l:
        summary_lines.append(f"  실패: ligand {failed_l}개")
    logger.info("\n".join(summary_lines))
    for r in ready_receptors:
        logger.debug("  Receptor > %s", r.name)
    for l in ready_ligands:
        logger.debug("  Ligand   > %s", l.name)
    for sf in style_files:
        logger.debug("  Style    > %s", sf.name)
    for s in skipped:
        logger.debug("  미인식   ? %s", s.name)

    return ready_receptors, ready_ligands


# ============================================================
# ANP 체크
# ============================================================
def check_anp(pdb_file: Path):
    anp_count = 0
    with open(pdb_file, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and "ANP" in line:
                anp_count += 1
    logger.info("ANP receptor에 ANP 원자 %d개 포함", anp_count)
    if anp_count == 0:
        logger.warning("ANP가 없습니다! ATP 포켓 차단 안 됨.")


# ============================================================
# PDB에서 box 중심 + 크기 자동 계산 (blind mode)
# ============================================================
def calc_box_from_pdb(pdb_file: Path, padding: float = 5.0, min_box: float = 70.0):
    coords = []
    with open(pdb_file, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue

    if not coords:
        logger.error("PDB에서 좌표를 읽을 수 없습니다.")
        sys.exit(1)

    coords = np.array(coords)
    center = coords.mean(axis=0).tolist()
    span = coords.max(axis=0) - coords.min(axis=0)
    box_size = (span + 2 * padding).tolist()
    box_size = [max(s, min_box) for s in box_size]

    logger.debug("Box 좌표 범위: %s", span)
    logger.info("Box 중심: [%.1f, %.1f, %.1f]", center[0], center[1], center[2])
    logger.info("Box 크기: [%.1f, %.1f, %.1f]", box_size[0], box_size[1], box_size[2])

    return center, box_size


# ============================================================
# Exclusion Zone: 포즈 필터링
# ============================================================
def parse_poses(pdbqt_path: Path):
    """멀티모델 PDBQT를 개별 포즈(텍스트 + 좌표)로 분리."""
    poses = []
    current_lines = []
    current_coords = []

    with open(pdbqt_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            current_lines.append(line)
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    current_coords.append([x, y, z])
                except ValueError:
                    pass
            if line.startswith("ENDMDL"):
                if current_coords:
                    centroid = np.mean(current_coords, axis=0)
                    poses.append({
                        "lines": current_lines[:],
                        "centroid": centroid,
                    })
                current_lines = []
                current_coords = []

    # 단일 모델 (ENDMDL 없는 경우)
    if current_lines and current_coords:
        centroid = np.mean(current_coords, axis=0)
        poses.append({"lines": current_lines[:], "centroid": centroid})

    return poses


def cluster_poses(poses, cluster_radius=5.0):
    """포즈를 거리 기반 클러스터링. 간단한 leader-follower 방식.

    Returns: list of clusters, 각 cluster는 (center, [pose_indices])
    """
    centroids = np.array([p["centroid"] for p in poses])
    assigned = [False] * len(poses)
    clusters = []

    # 에너지 순서(= 파일 순서)대로 클러스터 리더 선정
    for i in range(len(poses)):
        if assigned[i]:
            continue
        members = [i]
        assigned[i] = True
        for j in range(i + 1, len(poses)):
            if assigned[j]:
                continue
            dist = np.linalg.norm(centroids[i] - centroids[j])
            if dist <= cluster_radius:
                members.append(j)
                assigned[j] = True
        cluster_center = centroids[members].mean(axis=0)
        clusters.append((cluster_center, members))

    return clusters


def discover_pockets(poses, cluster_radius=5.0, max_per_pocket=3, n_pockets=5):
    """포즈를 에너지 순서대로 처리하면서 포켓을 발견/채우기.

    각 포즈를 가장 가까운 열린 포켓에 할당.
    포켓이 max_per_pocket개 차면 닫고 → 다음 포즈는 다른 포켓으로.
    n_pockets개 포켓이 모두 채워지면 종료.

    Returns: list of pocket dicts
    """
    pockets = []  # [{"center": np.array, "pose_indices": [], "closed": bool}]
    skipped = 0

    logger.info("Pocket Discovery 포켓 탐색 (목표 %d개, 포켓당 %d개, 반경 %.1fÅ)",
                n_pockets, max_per_pocket, cluster_radius)

    for i, pose in enumerate(poses):
        centroid = pose["centroid"]

        # 가장 가까운 기존 포켓 찾기
        assigned = None
        min_dist = float("inf")
        for pi, pocket in enumerate(pockets):
            dist = np.linalg.norm(centroid - pocket["center"])
            if dist <= cluster_radius and dist < min_dist:
                min_dist = dist
                assigned = pi

        if assigned is not None:
            pocket = pockets[assigned]
            if pocket["closed"]:
                skipped += 1
                continue  # 이 포켓은 가득 참 → 스킵

            pocket["pose_indices"].append(i)
            # 포켓 중심 업데이트
            pocket_centroids = [poses[j]["centroid"] for j in pocket["pose_indices"]]
            pocket["center"] = np.mean(pocket_centroids, axis=0)

            if len(pocket["pose_indices"]) >= max_per_pocket:
                pocket["closed"] = True
                c = pocket["center"]
                logger.debug("Pocket #%d 완성 (%d개, center=[%.1f, %.1f, %.1f])",
                             assigned+1, max_per_pocket, c[0], c[1], c[2])

                # 완성된 포켓이 목표 수에 도달하면 종료
                completed_count = sum(1 for p in pockets if p["closed"])
                if completed_count >= n_pockets:
                    logger.info("목표 %d개 포켓 모두 완성!", n_pockets)
                    break
        else:
            # 새 포켓 발견 (목표 수 이내만 허용)
            if len(pockets) >= n_pockets:
                skipped += 1
                continue  # 이미 충분한 포켓이 있으므로 스킵

            new_pocket = {
                "center": centroid.copy(),
                "pose_indices": [i],
                "closed": 1 >= max_per_pocket,
            }
            pockets.append(new_pocket)
            num = len(pockets)
            logger.debug("Pocket #%d 발견! (pose #%d, center=[%.1f, %.1f, %.1f])",
                         num, i+1, centroid[0], centroid[1], centroid[2])

            if new_pocket["closed"]:
                logger.debug("Pocket #%d 즉시 완성 (%d개)", num, max_per_pocket)
                completed_count = sum(1 for p in pockets if p["closed"])
                if completed_count >= n_pockets:
                    logger.info("목표 %d개 포켓 모두 완성!", n_pockets)
                    break

    # 결과 요약
    completed = [p for p in pockets if p["closed"]]
    partial = [p for p in pockets if not p["closed"]]

    logger.info("Pocket Discovery 결과: 총 포즈 %d개 | 스킵 %d개 | "
                "발견 포켓 %d개 (완성 %d, 미완성 %d)",
                len(poses), skipped, len(pockets), len(completed), len(partial))
    for pi, pocket in enumerate(pockets):
        c = pocket["center"]
        status = "완성" if pocket["closed"] else f"미완성({len(pocket['pose_indices'])}개)"
        best_idx = pocket["pose_indices"][0]
        logger.info("  Pocket #%d: [%.1f, %.1f, %.1f] | %s | best pose #%d",
                     pi+1, c[0], c[1], c[2], status, best_idx+1)

    return pockets


def save_pockets_info(pockets, poses, energies, lig_name, out_dir):
    """포켓 탐색 결과를 pockets.yaml에 누적 저장."""
    pockets_file = out_dir / "pockets.yaml"

    # 기존 데이터 로드 (여러 ligand 결과 누적)
    existing = {}
    if pockets_file.exists():
        existing = load_config(str(pockets_file)) or {}
    # JSON fallback 대응
    json_fallback = pockets_file.with_suffix(".json")
    if not existing and json_fallback.exists():
        existing = load_config(str(json_fallback)) or {}

    pocket_list = []
    for pi, pocket in enumerate(pockets):
        c = pocket["center"]
        best_idx = pocket["pose_indices"][0]
        best_energy = float(energies[best_idx][0]) if best_idx < len(energies) else None
        pocket_list.append({
            "id": pi + 1,
            "center": [round(float(c[0]), 1), round(float(c[1]), 1), round(float(c[2]), 1)],
            "n_poses": len(pocket["pose_indices"]),
            "completed": pocket["closed"],
            "best_pose": best_idx + 1,
            "best_energy": best_energy,
        })

    existing[lig_name] = pocket_list
    save_config(existing, pockets_file)


def load_previous_pockets(output_dir: Path):
    """output/ 하위에서 pockets.yaml 파일들을 찾아 포켓 정보 로드.

    Returns:
        list of dict: [{
            "dir": output_subdir_name,
            "ligand": ligand_name,
            "id": pocket_id,
            "center": [x, y, z],
            "best_energy": float,
            "n_poses": int,
            "completed": bool,
        }, ...]
    """
    results = []
    if not output_dir.exists():
        return results
    # YAML / JSON 둘 다 검색 (JSON fallback 대응)
    pocket_files = sorted(output_dir.rglob("pockets.yaml"))
    pocket_files += sorted(output_dir.rglob("pockets.json"))
    seen_dirs = set()
    for pf in pocket_files:
        # 같은 디렉토리에 yaml과 json 둘 다 있으면 첫 번째만
        if pf.parent in seen_dirs:
            continue
        seen_dirs.add(pf.parent)
        data = load_config(str(pf))
        if not data or not isinstance(data, dict):
            continue
        subdir = pf.parent.name
        for lig_name, pocket_list in data.items():
            if not isinstance(pocket_list, list):
                continue
            for pocket in pocket_list:
                if not isinstance(pocket, dict) or "center" not in pocket:
                    continue
                results.append({
                    "dir": subdir,
                    "ligand": lig_name,
                    "id": pocket.get("id", "?"),
                    "center": pocket["center"],
                    "best_energy": pocket.get("best_energy"),
                    "n_poses": pocket.get("n_poses", 0),
                    "completed": pocket.get("completed", False),
                })
    return results


def filter_poses_by_exclusion(poses, exclude_zones):
    """exclusion zone 내의 포즈를 제거 (수동 지정).

    exclude_zones: [(center_x, center_y, center_z, radius), ...]
    """
    kept = []
    excluded = 0
    for i, pose in enumerate(poses):
        centroid = pose["centroid"]
        in_zone = False
        for zone in exclude_zones:
            zc = np.array(zone[:3])
            zr = zone[3]
            dist = np.linalg.norm(centroid - zc)
            if dist <= zr:
                in_zone = True
                logger.debug("Exclude pose %d: centroid [%.1f, %.1f, %.1f], "
                             "zone center [%.1f, %.1f, %.1f]까지 %.1fÅ (반경 %.0fÅ 이내)",
                             i+1, centroid[0], centroid[1], centroid[2],
                             zc[0], zc[1], zc[2], dist, zr)
                break
        if in_zone:
            excluded += 1
        else:
            kept.append(pose)

    logger.info("Filter %d개 포즈 중 %d개 제외 → %d개 남음", len(poses), excluded, len(kept))
    return kept


def write_filtered_poses(poses, output_path: Path):
    """필터링된 포즈를 PDBQT로 저장."""
    with open(output_path, "w", encoding="utf-8") as f:
        for pose in poses:
            for line in pose["lines"]:
                f.write(line)


# ============================================================
# 도킹 실행
# ============================================================
def run_docking(
    receptor_pdbqt,
    ligand_pdbqt,
    output_pdbqt,
    center,
    box_size,
    exhaustiveness=128,
    n_poses=20,
    cpu=None,
    seed=None,
    energy_range=3.0,
):
    from vina import Vina

    global _VINA_CONSTRUCTOR_FALLBACK_WARNED
    global _VINA_WRITE_POSES_FALLBACK_WARNED
    global _VINA_ENERGIES_FALLBACK_WARNED

    vina_kwargs = {"sf_name": "vina"}
    if cpu is not None:
        vina_kwargs["cpu"] = int(cpu)
    if seed is not None:
        vina_kwargs["seed"] = int(seed)

    try:
        v = Vina(**vina_kwargs)
    except TypeError:
        if (cpu is not None or seed is not None) and not _VINA_CONSTRUCTOR_FALLBACK_WARNED:
            logger.warning(
                "Installed vina Python API does not expose cpu/seed constructor arguments; "
                "falling back to default Vina runtime settings."
            )
            _VINA_CONSTRUCTOR_FALLBACK_WARNED = True
        v = Vina(sf_name="vina")

    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    try:
        v.write_poses(
            str(output_pdbqt),
            n_poses=n_poses,
            energy_range=energy_range,
            overwrite=True,
        )
    except TypeError:
        if energy_range is not None and not _VINA_WRITE_POSES_FALLBACK_WARNED:
            logger.warning(
                "Installed vina Python API does not expose energy_range on write_poses; "
                "falling back to the library default pose filtering."
            )
            _VINA_WRITE_POSES_FALLBACK_WARNED = True
        v.write_poses(str(output_pdbqt), n_poses=n_poses, overwrite=True)

    try:
        energies = v.energies(n_poses=n_poses, energy_range=energy_range)
    except TypeError:
        if energy_range is not None and not _VINA_ENERGIES_FALLBACK_WARNED:
            logger.warning(
                "Installed vina Python API does not expose energy_range on energies(); "
                "returned pose metadata may reflect the library default filtering."
            )
            _VINA_ENERGIES_FALLBACK_WARNED = True
        energies = v.energies(n_poses=n_poses)
    return energies


# ============================================================
# 결과 출력
# ============================================================
def print_results(name, mode, energies, center, box_size, display_n=None):
    if not energies:
        logger.info(f"{name} {mode.capitalize()} Docking: no poses returned")
        return
    show_n = min(display_n or len(energies), len(energies))
    lines = [
        f"{name} {mode.capitalize()} Docking Results "
        f"Center=[{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}] "
        f"Box=[{box_size[0]:.1f}, {box_size[1]:.1f}, {box_size[2]:.1f}]",
    ]
    for i in range(show_n):
        row = energies[i]
        lines.append(f"  Mode {i+1:>3}: {row[0]:>10.2f}")
    if show_n < len(energies):
        lines.append(f"  ... (총 {len(energies)}개 포즈 중 상위 {show_n}개 표시)")
    logger.info("\n".join(lines))


def print_summary(all_results, out_dir):
    lines = ["Summary - Best Affinity per Ligand"]
    for name, energies in all_results.items():
        if not energies:
            lines.append(f"  {name:>12} {'N/A':>15}")
            continue
        best = energies[0][0]
        lines.append(f"  {name:>12} {best:>15.2f} kcal/mol")
    lines.append(f"결과 저장: {out_dir}")
    logger.info("\n".join(lines))


# ============================================================
# PyMOL 시각화 스크립트 생성
# ============================================================
# 리간드 색상 팔레트 (탄소 색상 util 함수)
LIGAND_COLORS = [
    ("cyan",    "util.cbac"),
    ("magenta", "util.cbam"),
    ("yellow",  "util.cbay"),
    ("salmon",  "util.cbas"),
    ("green",   "util.cbag"),
    ("orange",  "util.cbao"),
    ("pink",    "util.cbap"),
    ("wheat",   "util.cbaw"),
    ("slate",   "color slate,"),
    ("lime",    "color lime,"),
]


def find_style_pml(input_dir: Path, receptor_path: Path):
    """input/에서 receptor에 매칭되는 .pml 색상 스크립트 찾기."""
    pml_files = sorted(input_dir.glob("*.pml"))
    if not pml_files:
        return None
    if len(pml_files) == 1:
        return pml_files[0]
    # receptor stem으로 매칭
    r_stem = receptor_path.stem.replace("_receptor", "")
    for pml in pml_files:
        if r_stem in pml.stem or pml.stem in r_stem:
            return pml
    return pml_files[0]


def generate_pymol_script(out_dir, receptor_pdb, receptor_pdbqt,
                          ligand_results, mode, style_pml=None):
    """도킹 결과 PyMOL 시각화 스크립트(.pml) 생성.

    Args:
        out_dir: 출력 디렉토리
        receptor_pdb: receptor PDB 경로 (잔기 색칠용)
        receptor_pdbqt: receptor PDBQT 경로 (PDB 없을 때 대체)
        ligand_results: {lig_name: output_pdbqt_path} 딕셔너리
        mode: "blind" or "focused"
        style_pml: 사용자 색상 .pml 파일 경로 (optional)
    """
    from datetime import datetime as _dt

    # receptor 로드할 파일 결정 (PDB > PDBQT)
    receptor_file = receptor_pdb if receptor_pdb and Path(receptor_pdb).exists() else receptor_pdbqt
    receptor_file = Path(receptor_file)
    receptor_label = receptor_file.stem.replace("_receptor", "")

    lines = []
    lines.append("# " + "=" * 60)
    lines.append("# AutoDock Vina Docking Results Viewer (PyMOL)")
    lines.append(f"# Generated: {_dt.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"# Receptor: {receptor_label}")
    lines.append(f"# Mode: {mode}")
    lines.append("# " + "=" * 60)
    lines.append("")

    # ── 1. Receptor 로드 ──
    lines.append("# ── 1. Receptor ──")
    rpath = str(receptor_file).replace("\\", "/")
    lines.append(f'load {rpath}, receptor')
    lines.append("show cartoon, receptor")
    lines.append("hide lines, receptor")
    lines.append("color gray, receptor")
    lines.append("")

    # ── 2. 사용자 색상 적용 ──
    if style_pml and Path(style_pml).exists():
        lines.append("# ── 2. User Style ──")
        try:
            with open(style_pml, encoding="utf-8", errors="ignore") as f:
                user_lines = f.read().strip()
            lines.append(f"# source: {Path(style_pml).name}")
            lines.append(user_lines)
        except Exception:
            lines.append(f"# [WARNING] 색상 스크립트 읽기 실패: {style_pml}")
        lines.append("")
    else:
        lines.append("# ── 2. User Style (없음 — input/에 .pml 추가 가능) ──")
        lines.append("")

    # ── 3. 도킹 결과 로드 ──
    lines.append("# ── 3. Docking Results ──")
    lines.append("set stick_radius, 0.15")
    lines.append("")

    lig_names = []
    for i, (lig_name, output_pdbqt) in enumerate(ligand_results.items()):
        output_pdbqt = Path(output_pdbqt)
        if not output_pdbqt.exists():
            continue

        obj_name = lig_name.replace("-", "_").replace(" ", "_")
        lig_names.append(obj_name)
        lpath = str(output_pdbqt).replace("\\", "/")
        color_label, color_cmd = LIGAND_COLORS[i % len(LIGAND_COLORS)]

        lines.append(f"# {lig_name} ({color_label})")
        lines.append(f"load {lpath}, {obj_name}")
        lines.append(f"show sticks, {obj_name}")
        lines.append(f"hide cartoon, {obj_name}")
        lines.append(f"{color_cmd} {obj_name}")
        lines.append("")

    # ── 4. 포켓/필터 결과 ──
    pocket_loaded = False
    for lig_name in list(ligand_results.keys()):
        obj_name = lig_name.replace("-", "_").replace(" ", "_")

        # 완성 포켓 best pose
        pockets_file = out_dir / f"{lig_name}_{mode}_pockets.pdbqt"
        if pockets_file.exists():
            if not pocket_loaded:
                lines.append("# ── 4. Pocket Results ──")
                pocket_loaded = True
            pkt_name = f"{obj_name}_pockets"
            ppath = str(pockets_file).replace("\\", "/")
            lines.append(f"load {ppath}, {pkt_name}")
            lines.append(f"show sticks, {pkt_name}")
            lines.append(f"hide cartoon, {pkt_name}")
            lines.append(f"util.cbag {pkt_name}")
            lines.append(f"set stick_radius, 0.2, {pkt_name}")
            lines.append("")

        # 미완성 포켓 (비활성)
        partial_file = out_dir / f"{lig_name}_{mode}_partial_pockets.pdbqt"
        if partial_file.exists():
            prt_name = f"{obj_name}_partial"
            ppath = str(partial_file).replace("\\", "/")
            lines.append(f"load {ppath}, {prt_name}")
            lines.append(f"show sticks, {prt_name}")
            lines.append(f"color white, {prt_name}")
            lines.append(f"disable {prt_name}")
            lines.append("")

        # 필터 결과
        filtered_file = out_dir / f"{lig_name}_{mode}_filtered.pdbqt"
        if filtered_file.exists():
            if not pocket_loaded:
                lines.append("# ── 4. Filtered Results ──")
                pocket_loaded = True
            flt_name = f"{obj_name}_filtered"
            fpath = str(filtered_file).replace("\\", "/")
            lines.append(f"load {fpath}, {flt_name}")
            lines.append(f"show sticks, {flt_name}")
            lines.append(f"hide cartoon, {flt_name}")
            lines.append(f"util.cbag {flt_name}")
            lines.append("")

    # ── 5. 화면 설정 ──
    lines.append("# ── 5. Display Settings ──")
    lines.append("bg_color black")
    lines.append("set ray_shadow, 0")
    lines.append("set cartoon_transparency, 0.3")
    lines.append("set label_size, 20")
    lines.append("zoom receptor")
    lines.append("")

    # 파일 저장
    pml_out = out_dir / "view_results.pml"
    with open(pml_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    logger.info("PyMOL 시각화 스크립트 생성: %s", pml_out.name)


# ============================================================
# 출력 디렉토리 생성
# ============================================================
def setup_output_dir(args) -> Path:
    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        date_str = datetime.now().strftime("%Y-%m-%d")
        # label: receptor 이름에서 추출
        if args.label:
            label = args.label
        else:
            label = Path(args.receptor).name.replace("_receptor.pdbqt", "")
            if args.mode == "focused" and args.region:
                label = f"{args.region}_{label}"
        out_dir = OUTPUT_DIR / f"{date_str}_{label}"

    out_dir.mkdir(parents=True, exist_ok=True)

    # config 저장
    config = {
        "mode": args.mode,
        "receptor": str(args.receptor),
        "receptor_id": getattr(args, "receptor_id", None),
        "receptor_source_type": getattr(args, "receptor_source_type", None),
        "receptor_pdb": str(args.receptor_pdb) if args.receptor_pdb else None,
        "ligands": [str(l) for l in args.ligands],
        "center": args.center,
        "box_size": args.box_size,
        "exhaustiveness": args.exhaustiveness,
        "n_poses": args.n_poses,
        "energy_range": getattr(args, "energy_range", 3.0),
        "padding": args.padding,
        "min_box": args.min_box,
        "cpu": getattr(args, "cpu", 0),
        "seed": getattr(args, "seed", None),
        "parallel_receptors": getattr(args, "parallel_receptors", None),
        "region": args.region,
        "exclude_zone": args.exclude_zone,
        "n_pockets": getattr(args, "n_pockets", None),
        "max_per_pocket": getattr(args, "max_per_pocket", 3),
        "cluster_radius": getattr(args, "cluster_radius", 5.0),
        "smiles": getattr(args, "smiles", None),
        "project_name": getattr(args, "project_name", None),
        "output_root": getattr(args, "output_root", None),
        "max_workers": getattr(args, "max_workers", 16),
    }
    save_config(config, out_dir / "config.yaml")

    return out_dir


# ============================================================
# argparse
# ============================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="AutoDock Vina Docking Pipeline (Blind / Focused)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_docking.py --mode blind
  python run_docking.py --mode focused --region clobe
  python run_docking.py --mode focused --center -38.2 0.2 -73.8 --box-size 30 30 30
  python run_docking.py --config output/2026-02-26_clobe_3gt8_monomer_anp/config.yaml

Available regions: """ + ", ".join(REGION_PRESETS.keys()),
    )

    # Mode
    parser.add_argument("--mode", choices=["blind", "focused"],
                        help="Docking mode: 'blind' (whole protein) or 'focused' (specific pocket)")

    # Input files
    parser.add_argument("--receptor", type=str, default=None,
                        help="Receptor PDBQT path (default: auto-detect in input/)")
    parser.add_argument("--receptor-pdb", type=str, default=None,
                        help="Receptor PDB path (required for blind mode box calculation)")
    parser.add_argument("--ligands", nargs="+", default=None,
                        help="Ligand PDBQT file(s) (default: all *_ligand.pdbqt in input/)")
    parser.add_argument("--smiles", nargs="+", default=None,
                        metavar="SMILES",
                        help="SMILES 문자열 직접 입력 (자동 3D→PDBQT 변환 후 도킹)")

    # Box parameters
    parser.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"),
                        help="Box center coordinates (required for focused mode without --region)")
    parser.add_argument("--box-size", type=float, nargs=3, metavar=("X", "Y", "Z"),
                        default=None,
                        help="Box dimensions in Angstroms (default: 30 30 30 for focused)")
    parser.add_argument("--region", type=str, default=None,
                        choices=list(REGION_PRESETS.keys()),
                        help="Named region preset (overrides --center/--box-size)")

    # Docking parameters
    parser.add_argument("--exhaustiveness", type=int, default=128,
                        help="Search exhaustiveness (default: 128)")
    parser.add_argument("--n-poses", type=int, default=20,
                        help="Number of poses to generate (default: 20)")
    parser.add_argument("--energy-range", type=float, default=3.0,
                        help="Energy window retained in Vina output poses (default: 3.0)")
    parser.add_argument("--cpu", type=int, default=0,
                        help="CPU threads per Vina run (0 uses all visible cores)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Base random seed for deterministic per-job seeds")
    parser.add_argument("--padding", type=float, default=5.0,
                        help="Box padding for blind mode in Angstroms (default: 5.0)")
    parser.add_argument("--min-box", type=float, default=70.0,
                        help="Minimum box dimension for blind mode (default: 70.0)")

    # Output
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: output/<mode>/<date>_<label>/)")
    parser.add_argument("--label", type=str, default=None,
                        help="Label for run directory (default: derived from receptor name)")

    # Config
    parser.add_argument("--config", type=str, default=None,
                        help="YAML/JSON config file (overrides other args)")
    parser.add_argument("--max-workers", type=int, default=16,
                        help="Maximum ligand workers per receptor (default: 16)")
    parser.add_argument("--parse-results", action="store_true",
                        help="Parse Vina outputs into vina_pose_table.csv after docking")
    parser.add_argument("--extract-contacts", action="store_true",
                        help="Extract receptor contact residues after pose parsing")
    parser.add_argument("--contact-cutoff", type=float, default=4.0,
                        help="Heavy-atom contact cutoff in Angstrom (default: 4.0)")
    parser.add_argument("--cluster-pockets", action="store_true",
                        help="Assign receptor-level pocket ids after pose parsing")
    parser.add_argument("--pocket-cutoff", type=float, default=4.0,
                        help="Pocket clustering cutoff in Angstrom (default: 4.0)")
    parser.add_argument("--summarize-pockets", action="store_true",
                        help="Generate vina_pocket_table.csv and vina_drug_pocket_map.csv")
    parser.add_argument("--compare-pockets", action="store_true",
                        help="Cross-receptor pocket comparison (vina_pocket_comparison.csv)")
    parser.add_argument("--comparison-centroid-cutoff", type=float, default=15.0,
                        help="Max centroid distance for pocket comparison (default: 15.0 Å)")
    parser.add_argument("--extract-ppi-residues", action="store_true",
                        help="Extract PPI interface residues from PyRosetta results")
    parser.add_argument("--generate-report", action="store_true",
                        help="Generate project report after all postprocessing")
    parser.add_argument("--validate-outputs", action="store_true",
                        help="Run output validation checks")

    # Exclusion zone (수동)
    parser.add_argument("--exclude-zone", nargs=4, type=float, action="append",
                        metavar=("X", "Y", "Z", "RADIUS"),
                        help="Exclude poses within RADIUS Å of (X,Y,Z). 여러 번 사용 가능.")

    # Pocket discovery (자동 다중 포켓 탐색)
    parser.add_argument("--n-pockets", type=int, default=None, metavar="N",
                        help="다중 포켓 탐색: N개 포켓을 찾을 때까지 반복 (기본: 미사용)")
    parser.add_argument("--max-per-pocket", type=int, default=3,
                        help="포켓당 최대 포즈 수 (가득 차면 닫힘, 기본 3)")
    parser.add_argument("--cluster-radius", type=float, default=5.0,
                        help="같은 포켓 판정 반경 (기본 5.0 Å)")

    # Receptor preparation
    parser.add_argument("--prepare-receptor", action="store_true",
                        help="Prepare receptor PDBQT from PDB before docking")

    return parser.parse_args()


# ============================================================
# 인터랙티브 모드 헬퍼
# ============================================================
def ask_choice(prompt, options, allow_quit=True):
    """번호로 단일 선택. options: [(label, value), ...]"""
    print(f"\n--- {prompt} ---")
    for i, (label, _) in enumerate(options, 1):
        print(f"  [{i}] {label}")
    if allow_quit:
        print(f"  [q] 종료")

    while True:
        raw = input("\n선택: ").strip().lower()
        if allow_quit and raw == "q":
            print("종료합니다.")
            sys.exit(0)
        try:
            idx = int(raw) - 1
            if 0 <= idx < len(options):
                return options[idx][1]
        except ValueError:
            pass
        print(f"  1~{len(options)} 사이의 번호를 입력하세요.")


def ask_multi_choice(prompt, options):
    """번호로 복수 선택. 콤마 구분 또는 'all'. options: [(label, value), ...]"""
    print(f"\n--- {prompt} (복수 선택: 1,2 또는 all) ---")
    for i, (label, _) in enumerate(options, 1):
        print(f"  [{i}] {label}")

    while True:
        raw = input("\n선택: ").strip().lower()
        if raw == "all":
            return [v for _, v in options]
        try:
            indices = [int(x.strip()) - 1 for x in raw.split(",")]
            if all(0 <= i < len(options) for i in indices) and indices:
                return [options[i][1] for i in indices]
        except ValueError:
            pass
        print(f"  번호를 콤마로 구분하거나 'all'을 입력하세요. (예: 1,2)")


def ask_value(prompt, default=None, cast=None):
    """값 입력. Enter 시 기본값 사용."""
    if default is not None:
        raw = input(f"{prompt} [기본 {default}]: ").strip()
    else:
        raw = input(f"{prompt}: ").strip()

    if not raw and default is not None:
        return default
    if not raw:
        print("  값을 입력하세요.")
        return ask_value(prompt, default, cast)
    if cast:
        try:
            return cast(raw)
        except ValueError:
            print(f"  올바른 형식으로 입력하세요.")
            return ask_value(prompt, default, cast)
    return raw


def ask_confirm(prompt, default=True):
    """Y/N 확인."""
    suffix = "(Y/n)" if default else "(y/N)"
    raw = input(f"{prompt} {suffix}: ").strip().lower()
    if not raw:
        return default
    return raw in ("y", "yes")


def interactive_mode():
    """인터랙티브 메뉴로 설정을 수집하여 argparse.Namespace 반환."""
    print("=" * 60)
    print("  AutoDock Vina Docking Pipeline")
    print("  인터랙티브 모드")
    print("=" * 60)

    # Step 1: 모드 선택
    mode_choice = ask_choice("모드 선택", [
        ("Blind Docking    (단백질 전체 탐색)", "blind"),
        ("Focused Docking  (특정 포켓 집중)", "focused"),
        ("Config 재실행    (이전 설정 로드)", "config"),
    ])

    # Config 재실행
    if mode_choice == "config":
        configs = sorted(OUTPUT_DIR.rglob("config.*"))
        if not configs:
            print("[ERROR] output/ 에서 config 파일을 찾을 수 없습니다.")
            sys.exit(1)
        config_options = [
            (str(c.relative_to(SCRIPT_DIR)), c) for c in configs
        ]
        config_path = ask_choice("Config 파일 선택", config_options, allow_quit=True)
        config = load_config(str(config_path))

        # config 내용 표시
        print(f"\n--- 로드된 설정 ---")
        for k, v in config.items():
            print(f"  {k}: {v}")

        if not ask_confirm("\n이 설정으로 실행하시겠습니까?"):
            print("취소합니다.")
            sys.exit(0)

        return argparse.Namespace(
            mode=config.get("mode"),
            receptor=config.get("receptor"),
            receptors=None,  # config 재실행은 단일 receptor
            receptor_pdb=config.get("receptor_pdb"),
            ligands=config.get("ligands"),
            center=config.get("center"),
            box_size=config.get("box_size"),
            region=config.get("region"),
            exclude_zone=config.get("exclude_zone"),
            n_pockets=config.get("n_pockets"),
            max_per_pocket=config.get("max_per_pocket", 3),
            cluster_radius=config.get("cluster_radius", 5.0),
            exhaustiveness=config.get("exhaustiveness", 128),
            n_poses=config.get("n_poses", 20),
            energy_range=config.get("energy_range", 3.0),
            padding=config.get("padding", 5.0),
            min_box=config.get("min_box", 70.0),
            cpu=config.get("cpu", 0),
            seed=normalize_optional_int(config.get("seed")),
            parallel_receptors=config.get("parallel_receptors"),
            output_dir=None,
            label=None,
            config=None,
            prepare_receptor=False,
            smiles=config.get("smiles"),
        )

    mode = mode_choice

    # 통합 전처리 (스캔 → 분류 → 변환 → 요약)
    ready_receptors, ready_ligands = preprocess_inputs(INPUT_DIR)

    # Step 2: Receptor 선택 (복수 허용)
    if not ready_receptors:
        print("[ERROR] receptor 파일을 찾을 수 없습니다.")
        print("  input/에 PDB, MOL2, CIF, 또는 PDBQT 파일을 넣어주세요.")
        sys.exit(1)
    if len(ready_receptors) == 1:
        receptors = ready_receptors
        print(f"\n[Auto] Receptor: {receptors[0].name}")
    else:
        print(f"\n--- Receptor 선택 ---")
        for i, r in enumerate(ready_receptors):
            print(f"  [{i+1}] {r.name}")
        print(f"  [a] 전체 선택 ({len(ready_receptors)}개 × 모든 ligand)")
        choice = input("선택 (번호/a, 쉼표로 복수 가능): ").strip().lower()
        if choice == "a":
            receptors = ready_receptors
        else:
            try:
                indices = [int(x.strip()) - 1 for x in choice.split(",")]
                receptors = [ready_receptors[i] for i in indices
                             if 0 <= i < len(ready_receptors)]
            except (ValueError, IndexError):
                receptors = [ready_receptors[0]]
        if not receptors:
            print("[ERROR] 선택된 receptor가 없습니다.")
            sys.exit(1)
        for r in receptors:
            print(f"  > {r.name}")

    # Step 3: Ligand 선택
    ligands_all = ready_ligands  # preprocess_inputs()에서 반환된 목록 사용

    # SMILES 직접 입력 옵션
    while True:
        if ligands_all:
            print(f"\n--- Ligand 선택 ---")
            for i, l in enumerate(ligands_all):
                print(f"  [{i+1}] {l.name}")
            print(f"  [s] SMILES 직접 입력")
            print(f"  [a] 전체 선택")
            choice = input("선택 (번호/s/a): ").strip().lower()

            if choice == "s":
                smiles_input = input("SMILES 입력 (예: CCO ethanol, 쉼표로 복수 가능): ").strip()
                for entry in smiles_input.split(","):
                    entry = entry.strip()
                    if not entry:
                        continue
                    parts = entry.split()
                    smi_str = parts[0]
                    smi_name = parts[1] if len(parts) > 1 else f"interactive_{len(ligands_all)}"
                    sdf_f = INPUT_DIR / f"{smi_name}_ligand.sdf"
                    pdbqt_f = INPUT_DIR / f"{smi_name}_ligand.pdbqt"
                    if not pdbqt_f.exists():
                        if not sdf_f.exists():
                            smiles_to_sdf(smi_str, sdf_f, smi_name)
                        if sdf_f.exists():
                            prepare_ligand(sdf_f, pdbqt_f)
                    if pdbqt_f.exists():
                        ligands_all.append(pdbqt_f)
                ligands_all = sorted(set(ligands_all))
                continue  # 다시 선택 메뉴로

            elif choice == "a":
                ligands = ligands_all
                break
            else:
                # 번호 선택 (쉼표 구분)
                try:
                    indices = [int(x.strip()) - 1 for x in choice.split(",")]
                    ligands = [ligands_all[i] for i in indices if 0 <= i < len(ligands_all)]
                    if ligands:
                        break
                except (ValueError, IndexError):
                    pass
                print("잘못된 입력. 다시 선택하세요.")
        else:
            print(f"\n--- Ligand 없음 ---")
            print(f"  [s] SMILES 직접 입력")
            choice = input("선택: ").strip().lower()
            if choice == "s":
                smiles_input = input("SMILES 입력 (예: CCO ethanol, 쉼표로 복수 가능): ").strip()
                for entry in smiles_input.split(","):
                    entry = entry.strip()
                    if not entry:
                        continue
                    parts = entry.split()
                    smi_str = parts[0]
                    smi_name = parts[1] if len(parts) > 1 else f"interactive_{len(ligands_all)}"
                    sdf_f = INPUT_DIR / f"{smi_name}_ligand.sdf"
                    pdbqt_f = INPUT_DIR / f"{smi_name}_ligand.pdbqt"
                    if not pdbqt_f.exists():
                        if not sdf_f.exists():
                            smiles_to_sdf(smi_str, sdf_f, smi_name)
                        if sdf_f.exists():
                            prepare_ligand(sdf_f, pdbqt_f)
                    if pdbqt_f.exists():
                        ligands_all.append(pdbqt_f)
                ligands_all = sorted(set(ligands_all))
            else:
                print("[ERROR] ligand가 없습니다.")
                sys.exit(1)

    if not ligands:
        print("[ERROR] 선택된 ligand가 없습니다.")
        sys.exit(1)
    for l in ligands:
        print(f"  > {l.name}")

    # Step 4: Focused 모드 추가 설정
    region = None
    center = None
    box_size = None

    if mode == "focused":
        region_options = [
            (f"{name} — {info['description']}", name)
            for name, info in REGION_PRESETS.items()
        ]

        # 이전 포켓 결과 확인
        prev_pockets = load_previous_pockets(OUTPUT_DIR)
        if prev_pockets:
            region_options.append(("이전 결과에서 포켓 선택", "__prev_pocket__"))

        region_options.append(("직접 좌표 입력", "__custom__"))

        region_choice = ask_choice("Docking 영역 선택", region_options, allow_quit=True)

        if region_choice == "__prev_pocket__":
            # 이전 포켓 결과 표시
            print(f"\n--- 이전 포켓 결과 ({len(prev_pockets)}개) ---")
            print(f"  {'#':>3}  {'폴더':<28} {'리간드':<18} {'포켓':>3}  {'Center':>24}  {'Energy':>8}  {'상태'}")
            print(f"  {'-'*3}  {'-'*28} {'-'*18} {'-'*3}  {'-'*24}  {'-'*8}  {'-'*4}")
            for i, p in enumerate(prev_pockets, 1):
                c = p["center"]
                energy_str = f"{p['best_energy']:.1f}" if p["best_energy"] is not None else "N/A"
                status = "완료" if p["completed"] else "미완"
                print(f"  [{i:>2}] {p['dir']:<28} {p['ligand']:<18} #{p['id']:<2}  "
                      f"[{c[0]:>7.1f}, {c[1]:>7.1f}, {c[2]:>7.1f}]  "
                      f"{energy_str:>8}  {status}")

            while True:
                raw = input("\n포켓 번호 선택: ").strip()
                try:
                    idx = int(raw) - 1
                    if 0 <= idx < len(prev_pockets):
                        break
                except ValueError:
                    pass
                print(f"  1~{len(prev_pockets)} 사이의 번호를 입력하세요.")

            selected = prev_pockets[idx]
            center = selected["center"]
            raw = input(f"Box size (X Y Z) [기본 30 30 30]: ").strip()
            if raw:
                try:
                    parts = raw.replace(",", " ").split()
                    box_size = [float(x) for x in parts]
                except ValueError:
                    box_size = [30.0, 30.0, 30.0]
            else:
                box_size = [30.0, 30.0, 30.0]
            print(f"  → 선택된 포켓: {selected['dir']} / {selected['ligand']} #{selected['id']}")
            print(f"  → Center: [{center[0]}, {center[1]}, {center[2]}]")
            print(f"  → Box size: {box_size}")

        elif region_choice == "__custom__":
            while True:
                raw = input("\nCenter (X Y Z): ").strip()
                try:
                    parts = raw.replace(",", " ").split()
                    center = [float(x) for x in parts]
                    if len(center) == 3:
                        break
                except ValueError:
                    pass
                print("  3개의 좌표를 공백으로 구분하여 입력하세요. (예: -38.2 0.2 -73.8)")

            raw = input("Box size (X Y Z) [기본 30 30 30]: ").strip()
            if raw:
                try:
                    parts = raw.replace(",", " ").split()
                    box_size = [float(x) for x in parts]
                except ValueError:
                    box_size = [30.0, 30.0, 30.0]
            else:
                box_size = [30.0, 30.0, 30.0]
        else:
            region = region_choice

    # Step 5: 포켓 탐색 설정 (선택)
    n_pockets = None
    max_per_pocket = 3
    cluster_radius = 5.0
    exclude_zones = []

    exclude_choice = ask_choice("포켓 탐색 모드", [
        ("기본           (단일 결과)", "none"),
        ("다중 포켓 탐색  (포켓을 하나씩 채우며 다양한 부위 발견)", "pockets"),
        ("수동 제외       (좌표 직접 입력)", "manual"),
    ], allow_quit=False)

    if exclude_choice == "pockets":
        n_pockets = ask_value("  몇 개 포켓까지 탐색?", default=5, cast=int)
        max_per_pocket = ask_value("  포켓당 포즈 수 (가득 차면 다음 포켓으로)", default=3, cast=int)
        cluster_radius = ask_value("  같은 포켓 판정 반경 (Å)", default=5.0, cast=float)
        print(f"  → {n_pockets}개 포켓 탐색, 포켓당 {max_per_pocket}개, 반경 {cluster_radius}Å")
    elif exclude_choice == "manual":
        while True:
            raw = input("\nExclusion center (X Y Z): ").strip()
            try:
                parts = raw.replace(",", " ").split()
                ec = [float(x) for x in parts]
                if len(ec) != 3:
                    raise ValueError
            except ValueError:
                print("  3개의 좌표를 입력하세요. (예: -38.2 0.2 -73.8)")
                continue

            radius = ask_value("  반경 (Å)", default=10.0, cast=float)
            exclude_zones.append(ec + [radius])
            print(f"  → Zone 추가: center=[{ec[0]}, {ec[1]}, {ec[2]}], radius={radius}Å")

            if not ask_confirm("  추가 Exclusion Zone 설정?", default=False):
                break

    # Step 6: 파라미터
    exhaustiveness = 128
    n_poses = 20
    energy_range = 5.0

    print(f"\n--- 파라미터 ---")
    print(f"  Exhaustiveness: {exhaustiveness}")
    print(f"  Poses: {n_poses}")
    print(f"  Energy range:   {energy_range}")

    if ask_confirm("변경하시겠습니까?", default=False):
        exhaustiveness = ask_value("  Exhaustiveness", default=128, cast=int)
        n_poses = ask_value("  Poses", default=20, cast=int)
        energy_range = ask_value("  Energy range", default=5.0, cast=float)

    # Step 7: 최종 확인
    lig_names = [l.name.replace("_ligand.pdbqt", "") for l in ligands]
    rec_names = [r.name.replace("_receptor.pdbqt", "") for r in receptors]
    print(f"\n{'='*60}")
    print(f"  설정 요약")
    print(f"{'='*60}")
    print(f"  Mode:           {mode}")
    print(f"  Receptor:       {', '.join(rec_names)} ({len(receptors)}개)")
    print(f"  Ligands:        {', '.join(lig_names)} ({len(ligands)}개)")
    if len(receptors) > 1:
        print(f"  조합:           {len(receptors)} × {len(ligands)} = {len(receptors)*len(ligands)}개 도킹")
    if mode == "focused":
        if region:
            preset = REGION_PRESETS[region]
            print(f"  Region:         {region}")
            c = preset['center']
            print(f"  Center:         [{c[0]}, {c[1]}, {c[2]}]")
            print(f"  Box:            {preset['box_size']}")
        else:
            print(f"  Center:         {center}")
            print(f"  Box:            {box_size}")
    if n_pockets:
        print(f"  Pocket 탐색:    {n_pockets}개 포켓 (포켓당 {max_per_pocket}개, 반경 {cluster_radius}Å)")
    if exclude_zones:
        print(f"  Exclude Zones:  {len(exclude_zones)}개")
        for ez in exclude_zones:
            print(f"    center=[{ez[0]}, {ez[1]}, {ez[2]}], radius={ez[3]}Å")
    print(f"  Exhaustiveness: {exhaustiveness}")
    print(f"  Poses:          {n_poses}")
    print(f"  Energy range:   {energy_range}")
    print(f"{'='*60}")

    if not ask_confirm("실행하시겠습니까?"):
        print("취소합니다.")
        sys.exit(0)

    return argparse.Namespace(
        mode=mode,
        receptor=None,  # multi-receptor: receptors 사용
        receptors=[str(r) for r in receptors],
        receptor_pdb=None,
        ligands=[str(l) for l in ligands],
        center=center,
        box_size=box_size,
        region=region,
        exclude_zone=exclude_zones if exclude_zones else None,
        n_pockets=n_pockets,
        max_per_pocket=max_per_pocket,
        cluster_radius=cluster_radius,
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
        energy_range=energy_range,
        padding=5.0,
        min_box=70.0,
        cpu=0,
        seed=None,
        parallel_receptors=len(receptors),
        output_dir=None,
        label=None,
        config=None,
        prepare_receptor=False,
        smiles=None,
    )


# ============================================================
# main
# ============================================================
def main():
    # 인자 없으면 인터랙티브 모드
    if len(sys.argv) == 1:
        args = interactive_mode()
    else:
        args = parse_args()

    # --config 로드
    if args.config:
        config = load_config(args.config)
        apply_config_to_args(args, config)
        if is_project_config(config):
            ensure_project_config_inputs_ready(config)

    # mode 필수 확인
    if args.mode is None:
        logger.error("--mode (blind/focused) 지정 필수")
        sys.exit(1)

    # 헤더 출력
    logger.info(
        "AutoDock Vina %s Docking Pipeline | exhaustiveness=%d | poses=%d | "
        "energy_range=%s | cpu=%s | seed=%s",
        args.mode.capitalize(), args.exhaustiveness, args.n_poses,
        args.energy_range, args.cpu, args.seed,
    )

    # --smiles CLI 입력 처리 (SMILES → SDF → PDBQT)
    if hasattr(args, "smiles") and args.smiles:
        logger.info("SMILES CLI에서 %d개 SMILES 입력 처리...", len(args.smiles))
        for i, smi in enumerate(args.smiles):
            name = f"cli_smiles_{i}"
            sdf_file = INPUT_DIR / f"{name}_ligand.sdf"
            pdbqt_file = INPUT_DIR / f"{name}_ligand.pdbqt"
            if not pdbqt_file.exists():
                if not sdf_file.exists():
                    smiles_to_sdf(smi, sdf_file, name)
                if sdf_file.exists():
                    prepare_ligand(sdf_file, pdbqt_file)

    # ── Receptor / Ligand 감지 ──
    input_dir = INPUT_DIR
    # Project config 모드: config에서 경로 직접 지정 → auto-detect 불필요
    if hasattr(args, "receptors") and args.receptors:
        all_receptors = [Path(r) for r in args.receptors]
        if args.ligands:
            args.ligands = [Path(l) for l in args.ligands]
        else:
            args.ligands = auto_detect_ligands(INPUT_DIR)
    else:
        # 비-project config: input/ 스캔
        ready_receptors, ready_ligands = preprocess_inputs(INPUT_DIR)
        input_dir = INPUT_DIR
        if args.receptor:
            all_receptors = [Path(args.receptor)]
        else:
            all_receptors = auto_detect_receptor(input_dir)
        if args.ligands:
            args.ligands = [Path(l) for l in args.ligands]
        else:
            args.ligands = auto_detect_ligands(input_dir)

    # 파일 존재 확인 (ligands)
    for lig in args.ligands:
        if not lig.exists():
            logger.error("Ligand 파일 없음: %s", lig)
            sys.exit(1)

    # Region 프리셋 적용
    if args.region:
        preset = REGION_PRESETS[args.region]
        args.center = preset["center"]
        args.box_size = preset["box_size"]
        logger.info("Region %s: %s", args.region, preset['description'])

    # 필터링 설정 (모든 receptor에 공통)
    exclude_zones = getattr(args, "exclude_zone", None) or []
    n_pockets_val = getattr(args, "n_pockets", None)
    max_per_pocket = getattr(args, "max_per_pocket", 3)
    c_radius = getattr(args, "cluster_radius", 5.0)

    if exclude_zones:
        logger.info("Exclusion 수동 제외 영역 %d개:", len(exclude_zones))
        for ez in exclude_zones:
            logger.info("  center=[%.1f, %.1f, %.1f], radius=%.0fÅ", ez[0], ez[1], ez[2], ez[3])
    if n_pockets_val:
        logger.info("Pocket Discovery %d개 포켓 탐색 (포켓당 %d개, 반경 %.1fÅ)",
                    n_pockets_val, max_per_pocket, c_radius)

    # 다중 포켓 탐색 시 포즈 수 자동 증가
    dock_n_poses = args.n_poses
    if n_pockets_val:
        min_needed = n_pockets_val * max_per_pocket * 5
        dock_n_poses = max(args.n_poses, min_needed)
        if dock_n_poses > args.n_poses:
            logger.info("Pocket Discovery 포즈 수 자동 증가: %d → %d", args.n_poses, dock_n_poses)

    if len(all_receptors) > 1:
        logger.info("다중 Receptor 도킹: %d개 × %d개 ligand", len(all_receptors), len(args.ligands))

    # ══════════════════════════════════════════════════════════
    # Receptor 루프: 각 receptor에 대해 전체 ligand 도킹
    # ══════════════════════════════════════════════════════════
    for ri, receptor_path in enumerate(all_receptors):
        if len(all_receptors) > 1:
            logger.info("Receptor [%d/%d]: %s", ri+1, len(all_receptors), receptor_path.name)

        args.receptor = receptor_path

        # 파일 존재 확인
        if not receptor_path.exists():
            logger.error("Receptor 파일 없음: %s", receptor_path)
            continue

        # Receptor PDB 감지 (per-receptor)
        receptor_pdb = None
        if args.mode == "blind":
            # 1) config의 receptor_pdb_map에서 찾기
            pdb_map = getattr(args, "receptor_pdb_map", {}) or {}
            mapped_pdb = pdb_map.get(str(receptor_path))
            if mapped_pdb and Path(mapped_pdb).exists():
                receptor_pdb = Path(mapped_pdb)
            elif args.receptor_pdb and len(all_receptors) == 1:
                receptor_pdb = Path(args.receptor_pdb)
            else:
                receptor_pdb = auto_detect_receptor_pdb(receptor_path, input_dir)
            if receptor_pdb is None:
                logger.error("Blind 모드에는 receptor PDB가 필요합니다. %s 건너뜁니다.", receptor_path.name)
                continue
            center, box_size = calc_box_from_pdb(receptor_pdb, args.padding, args.min_box)
            if "anp" in receptor_pdb.name.lower():
                check_anp(receptor_pdb)
        else:
            # focused mode
            if args.center is None:
                logger.error("Focused 모드에는 --center X Y Z 또는 --region 지정 필수")
                sys.exit(1)
            if args.box_size is None:
                args.box_size = [30.0, 30.0, 30.0]
            center = args.center
            box_size = args.box_size
            # PDB 감지 시도
            if args.receptor_pdb and len(all_receptors) == 1:
                receptor_pdb = Path(args.receptor_pdb)
            else:
                stem = receptor_path.name.replace("_receptor.pdbqt", "")
                pdb_cand = input_dir / f"{stem}.pdb"
                if pdb_cand.exists():
                    receptor_pdb = pdb_cand
            if receptor_pdb and "anp" in receptor_pdb.name.lower():
                check_anp(receptor_pdb)

        # Receptor 변환 (--prepare-receptor)
        if args.prepare_receptor and ri == 0:
            if receptor_pdb is None:
                receptor_pdb = auto_detect_receptor_pdb(receptor_path, input_dir)
            if receptor_pdb is None:
                logger.error("--prepare-receptor: PDB 파일을 찾을 수 없습니다.")
                continue
            pdbqt_out = input_dir / (receptor_pdb.stem + "_receptor.pdbqt")
            if not prepare_receptor(receptor_pdb, pdbqt_out):
                logger.error("Receptor 변환 실패: %s", receptor_pdb)
                continue
            args.receptor = pdbqt_out
            receptor_path = pdbqt_out

        # 출력 디렉토리 (receptor별)
        args.receptor_pdb = str(receptor_pdb) if receptor_pdb else None
        receptor_id = resolve_receptor_id(args, receptor_path)
        args.receptor_id = receptor_id
        receptor_output_dir = build_receptor_output_dir(args, receptor_id)
        if receptor_output_dir:
            args.output_dir = str(receptor_output_dir)
            args.label = receptor_id
        out_dir = setup_output_dir(args)

        # 도킹 실행
        all_results = {}
        for lig_path in args.ligands:
            lig_name = lig_path.name.replace("_ligand.pdbqt", "")
            output_file = out_dir / f"{lig_name}_{args.mode}.pdbqt"

            logger.info("%s %s docking 시작...", lig_name, args.mode)
            energies = run_docking(
                receptor_path, lig_path, output_file,
                center, box_size,
                args.exhaustiveness, dock_n_poses,
            )
            all_results[lig_name] = energies
            print_results(lig_name, args.mode, energies, center, box_size,
                          display_n=args.n_poses)

            # 다중 포켓 탐색
            if n_pockets_val:
                logger.info("Pocket Discovery %s 포켓 탐색 중...", lig_name)
                poses = parse_poses(output_file)
                pockets = discover_pockets(
                    poses, cluster_radius=c_radius,
                    max_per_pocket=max_per_pocket, n_pockets=n_pockets_val)

                completed_pockets = [p for p in pockets if p["closed"]]
                partial_pockets = [p for p in pockets if not p["closed"]]

                # 포켓 좌표 저장 (focused 재실행용)
                save_pockets_info(pockets, poses, energies, lig_name, out_dir)

                best_poses = []
                for pocket in completed_pockets:
                    best_idx = pocket["pose_indices"][0]
                    best_poses.append(poses[best_idx])

                if best_poses:
                    pockets_file = out_dir / f"{lig_name}_{args.mode}_pockets.pdbqt"
                    write_filtered_poses(best_poses, pockets_file)
                    logger.info("완성 포켓 best pose: %s (%d개)", pockets_file.name, len(best_poses))

                all_pocket_poses = []
                for pocket in completed_pockets:
                    for idx in pocket["pose_indices"]:
                        all_pocket_poses.append(poses[idx])
                if all_pocket_poses:
                    all_file = out_dir / f"{lig_name}_{args.mode}_all_pockets.pdbqt"
                    write_filtered_poses(all_pocket_poses, all_file)
                    logger.info("완성 포켓 전체 포즈: %s (%d개)", all_file.name, len(all_pocket_poses))

                if partial_pockets:
                    partial_poses = []
                    for pocket in partial_pockets:
                        best_idx = pocket["pose_indices"][0]
                        partial_poses.append(poses[best_idx])
                    partial_file = out_dir / f"{lig_name}_{args.mode}_partial_pockets.pdbqt"
                    write_filtered_poses(partial_poses, partial_file)
                    logger.info("미완성 포켓 (참고): %s (%d개, 각 %s개 포즈)",
                               partial_file.name, len(partial_poses),
                               [len(p['pose_indices']) for p in partial_pockets])

            elif exclude_zones:
                logger.info("Filter %s 포즈 필터링 중...", lig_name)
                poses = parse_poses(output_file)
                filtered = filter_poses_by_exclusion(poses, exclude_zones)
                if filtered:
                    filtered_file = out_dir / f"{lig_name}_{args.mode}_filtered.pdbqt"
                    write_filtered_poses(filtered, filtered_file)
                    logger.info("필터링 결과: %s (%d개 포즈)", filtered_file.name, len(filtered))

        # 요약
        print_summary(all_results, out_dir)

        # PyMOL 시각화 스크립트 생성
        style_pml = find_style_pml(INPUT_DIR, receptor_path)
        ligand_outputs = {name: out_dir / f"{name}_{args.mode}.pdbqt"
                          for name in all_results}
        r_pdb = receptor_pdb if receptor_pdb else receptor_path
        generate_pymol_script(
            out_dir, r_pdb, receptor_path,
            ligand_outputs, args.mode, style_pml,
        )


def run_taskgroup12_main():
    if len(sys.argv) == 1:
        args = interactive_mode()
    else:
        args = parse_args()

    if args.config:
        config = load_config(args.config)
        args = apply_config_to_args(args, config)
        if is_project_config(config):
            ensure_project_config_inputs_ready(config)

    if args.mode is None:
        logger.error("--mode (blind/focused) is required")
        sys.exit(1)

    logger.info(
        "AutoDock Vina %s Docking Pipeline | exhaustiveness=%d | poses=%d | "
        "energy_range=%s | cpu=%s | max_workers=%s",
        args.mode.capitalize(), args.exhaustiveness, args.n_poses,
        args.energy_range, args.cpu, args.max_workers,
    )

    if hasattr(args, "smiles") and args.smiles:
        logger.info("SMILES Processing %d CLI SMILES entries...", len(args.smiles))
        for i, smi in enumerate(args.smiles):
            name = f"cli_smiles_{i}"
            sdf_file = INPUT_DIR / f"{name}_ligand.sdf"
            pdbqt_file = INPUT_DIR / f"{name}_ligand.pdbqt"
            if not pdbqt_file.exists():
                if not sdf_file.exists():
                    smiles_to_sdf(smi, sdf_file, name)
                if sdf_file.exists():
                    prepare_ligand(sdf_file, pdbqt_file)

    preprocess_inputs(INPUT_DIR)
    input_dir = INPUT_DIR

    if hasattr(args, "receptors") and args.receptors:
        all_receptors = [Path(r) for r in args.receptors]
    elif args.receptor:
        all_receptors = [Path(args.receptor)]
    else:
        all_receptors = auto_detect_receptor(input_dir)

    if args.ligands:
        args.ligands = [Path(l) for l in args.ligands]
    else:
        args.ligands = auto_detect_ligands(input_dir)

    for lig in args.ligands:
        if not lig.exists():
            logger.error("Ligand file not found: %s", lig)
            sys.exit(1)

    if args.region:
        preset = REGION_PRESETS[args.region]
        args.center = preset["center"]
        args.box_size = preset["box_size"]
        logger.info("Region %s: %s", args.region, preset['description'])

    exclude_zones = getattr(args, "exclude_zone", None) or []
    n_pockets_val = getattr(args, "n_pockets", None)
    max_per_pocket = getattr(args, "max_per_pocket", 3)
    c_radius = getattr(args, "cluster_radius", 5.0)

    dock_n_poses = args.n_poses
    if n_pockets_val:
        min_needed = n_pockets_val * max_per_pocket * 5
        dock_n_poses = max(args.n_poses, min_needed)

    if len(all_receptors) > 1:
        logger.info("Multi-receptor docking: %d x %d ligands", len(all_receptors), len(args.ligands))

    for ri, receptor_path in enumerate(all_receptors):
        if len(all_receptors) > 1:
            logger.info("Receptor [%d/%d]: %s", ri+1, len(all_receptors), receptor_path.name)

        args.receptor = receptor_path
        receptor_id = resolve_receptor_id(args, receptor_path)
        args.receptor_id = receptor_id
        args.receptor_source_type = (getattr(args, "receptor_source_types", {}) or {}).get(str(receptor_path))

        if not receptor_path.exists():
            logger.error("Receptor file not found: %s", receptor_path)
            continue

        receptor_pdb = None
        receptor_pdb_override = resolve_receptor_pdb_override(args, receptor_path)
        if args.mode == "blind":
            if receptor_pdb_override is not None:
                receptor_pdb = receptor_pdb_override
            elif args.receptor_pdb and len(all_receptors) == 1:
                receptor_pdb = Path(args.receptor_pdb)
            else:
                receptor_pdb = auto_detect_receptor_pdb(receptor_path, input_dir)
            if receptor_pdb is None:
                logger.error("Blind mode requires receptor PDB: %s", receptor_path.name)
                continue
            center, box_size = calc_box_from_pdb(receptor_pdb, args.padding, args.min_box)
            if "anp" in receptor_pdb.name.lower():
                check_anp(receptor_pdb)
        else:
            if args.center is None:
                logger.error("Focused mode requires --center or --region")
                sys.exit(1)
            if args.box_size is None:
                args.box_size = [30.0, 30.0, 30.0]
            center = args.center
            box_size = args.box_size
            if receptor_pdb_override is not None:
                receptor_pdb = receptor_pdb_override
            elif args.receptor_pdb and len(all_receptors) == 1:
                receptor_pdb = Path(args.receptor_pdb)
            else:
                stem = receptor_path.name.replace("_receptor.pdbqt", "")
                pdb_cand = input_dir / f"{stem}.pdb"
                if pdb_cand.exists():
                    receptor_pdb = pdb_cand
            if receptor_pdb and "anp" in receptor_pdb.name.lower():
                check_anp(receptor_pdb)

        if args.prepare_receptor and ri == 0:
            if receptor_pdb is None:
                receptor_pdb = auto_detect_receptor_pdb(receptor_path, input_dir)
            if receptor_pdb is None:
                logger.error("--prepare-receptor requires receptor PDB")
                continue
            pdbqt_out = input_dir / (receptor_pdb.stem + "_receptor.pdbqt")
            if not prepare_receptor(receptor_pdb, pdbqt_out):
                logger.error("Receptor preparation failed: %s", receptor_pdb)
                continue
            args.receptor = pdbqt_out
            receptor_path = pdbqt_out

        args.receptor_pdb = str(receptor_pdb) if receptor_pdb else None
        receptor_output_dir = build_receptor_output_dir(args, receptor_id)
        args.output_dir = str(receptor_output_dir) if receptor_output_dir else None
        args.label = receptor_id
        out_dir = setup_output_dir(args)

        ligand_jobs: List[Dict[str, Any]] = []
        for lig_path in args.ligands:
            lig_name = lig_path.name.replace("_ligand.pdbqt", "")
            ligand_id = resolve_ligand_id(args, lig_path)
            output_file = out_dir / f"{lig_name}_{args.mode}.pdbqt"
            ligand_jobs.append({
                "receptor_id": receptor_id,
                "ligand_id": ligand_id,
                "ligand_name": lig_name,
                "receptor_path": receptor_path,
                "ligand_path": lig_path,
                "output_file": output_file,
                "center": center,
                "box_size": box_size,
                "exhaustiveness": args.exhaustiveness,
                "n_poses": dock_n_poses,
                "energy_range": args.energy_range,
                "cpu": args.cpu,
                "seed": derive_docking_seed(args.seed, receptor_id, ligand_id),
            })

        worker_count = resolve_ligand_worker_count(
            args.max_workers,
            len(ligand_jobs),
            cpu_per_job=args.cpu,
        )
        logger.info("Dispatch receptor=%s ligands=%d workers=%d", receptor_id, len(ligand_jobs), worker_count)

        if len(ligand_jobs) == 1 or worker_count <= 1:
            job_results = [execute_ligand_job(job) for job in ligand_jobs]
        else:
            job_results = []
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                future_map = {executor.submit(execute_ligand_job, job): job for job in ligand_jobs}
                for future in as_completed(future_map):
                    job_results.append(future.result())

        all_results = {}
        for result in job_results:
            append_job_status(out_dir, result)
            if result["status"] != "ok":
                logger.error(
                    "receptor=%s ligand=%s output=%s error=%s",
                    result['receptor_id'], result['ligand_id'],
                    result['output_file'], result['error'],
                )
                continue

            lig_name = result["ligand_name"]
            output_file = Path(result["output_file"])
            energies = result["energies"]
            logger.info("receptor=%s ligand=%s output=%s", result['receptor_id'], result['ligand_id'], output_file.name)
            all_results[lig_name] = energies
            print_results(lig_name, args.mode, energies, center, box_size, display_n=args.n_poses)

            if n_pockets_val:
                logger.info("Pocket Discovery %s pocket scan...", lig_name)
                poses = parse_poses(output_file)
                pockets = discover_pockets(
                    poses, cluster_radius=c_radius,
                    max_per_pocket=max_per_pocket, n_pockets=n_pockets_val)

                completed_pockets = [p for p in pockets if p["closed"]]
                partial_pockets = [p for p in pockets if not p["closed"]]
                save_pockets_info(pockets, poses, energies, lig_name, out_dir)

                best_poses = []
                for pocket in completed_pockets:
                    best_idx = pocket["pose_indices"][0]
                    best_poses.append(poses[best_idx])

                if best_poses:
                    pockets_file = out_dir / f"{lig_name}_{args.mode}_pockets.pdbqt"
                    write_filtered_poses(best_poses, pockets_file)
                    logger.info("Saved best pocket poses: %s (%d)", pockets_file.name, len(best_poses))

                all_pocket_poses = []
                for pocket in completed_pockets:
                    for idx in pocket["pose_indices"]:
                        all_pocket_poses.append(poses[idx])
                if all_pocket_poses:
                    all_file = out_dir / f"{lig_name}_{args.mode}_all_pockets.pdbqt"
                    write_filtered_poses(all_pocket_poses, all_file)
                    logger.info("Saved all completed pocket poses: %s (%d)", all_file.name, len(all_pocket_poses))

                if partial_pockets:
                    partial_poses = []
                    for pocket in partial_pockets:
                        best_idx = pocket["pose_indices"][0]
                        partial_poses.append(poses[best_idx])
                    partial_file = out_dir / f"{lig_name}_{args.mode}_partial_pockets.pdbqt"
                    write_filtered_poses(partial_poses, partial_file)
                    logger.info("Saved partial pocket poses: %s (%d)", partial_file.name, len(partial_poses))

            elif exclude_zones:
                logger.info("Filter %s exclusion filtering...", lig_name)
                poses = parse_poses(output_file)
                filtered = filter_poses_by_exclusion(poses, exclude_zones)
                if filtered:
                    filtered_file = out_dir / f"{lig_name}_{args.mode}_filtered.pdbqt"
                    write_filtered_poses(filtered, filtered_file)
                    logger.info("Saved filtered poses: %s (%d)", filtered_file.name, len(filtered))

        if all_results:
            print_summary(all_results, out_dir)
        else:
            logger.warning("No successful docking results for receptor=%s", receptor_id)

        style_pml = find_style_pml(INPUT_DIR, receptor_path)
        ligand_outputs = {name: out_dir / f"{name}_{args.mode}.pdbqt" for name in all_results}
        if ligand_outputs:
            r_pdb = receptor_pdb if receptor_pdb else receptor_path
            generate_pymol_script(
                out_dir, r_pdb, receptor_path,
                ligand_outputs, args.mode, style_pml,
            )

    if args.config and args.parse_results:
        from egfr_pipeline.vina.parse_poses import build_pose_table_from_config

        pose_table = build_pose_table_from_config(args.config)
        logger.info("Parse Wrote pose table: %s", pose_table)

        if args.extract_contacts:
            from egfr_pipeline.vina.pose_contacts import enrich_pose_table_with_contacts

            pose_table = enrich_pose_table_with_contacts(
                args.config,
                str(pose_table),
                args.contact_cutoff,
            )
            logger.info("Contacts Updated pose table: %s", pose_table)

        if args.cluster_pockets:
            from egfr_pipeline.vina.pocket_cluster import cluster_pose_table

            pose_table = cluster_pose_table(
                args.config,
                str(pose_table),
                cutoff=args.pocket_cutoff,
            )
            logger.info("Cluster Updated pose table with pocket ids: %s", pose_table)

        if args.summarize_pockets:
            from egfr_pipeline.vina.pocket_summary import summarize_from_config

            pocket_csv, drug_csv, occupancy_csv = summarize_from_config(args.config, str(pose_table))
            logger.info("Summary Wrote pocket table: %s", pocket_csv)
            logger.info("Summary Wrote drug pocket map: %s", drug_csv)
            logger.info("Summary Wrote residue occupancy: %s", occupancy_csv)

        if args.compare_pockets:
            from egfr_pipeline.vina.cross_receptor import compare_from_config

            comparison_csv = compare_from_config(
                args.config,
                centroid_cutoff=args.comparison_centroid_cutoff,
            )
            logger.info("Compare Wrote pocket comparison: %s", comparison_csv)

        if args.extract_ppi_residues:
            from egfr_pipeline.ppi.pyrosetta_extract import extract_pyrosetta_batch
            from egfr_pipeline.ppi.afm_extract import extract_afm_batch

            try:
                r, s = extract_pyrosetta_batch(args.config)
                logger.info("PPI PyRosetta residues: %s", r)
                logger.info("PPI PyRosetta summary: %s", s)
            except Exception as e:
                logger.warning("PPI PyRosetta extraction skipped: %s", e)
            try:
                r, s = extract_afm_batch(args.config)
                logger.info("PPI AFM residues: %s", r)
                logger.info("PPI AFM summary: %s", s)
            except Exception as e:
                logger.warning("PPI AFM extraction skipped: %s", e)

        if args.generate_report:
            from egfr_pipeline.report import generate_report

            report_path, combined_csv = generate_report(args.config)
            logger.info("Report Project report: %s", report_path)
            logger.info("Report Combined residue evidence: %s", combined_csv)

        if args.validate_outputs:
            from egfr_pipeline.validate import run_validation

            vresult = run_validation(args.config)
            logger.info(vresult.summary())


# ============================================================
# 단일 receptor 도킹 워커 (ProcessPoolExecutor 워커용)
# ============================================================

def dock_one_receptor(receptor_entry, ligand_entries, config):
    """단일 receptor에 대해 모든 ligand 도킹 (subprocess 워커용).

    Args:
        receptor_entry: config receptors 항목 (id, pdb, pdbqt 키 포함)
        ligand_entries: config ligands 리스트 (pdbqt 키 포함)
        config: 프로젝트 config dict (vina, output_root, project_name 등)

    Returns:
        (receptor_id, {ligand_name: energies}) 튜플
    """
    receptor_id = receptor_entry["id"]
    receptor_pdbqt = Path(receptor_entry["pdbqt"])
    receptor_pdb = Path(receptor_entry["pdb"])
    mode = config.get("mode", "blind")

    # box 계산
    padding = config.get("vina", {}).get("padding", 5.0)
    min_box = config.get("vina", {}).get("min_box", 70.0)
    center, box_size = calc_box_from_pdb(receptor_pdb, padding, min_box)

    # 출력 디렉토리
    out_dir = paths.wa_phase1_vina_receptor(config, receptor_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    resources_cfg, _warnings = resolve_resource_config(config)
    vina_cfg = config.get("vina", {}) if isinstance(config.get("vina"), dict) else {}
    exhaustiveness = vina_cfg.get("exhaustiveness", 8)
    n_poses = vina_cfg.get("n_poses", 5)
    energy_range = float(vina_cfg.get("energy_range", 3.0))
    cpu = int(resources_cfg["vina"].get("cpu_per_job", vina_cfg.get("cpu", 0)))
    base_seed = normalize_optional_int(vina_cfg.get("seed"))

    results = {}
    for lig_entry in ligand_entries:
        lig_path = Path(lig_entry["pdbqt"])
        lig_name = lig_path.name.replace("_ligand.pdbqt", "")
        output_file = out_dir / f"{lig_name}_{mode}.pdbqt"
        ligand_id = lig_entry.get("id", lig_name)
        job_seed = derive_docking_seed(base_seed, receptor_id, ligand_id)

        logger.info(
            "[%s] %s %s docking (cpu=%s, seed=%s, energy_range=%s) ...",
            receptor_id, lig_name, mode, cpu, job_seed, energy_range,
        )
        energies = run_docking(
            receptor_pdbqt, lig_path, output_file,
            center, box_size, exhaustiveness, n_poses,
            cpu=cpu,
            seed=job_seed,
            energy_range=energy_range,
        )
        results[lig_name] = energies
        print_results(lig_name, mode, energies, center, box_size, display_n=n_poses)

    # 요약
    lines = [f"[{receptor_id}] Summary - Best Affinity per Ligand"]
    for lig_name, energies in results.items():
        best = energies[0][0] if energies is not None and len(energies) > 0 else float("nan")
        lines.append(f"  {lig_name:>12}  {best:>15.2f} kcal/mol")
    lines.append(f"[done] {receptor_id} results: {out_dir}")
    logger.info("\n".join(lines))
    return receptor_id, results

