"""Confidentiality and claim-safety scans for public report outputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from egfr_myo1d.compound.confidentiality import PrivateMapEntry, scan_internal_id_leaks
from egfr_myo1d.core.run_context import RunContext


FORBIDDEN_CLAIM_PATTERNS = [
    r"\bvalidated inhibitor\b",
    r"\bproven binder\b",
    r"\bproven PPI inhibitor\b",
    r"\bclinically relevant drug candidate\b",
    r"\bvalidated EGFR inhibitor\b",
    r"\bvalidated EGFR-MYO1D PPI inhibitor\b",
    r"\bproven EGFR binder\b",
    r"\bclinically validated\b",
    r"\brecommended final candidate\b",
    r"\bfinal candidate recommendation\b",
    r"\bfinal candidate selected\b",
]

PRIVATE_STRUCTURE_PATTERN = re.compile(
    r"\b(?:SMILES|canonical_smiles|isomeric_smiles)\b\s*[:=,]\s*[A-Za-z0-9@+\-\[\]\(\)=#$\\/]{3,}",
    re.IGNORECASE,
)


@dataclass
class SafetyScanResult:
    leaks_detected: int
    smiles_logged: bool
    ligand_coordinates_logged: bool
    receptor_coordinates_logged: bool
    pose_coordinates_logged_outside_pose_file: bool
    candidate_overclaims_detected: bool
    coordinate_leak_files: list[str]
    public_outputs_scanned: list[str]

    @property
    def passed(self) -> bool:
        return (
            self.leaks_detected == 0
            and not self.smiles_logged
            and not self.ligand_coordinates_logged
            and not self.receptor_coordinates_logged
            and not self.pose_coordinates_logged_outside_pose_file
            and not self.candidate_overclaims_detected
        )


def public_report_scan_paths(ctx: RunContext) -> list[Path]:
    roots = [
        ctx.run_dir / "phase3_compounds" / "manifests",
        ctx.run_dir / "phase3_compounds" / "qc",
        ctx.run_dir / "phase3_compounds" / "reports",
        ctx.run_dir / "phase3_compounds" / "qsub",
        ctx.run_dir / "phase3_compounds" / "tables",
        ctx.run_dir / "report",
        ctx.logs_dir,
    ]
    files: list[Path] = []
    for root in roots:
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if path.is_file():
                files.append(ctx.require_within_run_dir(path))
    return files


def scan_public_reports(ctx: RunContext, private_entries: Iterable[PrivateMapEntry]) -> SafetyScanResult:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    scanned: list[str] = []
    coord_hits: list[str] = []
    smiles_logged = False
    ligand_coords = False
    receptor_coords = False
    pose_coords = False
    overclaim = False
    claim_regex = re.compile("|".join(FORBIDDEN_CLAIM_PATTERNS), re.IGNORECASE)
    for path in public_report_scan_paths(ctx):
        if "docking_outputs" in path.parts and path.suffix.lower() == ".pdbqt":
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        scanned.append(ctx.relative_to_repo(path))
        if PRIVATE_STRUCTURE_PATTERN.search(text):
            smiles_logged = True
        if re.search(r"^(ATOM|HETATM)\s+\d+", text, re.MULTILINE):
            rel = ctx.relative_to_repo(path)
            coord_hits.append(rel)
            lower = rel.lower()
            ligand_coords = ligand_coords or "ligand" in lower
            receptor_coords = receptor_coords or "receptor" in lower
            pose_coords = pose_coords or any(token in lower for token in ["pose", "candidate", "evidence", "report", "qc", "table"])
        if claim_regex.search(text):
            overclaim = True
    return SafetyScanResult(
        leaks_detected=len(leaks),
        smiles_logged=smiles_logged,
        ligand_coordinates_logged=ligand_coords,
        receptor_coordinates_logged=receptor_coords,
        pose_coordinates_logged_outside_pose_file=pose_coords,
        candidate_overclaims_detected=overclaim,
        coordinate_leak_files=coord_hits,
        public_outputs_scanned=sorted(set(scanned)),
    )
