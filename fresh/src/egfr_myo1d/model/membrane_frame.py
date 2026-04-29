"""M1 Phase 5: Membrane frame generation per handoff §15.

Coordinate convention (handoff §15.1):
    Z+ = C2 symmetry axis = membrane normal (extracellular -> cytosolic)
    X+ = chain A -> chain B
    lower    = TM / JM / membrane-proximal kinase-domain side
    lateral  = outward-facing surface away from central dimer interface

Vectors are computed from input C-alpha coordinates. No numerical fallback
to literal axis vectors is allowed in computation paths.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.structure.pdb_parser import (
    PDBParseError,
    PDBStructure,
    parse_pdb,
)


# Residue range for membrane-normal computation per project convention
TM_JM_RESIDUE_START = 634
TM_JM_RESIDUE_END = 674


PRIMARY_STATES = ("EGFR_160-185", "EGFR_170-200")
REFERENCE_CONTROL_STATES = ("3GT8_raw",)
ALL_STATES = PRIMARY_STATES + REFERENCE_CONTROL_STATES


MEMBRANE_FRAME_QC_COLUMNS = [
    "state",
    "role",
    "frame_source",
    "status",
    "n_membrane_norm",
    "x_dimer_axis_norm",
    "centroid_distance",
    "n_tm_jm_residues",
    "warnings",
    "notes",
]


@dataclass
class StateMembraneFrame:
    state_id: str
    role: str
    frame_source: str
    n_membrane: tuple[float, float, float] | None = None
    x_dimer_axis: tuple[float, float, float] | None = None
    protomer_a_centroid: tuple[float, float, float] | None = None
    protomer_b_centroid: tuple[float, float, float] | None = None
    p_jm_anchor: tuple[float, float, float] | None = None
    n_membrane_norm: float | None = None
    x_dimer_axis_norm: float | None = None
    centroid_distance: float | None = None
    n_tm_jm_residues: int = 0
    status: str = "PASS"
    warnings: list[str] = field(default_factory=list)
    notes: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "frame_source": self.frame_source,
            "n_membrane": list(self.n_membrane) if self.n_membrane else None,
            "x_dimer_axis": list(self.x_dimer_axis) if self.x_dimer_axis else None,
            "protomer_a_centroid": list(self.protomer_a_centroid) if self.protomer_a_centroid else None,
            "protomer_b_centroid": list(self.protomer_b_centroid) if self.protomer_b_centroid else None,
            "p_jm_anchor": list(self.p_jm_anchor) if self.p_jm_anchor else None,
            "n_membrane_norm": self.n_membrane_norm,
            "x_dimer_axis_norm": self.x_dimer_axis_norm,
            "centroid_distance": self.centroid_distance,
            "n_tm_jm_residues": self.n_tm_jm_residues,
            "status": self.status,
            "warnings": self.warnings,
            "notes": self.notes,
        }

    def to_csv_row(self) -> dict[str, Any]:
        return {
            "state": self.state_id,
            "role": self.role,
            "frame_source": self.frame_source,
            "status": self.status,
            "n_membrane_norm": "" if self.n_membrane_norm is None else "{0:.6f}".format(self.n_membrane_norm),
            "x_dimer_axis_norm": "" if self.x_dimer_axis_norm is None else "{0:.6f}".format(self.x_dimer_axis_norm),
            "centroid_distance": "" if self.centroid_distance is None else "{0:.6f}".format(self.centroid_distance),
            "n_tm_jm_residues": self.n_tm_jm_residues,
            "warnings": ";".join(self.warnings),
            "notes": self.notes,
        }


# ---------------------------------------------------------------------------
# Geometry primitives (numpy-backed)
# ---------------------------------------------------------------------------


def _ca_atoms(structure: PDBStructure):
    return [a for a in structure.atoms if a.atom_name == "CA"]


def _ca_coords(atoms) -> np.ndarray:
    if not atoms:
        return np.zeros((0, 3), dtype=float)
    return np.asarray([(a.x, a.y, a.z) for a in atoms], dtype=float)


def _normalize_vector(vec: np.ndarray) -> tuple[tuple[float, float, float] | None, float]:
    """Return (unit_vector, post_normalization_norm).

    post_normalization_norm is 1.0 for any non-zero input (sanity check) and
    0.0 if the input has zero magnitude. Use ``centroid_distance`` and
    related raw-distance fields elsewhere in the report to capture
    pre-normalization magnitudes.
    """
    norm = float(np.linalg.norm(vec))
    if norm == 0.0:
        return None, 0.0
    unit = vec / norm
    unit_norm = float(np.linalg.norm(unit))
    return tuple(float(c) for c in unit), unit_norm


def _principal_axis(coords: np.ndarray) -> tuple[np.ndarray | None, str | None]:
    """Compute the largest-variance direction via SVD on centered coords."""
    if coords.shape[0] < 4:
        return None, "insufficient_atoms"
    centered = coords - coords.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return None, "svd_failed"
    return vt[0], None


# ---------------------------------------------------------------------------
# Per-state frame computation
# ---------------------------------------------------------------------------


def _resolve_role(state_id: str) -> str:
    if state_id in REFERENCE_CONTROL_STATES:
        return "crystallographic_reference_control"
    if state_id in PRIMARY_STATES:
        return "primary_membrane_validated_state"
    return "unknown_state"


def compute_membrane_frame(
    state_full_frame_pdb: Path | None,
    plus10_full_frame_pdb: Path | None,
    state_id: str,
    profile: str = "codex_dev",
) -> StateMembraneFrame:
    """Compute a state's membrane frame from coordinates.

    For ``state_id == "3GT8_raw"`` returns a reference-control marker without
    computing any vectors.

    For primary states the function tries (in order):
      1. ``state_full_frame_pdb`` if provided and present
      2. ``plus10_full_frame_pdb`` if provided and present (fallback)

    If neither is available the result has ``status="missing_frame_source"``
    in ``codex_dev`` (WARN) or FAIL in ``hpc_strict``; vectors are ``None``
    (no fabricated values).
    """
    role = _resolve_role(state_id)

    if state_id in REFERENCE_CONTROL_STATES:
        return StateMembraneFrame(
            state_id=state_id,
            role=role,
            frame_source="not_primary; optional alignment-derived only",
            status="reference_control",
            notes="3GT8_raw is not a primary membrane-frame source",
        )

    source_path: Path | None = None
    frame_source = "missing"
    if state_full_frame_pdb is not None and Path(state_full_frame_pdb).is_file():
        source_path = Path(state_full_frame_pdb)
        frame_source = "state_full_frame"
    elif plus10_full_frame_pdb is not None and Path(plus10_full_frame_pdb).is_file():
        source_path = Path(plus10_full_frame_pdb)
        frame_source = "plus10_inherited"

    if source_path is None:
        return StateMembraneFrame(
            state_id=state_id,
            role=role,
            frame_source="missing",
            status="FAIL" if profile == "hpc_strict" else "WARN",
            warnings=["missing_frame_source"],
            notes="no state full_frame and no plus10_full_frame available",
        )

    try:
        structure = parse_pdb(source_path)
    except PDBParseError as exc:
        return StateMembraneFrame(
            state_id=state_id,
            role=role,
            frame_source="missing",
            status="FAIL",
            warnings=["parse_error: {0}".format(exc)],
        )

    chains = set(structure.chain_ids)
    if not {"A", "B"}.issubset(chains):
        return StateMembraneFrame(
            state_id=state_id,
            role=role,
            frame_source=frame_source,
            status="FAIL" if profile == "hpc_strict" else "WARN",
            warnings=["chains_A_B_required_observed: {0}".format(sorted(chains))],
            notes="centroid computation requires explicit A and B chains",
        )

    cas = _ca_atoms(structure)
    cas_a = [a for a in cas if a.chain_id == "A"]
    cas_b = [a for a in cas if a.chain_id == "B"]
    if not cas_a or not cas_b:
        return StateMembraneFrame(
            state_id=state_id,
            role=role,
            frame_source=frame_source,
            status="FAIL",
            warnings=[
                "chain_CA_atoms_missing: A={0} B={1}".format(len(cas_a), len(cas_b))
            ],
        )

    coords_a = _ca_coords(cas_a)
    coords_b = _ca_coords(cas_b)
    centroid_a = coords_a.mean(axis=0)
    centroid_b = coords_b.mean(axis=0)
    centroid_distance = float(np.linalg.norm(centroid_b - centroid_a))

    x_unit, x_norm = _normalize_vector(centroid_b - centroid_a)

    tm_jm_cas = [
        a for a in cas
        if TM_JM_RESIDUE_START <= a.residue_number <= TM_JM_RESIDUE_END
    ]
    n_unit: tuple[float, float, float] | None = None
    n_norm: float | None = None
    p_jm_anchor: tuple[float, float, float] | None = None
    warnings: list[str] = []

    if tm_jm_cas:
        coords_tm = _ca_coords(tm_jm_cas)
        principal, err = _principal_axis(coords_tm)
        if principal is None:
            warnings.append("n_membrane_computation_{0}".format(err or "failed"))
        else:
            # Sign convention: prefer principal direction with positive Z component
            # (extracellular -> cytosolic axis). When Z is exactly zero we keep
            # the SVD orientation as-is.
            if principal[2] < 0:
                principal = -principal
            n_unit, n_norm = _normalize_vector(principal)
        anchor = coords_tm.mean(axis=0)
        p_jm_anchor = tuple(float(c) for c in anchor)
    else:
        warnings.append("missing_TM_JM_residues_634_674")

    status = "PASS"
    if n_unit is None or x_unit is None:
        status = "FAIL" if profile == "hpc_strict" else "WARN"
    if x_norm == 0.0 or (n_norm is not None and n_norm == 0.0):
        status = "FAIL"
        warnings.append("zero_norm_vector")

    if n_unit is not None and x_unit is not None:
        dot = float(np.dot(np.asarray(n_unit), np.asarray(x_unit)))
        if abs(dot) > 0.95:
            warnings.append("n_and_x_nearly_parallel:dot={0:.3f}".format(dot))
            if status == "PASS":
                status = "WARN"

    if frame_source == "plus10_inherited" and status == "PASS":
        status = "WARN"
        warnings.append("plus10_inherited_frame_used_as_fallback")

    return StateMembraneFrame(
        state_id=state_id,
        role=role,
        frame_source=frame_source,
        n_membrane=n_unit,
        x_dimer_axis=x_unit,
        protomer_a_centroid=tuple(float(c) for c in centroid_a),
        protomer_b_centroid=tuple(float(c) for c in centroid_b),
        p_jm_anchor=p_jm_anchor,
        n_membrane_norm=n_norm,
        x_dimer_axis_norm=x_norm,
        centroid_distance=centroid_distance,
        n_tm_jm_residues=len(tm_jm_cas),
        status=status,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


COORDINATE_CONVENTION = (
    "Z+ is membrane normal from extracellular to intracellular/cytosolic; "
    "X+ is chain A to chain B"
)

FRAME_SOURCE_POLICY = {
    "primary": "state full-frame coordinates when available",
    "fallback": "plus10_full_frame coordinate frame or TM/JM 634-674 axis",
    "reference_only": "3GT8_raw is not a primary membrane-frame source",
}


def write_state_aware_membrane_frame_json(
    ctx: RunContext,
    frames: Iterable[StateMembraneFrame],
) -> Path:
    payload: dict[str, Any] = {
        "coordinate_convention": COORDINATE_CONVENTION,
        "frame_source_policy": FRAME_SOURCE_POLICY,
        "states": {f.state_id: f.to_json() for f in frames},
        "timestamp": now_iso(),
    }
    target = ctx.manifest_dir / "membrane_frame.json"
    write_json(target, payload, ctx)
    return target


def write_membrane_frame_qc_csv(
    ctx: RunContext,
    frames: Iterable[StateMembraneFrame],
) -> Path:
    target = ctx.qc_dir / "membrane_frame_qc.csv"
    safe = ctx.require_within_run_dir(target)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MEMBRANE_FRAME_QC_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        for frame in frames:
            writer.writerow(frame.to_csv_row())
    return safe


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def _resolve_state_paths(
    ctx: RunContext,
    full_frame_source_override: Path | None,
) -> tuple[dict[str, Path | None], Path | None]:
    """Resolve per-state input paths from receptor_states.yaml plus an optional
    user-supplied ``--full-frame-source`` override. Returns (state -> path|None,
    plus10 fallback path|None)."""
    config = load_yaml_config(
        ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml"
    )
    input_files = config.get("input_files", {})
    state_paths: dict[str, Path | None] = {}
    for state_id in ALL_STATES:
        rel = input_files.get(state_id)
        candidate = (ctx.repo_root / rel) if rel else None
        if candidate is not None and candidate.is_file():
            state_paths[state_id] = candidate
        elif full_frame_source_override is not None:
            state_paths[state_id] = Path(full_frame_source_override)
        else:
            state_paths[state_id] = None
    plus10_rel = input_files.get("plus10_full_frame")
    plus10_path = (ctx.repo_root / plus10_rel) if plus10_rel else None
    if full_frame_source_override is not None:
        plus10_path = Path(full_frame_source_override)
    if plus10_path is not None and not plus10_path.is_file():
        plus10_path = None
    return state_paths, plus10_path


def run_membrane_frame_computation(
    ctx: RunContext,
    state_ids: Iterable[str] | None = None,
    full_frame_source: Path | None = None,
    profile: str = "codex_dev",
) -> tuple[list[StateMembraneFrame], str]:
    """Iterate states, compute frames, and emit consolidated outputs.

    Parameters
    ----------
    ctx
        run context; outputs go under ``ctx.run_dir``.
    state_ids
        states to compute. If ``None``, all known states are computed
        (``EGFR_160-185``, ``EGFR_170-200``, ``3GT8_raw``).
    full_frame_source
        optional override path used when a state's own input PDB is missing
        (acts as the plus10 fallback). Accepts ``Path`` or ``str``; if not a
        file, treated as missing.
    profile
        ``codex_dev`` or ``hpc_strict``.

    Returns
    -------
    (frames, overall_status)
        list of per-state results and the aggregate status
        (``PASS`` | ``PASS_WITH_WARNINGS`` | ``FAIL``).
    """
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))

    selected = list(state_ids) if state_ids else list(ALL_STATES)
    for state_id in selected:
        if state_id not in ALL_STATES:
            raise ValueError("unknown state_id: {0}".format(state_id))

    state_paths, plus10_path = _resolve_state_paths(ctx, full_frame_source)

    frames: list[StateMembraneFrame] = []
    for state_id in selected:
        frame = compute_membrane_frame(
            state_paths.get(state_id),
            plus10_path,
            state_id,
            profile=profile,
        )
        frames.append(frame)

    fail_count = sum(1 for f in frames if f.status == "FAIL")
    warn_count = sum(1 for f in frames if f.status == "WARN")
    if fail_count:
        overall = "FAIL"
    elif warn_count:
        overall = "PASS_WITH_WARNINGS"
    else:
        overall = "PASS"

    json_path = write_state_aware_membrane_frame_json(ctx, frames)
    csv_path = write_membrane_frame_qc_csv(ctx, frames)

    append_phase_status(
        ctx,
        phase="compute-membrane-frame",
        status="WARN" if overall == "PASS_WITH_WARNINGS" else overall,
        message=(
            "compute-membrane-frame overall={0} states={1} fail={2} warn={3}".format(
                overall, len(frames), fail_count, warn_count
            )
        ),
        details={
            "states": [f.state_id for f in frames],
            "fail_count": fail_count,
            "warn_count": warn_count,
            "manifest": str(json_path),
            "qc_csv": str(csv_path),
            "profile": profile,
            "score_bonus_allowed": False,
        },
    )

    return frames, overall
