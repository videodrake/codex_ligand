"""Tests for PDB-based orientation filter in lightdock_validation."""

import csv
import json
import textwrap
from pathlib import Path

import pytest

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

from egfr_pipeline.phase1.lightdock_validation import (
    _parse_pdb_atoms,
    check_lightdock_available,
    compute_orientation_score_from_pdb,
    extract_lightdock_interfaces,
    generate_lightdock_setup,
)

needs_numpy = pytest.mark.skipif(not _HAS_NUMPY, reason="numpy not available")


def _make_pdb_line(atom_name, resname, chain, resnum, x, y, z):
    """Format a single ATOM line in PDB format."""
    return (
        f"ATOM  {1:5d} {atom_name:<4s} {resname:3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
    )


def _build_oriented_pdb(tmp_path, orientation="pass"):
    """Build a synthetic PDB with controlled orientation geometry.

    orientation="pass": active face (sheets 8+9) CA->CB vectors point
        toward receptor → dot product > 0
    orientation="fail": active face points away from receptor
    """
    lines = []
    # Receptor chain A: cluster of CA atoms at origin
    for i, rn in enumerate(range(800, 810)):
        lines.append(_make_pdb_line("CA", "ALA", "A", rn, 0.0, 0.0, float(i)))

    # Partner chain B: active face residues (961-964, 968-972)
    # Place sheet at (20, 0, 0) — 20A from receptor
    active_resnums = [961, 962, 963, 964, 968, 969, 970, 971, 972]
    for i, rn in enumerate(active_resnums):
        y_offset = float(i) * 0.5
        lines.append(_make_pdb_line("CA", "VAL", "B", rn, 20.0, y_offset, 0.0))

        # CB placement controls orientation:
        if orientation == "pass":
            # CB points toward receptor (negative x direction from CA)
            lines.append(_make_pdb_line("CB", "VAL", "B", rn, 19.0, y_offset, 0.0))
        else:
            # CB points away from receptor (positive x direction from CA)
            lines.append(_make_pdb_line("CB", "VAL", "B", rn, 21.0, y_offset, 0.0))

    # Back face residues (977-997)
    for rn in [977, 978, 979, 980, 985, 986, 987, 988, 993, 994, 995, 996, 997]:
        lines.append(_make_pdb_line("CA", "ALA", "B", rn, 25.0, 0.0, 0.0))

    lines.append("END")
    pdb_path = tmp_path / f"test_{orientation}.pdb"
    pdb_path.write_text("\n".join(lines) + "\n")
    return pdb_path


def test_parse_pdb_atoms_extracts_ca_and_cb(tmp_path):
    pdb_content = "\n".join([
        _make_pdb_line("CA", "ALA", "A", 100, 1.0, 2.0, 3.0),
        _make_pdb_line("CB", "ALA", "A", 100, 1.5, 2.5, 3.5),
        _make_pdb_line("CA", "GLY", "B", 200, 4.0, 5.0, 6.0),
        _make_pdb_line("N", "ALA", "A", 100, 0.0, 0.0, 0.0),  # Should be skipped
        "END",
    ]) + "\n"
    pdb_path = tmp_path / "test.pdb"
    pdb_path.write_text(pdb_content)

    atoms = _parse_pdb_atoms(pdb_path, ("CA", "CB"))

    assert ("A", 100, "CA") in atoms
    assert ("A", 100, "CB") in atoms
    assert ("B", 200, "CA") in atoms
    assert ("A", 100, "N") not in atoms  # N atom not requested
    assert len(atoms) == 3

    x, y, z, resname = atoms[("A", 100, "CA")]
    assert abs(x - 1.0) < 0.01
    assert resname == "ALA"


@needs_numpy
def test_orientation_score_pass(tmp_path):
    pdb_path = _build_oriented_pdb(tmp_path, "pass")
    result = compute_orientation_score_from_pdb(pdb_path)

    assert result["orientation_class"] == "pass"
    assert result["orientation_score"] > 0
    assert result["n_active_face_ca"] == 9


@needs_numpy
def test_orientation_score_fail(tmp_path):
    pdb_path = _build_oriented_pdb(tmp_path, "fail")
    result = compute_orientation_score_from_pdb(pdb_path)

    assert result["orientation_class"] == "fail"
    assert result["orientation_score"] < 0
    assert result["n_active_face_ca"] == 9


@needs_numpy
def test_orientation_score_insufficient_data(tmp_path):
    """Only 1 active-face CA atom → insufficient_data."""
    lines = [
        _make_pdb_line("CA", "ALA", "A", 800, 0.0, 0.0, 0.0),
        _make_pdb_line("CA", "VAL", "B", 962, 10.0, 0.0, 0.0),
        "END",
    ]
    pdb_path = tmp_path / "test_insufficient.pdb"
    pdb_path.write_text("\n".join(lines) + "\n")

    result = compute_orientation_score_from_pdb(pdb_path)
    assert result["orientation_class"] == "insufficient_data"
    assert result["n_active_face_ca"] == 1


@needs_numpy
def test_orientation_score_no_receptor(tmp_path):
    """No receptor atoms → receptor_error."""
    lines = []
    for rn in [961, 962, 963, 964, 968, 969, 970, 971, 972]:
        lines.append(_make_pdb_line("CA", "VAL", "B", rn, 10.0, 0.0, 0.0))
        lines.append(_make_pdb_line("CB", "VAL", "B", rn, 11.0, 0.0, 0.0))
    lines.append("END")
    pdb_path = tmp_path / "test_no_receptor.pdb"
    pdb_path.write_text("\n".join(lines) + "\n")

    result = compute_orientation_score_from_pdb(pdb_path)
    assert result["orientation_class"] == "receptor_error"


def test_check_lightdock_available():
    available, msg = check_lightdock_available()
    # In this workspace LightDock is not installed
    assert isinstance(available, bool)
    assert isinstance(msg, str)


def test_generate_lightdock_setup_script_contains_per_swarm_loop(tmp_path):
    """Verify generated script contains per-swarm lgd_generate_conformations loop."""
    script_path = generate_lightdock_setup("3GT8_raw", tmp_path)
    script_content = script_path.read_text()

    # Must iterate per swarm
    assert "for swarm_dir in swarm_*/" in script_content
    assert "lgd_generate_conformations.py" in script_content

    # Must include lgd_rank.py
    assert "lgd_rank.py" in script_content

    # Must have availability check
    assert "command -v" in script_content

    # Must have completion marker
    assert ".lightdock_complete" in script_content

    # Must NOT have hardcoded CPU count (should use variable)
    assert "-c $N_CPUS" in script_content


def test_generate_lightdock_setup_absolute_workdir(tmp_path):
    """Verify generated script uses absolute WORKDIR path."""
    script_path = generate_lightdock_setup("3GT8_raw", tmp_path)
    script_content = script_path.read_text()

    # WORKDIR should be an absolute path
    for line in script_content.splitlines():
        if line.startswith("WORKDIR="):
            workdir_val = line.split("=", 1)[1].strip('"')
            assert workdir_val.startswith("/"), f"WORKDIR is relative: {workdir_val}"
            break


@needs_numpy
def test_extract_with_orientation_columns(tmp_path):
    """Verify extraction writes orientation columns to output CSVs."""
    state_name = "3GT8_raw"
    lightdock_dir = tmp_path / state_name / "lightdock"
    swarm_dir = lightdock_dir / "swarm_1"
    swarm_dir.mkdir(parents=True, exist_ok=True)

    # Create ranking file
    (lightdock_dir / "rank_by_scoring.list").write_text(
        "swarm_1/lightdock_1.pdb 12.0\n"
    )

    # Create a PDB with enough active-face residues for orientation scoring
    pdb_path = _build_oriented_pdb(swarm_dir, "pass")
    # Rename to match expected path
    pdb_path.rename(swarm_dir / "lightdock_1.pdb")

    # Create metadata
    with open(lightdock_dir / "lightdock_run_metadata.json", "w") as f:
        json.dump({"construct_type": "full_kinase_domain"}, f)

    iface_path, model_path = extract_lightdock_interfaces(
        state_name, tmp_path, max_poses=10
    )

    assert iface_path is not None
    assert model_path is not None

    # Verify orientation columns exist
    with open(iface_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert "orientation_validation_status" in reader.fieldnames
        assert "orientation_score" in reader.fieldnames
        # Should have pass orientation (built with "pass" geometry)
        orient_statuses = set(r["orientation_validation_status"] for r in rows)
        assert "pass" in orient_statuses

    with open(model_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert "orientation_score" in reader.fieldnames
        assert rows[0]["orientation_validation_status"] == "pass"
