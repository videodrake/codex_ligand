"""Group G tests: Dimer input recognition and preprocessing.

Tests for:
- detect_dimer() function
- extract_dimer_chains() function
- RECEPTOR_STATES metadata (is_dimer, dimer_chain)
- prepare_docking_pair() mapping CSV generation
- excluded_residues dimer coverage
- PPI_TARGETS dimer configuration
"""
from __future__ import annotations

import csv
from pathlib import Path

from egfr_pipeline.phase1 import prepare_inputs as module


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resnum: int,
) -> str:
    element = atom_name.strip().lstrip("0123456789")[0]
    x = float(serial)
    y = float(serial) + 0.5
    z = float(serial) + 1.0
    return (
        f"ATOM  {serial:5d} {atom_name:>4} {resname:>3} {chain}{resnum:4d}"
        f"    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}\n"
    )


def _write_pdb(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "TER\nEND\n", encoding="utf-8")


def _make_monomer_lines(chain: str, start: int = 699, count: int = 5) -> list[str]:
    """Create minimal monomer PDB lines (N, CA, C, O per residue)."""
    lines = []
    serial = 1
    for i in range(count):
        rn = start + i
        for atom in [" N  ", " CA ", " C  ", " O  "]:
            lines.append(_atom_line(serial, atom, "ALA", chain, rn))
            serial += 1
    return lines


def _make_dimer_single_chain(chain: str, start: int = 699, count: int = 5) -> list[str]:
    """Create dimer in single chain with duplicate residue numbers."""
    mono1 = _make_monomer_lines(chain, start, count)
    mono2 = _make_monomer_lines(chain, start, count)
    return mono1 + mono2


def _make_dimer_two_chains(
    chain_a: str, chain_b: str, start: int = 699, count: int = 5,
) -> list[str]:
    """Create dimer with two separate chains."""
    return (
        _make_monomer_lines(chain_a, start, count)
        + _make_monomer_lines(chain_b, start, count)
    )


# ---------------------------------------------------------------------------
# Tests: detect_dimer
# ---------------------------------------------------------------------------

def test_detect_dimer_monomer(tmp_path: Path) -> None:
    """Single monomer should NOT be detected as dimer."""
    pdb = tmp_path / "mono.pdb"
    _write_pdb(pdb, _make_monomer_lines("A", 699, 5))
    result = module.detect_dimer(pdb, "A")
    assert result["is_dimer"] is False
    assert result["n_duplicated_atoms"] == 0


def test_detect_dimer_single_chain_dimer(tmp_path: Path) -> None:
    """Dimer in single chain (duplicate residues) should be detected."""
    pdb = tmp_path / "dimer_x.pdb"
    # 100 residues × 4 atoms × 2 copies = duplicate_atoms >> 50
    lines = _make_dimer_single_chain("X", 699, 100)
    _write_pdb(pdb, lines)
    result = module.detect_dimer(pdb, "X")
    assert result["is_dimer"] is True
    assert result["n_duplicated_atoms"] > 50
    assert result["detection_method"] == "duplicate_atoms"


def test_detect_dimer_two_chain(tmp_path: Path) -> None:
    """Dimer with two separate chains should be detected."""
    pdb = tmp_path / "dimer_ab.pdb"
    lines = _make_dimer_two_chains("A", "B", 699, 309)
    _write_pdb(pdb, lines)
    result = module.detect_dimer(pdb, "A")
    assert result["is_dimer"] is True
    assert result["detection_method"] == "multi_chain"
    assert "A" in result["chains"]
    assert "B" in result["chains"]


# ---------------------------------------------------------------------------
# Tests: extract_dimer_chains
# ---------------------------------------------------------------------------

def test_extract_dimer_two_chains(tmp_path: Path) -> None:
    """Two-chain dimer extracts both monomers correctly."""
    pdb = tmp_path / "dimer_ab.pdb"
    _write_pdb(pdb, _make_dimer_two_chains("A", "B", 699, 10))
    state_info = {"source_chain": "A", "dimer_chain": "B"}
    a_lines, a_res, b_lines, b_res = module.extract_dimer_chains(
        pdb, state_info, 699, 708
    )
    assert len(a_res) == 10
    assert len(b_res) == 10
    # Chains should be renamed
    assert all(module.get_chain(l) == "A" for l in a_lines)
    assert all(module.get_chain(l) == "B" for l in b_lines)


def test_extract_dimer_single_chain(tmp_path: Path) -> None:
    """Single-chain dimer splits into two monomers by occurrence."""
    pdb = tmp_path / "dimer_x.pdb"
    _write_pdb(pdb, _make_dimer_single_chain("X", 699, 10))
    state_info = {"source_chain": "X", "dimer_chain": "X"}
    a_lines, a_res, b_lines, b_res = module.extract_dimer_chains(
        pdb, state_info, 699, 708
    )
    assert len(a_res) == 10
    assert len(b_res) == 10
    # First occurrence → A, second → B
    assert all(module.get_chain(l) == "A" for l in a_lines)
    assert all(module.get_chain(l) == "B" for l in b_lines)


# ---------------------------------------------------------------------------
# Tests: RECEPTOR_STATES metadata
# ---------------------------------------------------------------------------

def test_receptor_states_all_dimer() -> None:
    """All RECEPTOR_STATES should have is_dimer=True."""
    for name, info in module.RECEPTOR_STATES.items():
        assert info["is_dimer"] is True, f"{name} should be dimer"
        assert "dimer_chain" in info, f"{name} missing dimer_chain"


def test_receptor_states_descriptions_no_monomer() -> None:
    """Descriptions should not claim structure is a monomer."""
    for name, info in module.RECEPTOR_STATES.items():
        desc = info["description"].lower()
        assert "monomer" not in desc, f"{name}: desc still says monomer"
        assert "dimer" in desc, f"{name}: desc should mention dimer"


# ---------------------------------------------------------------------------
# Tests: +1000 offset constants
# ---------------------------------------------------------------------------

def test_chain_b_offset() -> None:
    assert module.CHAIN_B_OFFSET == 1000


def test_membrane_proximal_b_defined() -> None:
    assert module.MEMBRANE_PROXIMAL_B
    nums = module.parse_ranges(module.MEMBRANE_PROXIMAL_B)
    assert len(nums) > 0
    # All should be in kinase domain range
    assert all(699 <= n <= 1007 for n in nums)


# ---------------------------------------------------------------------------
# Tests: excluded_residues in generate_configs
# ---------------------------------------------------------------------------

def test_generate_configs_excluded_includes_dimer() -> None:
    """EXCLUDED_RESIDUES_A should include +1000 offset residues."""
    from egfr_pipeline.phase1.generate_configs import EXCLUDED_RESIDUES_A
    nums = module.parse_ranges(EXCLUDED_RESIDUES_A)
    # Should have both monomer A and monomer B (+1000) residues
    monomer_a = [n for n in nums if n <= 1007]
    monomer_b = [n for n in nums if n > 1007]
    assert len(monomer_a) > 0, "Missing monomer A excluded residues"
    assert len(monomer_b) > 0, "Missing monomer B excluded residues (offset)"
    # Monomer B should be in 1700-2000 range
    assert all(1700 <= n <= 2000 for n in monomer_b)


# ---------------------------------------------------------------------------
# Tests: PPI_TARGETS configuration
# ---------------------------------------------------------------------------

def test_ppi_targets_dimer_offset() -> None:
    """All PPI_TARGETS should have construct_type=dimer_offset."""
    from run_production import PPI_TARGETS
    for t in PPI_TARGETS:
        assert t["construct_type"] == "dimer_offset", f"{t['name']}: wrong construct_type"


def test_ppi_targets_mapping_csv() -> None:
    """All PPI_TARGETS should have non-empty mapping_csv."""
    from run_production import PPI_TARGETS
    for t in PPI_TARGETS:
        assert t["mapping_csv"], f"{t['name']}: missing mapping_csv"
        assert t["mapping_csv"].endswith("_mapping.csv")


def test_ppi_targets_count() -> None:
    """Should be 3 states × 5 seeds = 15 targets."""
    from run_production import PPI_TARGETS
    assert len(PPI_TARGETS) == 15


# ---------------------------------------------------------------------------
# Tests: patch_ingestion accepts dimer_offset
# ---------------------------------------------------------------------------

def test_valid_construct_types_includes_dimer() -> None:
    from egfr_pipeline.phase2.patch_ingestion import VALID_CONSTRUCT_TYPES
    assert "dimer_offset" in VALID_CONSTRUCT_TYPES
    assert "full_kinase_domain" in VALID_CONSTRUCT_TYPES  # backwards compat


# ---------------------------------------------------------------------------
# Tests: extract_chain_residue_range (monomer extraction)
# ---------------------------------------------------------------------------

def test_extract_chain_first_occurrence_only(tmp_path: Path) -> None:
    """extract_chain_residue_range should pick first occurrence only."""
    pdb = tmp_path / "dimer_x.pdb"
    _write_pdb(pdb, _make_dimer_single_chain("X", 699, 10))
    lines, resnums = module.extract_chain_residue_range(
        pdb, "X", 699, 708, target_chain="A"
    )
    assert len(resnums) == 10
    # Should have 4 atoms per residue (not 8 = 4 × 2 copies)
    ca_count = sum(1 for l in lines if " CA " in l)
    assert ca_count == 10, f"Expected 10 CA atoms (monomer), got {ca_count}"
