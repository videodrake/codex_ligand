"""Group 0 setup verification tests.

Checks:
1. All 12 modification-target files exist.
2. region_definitions.py constants load correctly.
3. get_region() returns expected values.
4. Four regions are mutually exclusive and cover the full kinase domain.
5. Ko et al. sheet constants are consistent with orientation_filter.py.
"""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ── 1. Modification-target files exist ────────────────────────────────────

_TARGET_FILES = [
    "egfr_pipeline/verdict.py",
    "egfr_pipeline/report.py",
    "egfr_pipeline/validate.py",
    "egfr_pipeline/paths.py",
    "egfr_pipeline/vina/pocket_summary.py",
    "egfr_pipeline/vina/pocket_stability.py",
    "egfr_pipeline/phase1/orientation_filter.py",
    "egfr_pipeline/phase1/cluster_consensus.py",
    "egfr_pipeline/phase1/compare_states.py",
    "egfr_pipeline/phase1/review_report.py",
    "egfr_pipeline/phase1/lightdock_validation.py",
    "egfr_pipeline/pyrosetta_docking/pipeline_manager.py",
]


@pytest.mark.parametrize("relpath", _TARGET_FILES)
def test_target_file_exists(relpath):
    assert (PROJECT_ROOT / relpath).is_file(), f"Missing: {relpath}"


# ── 2. region_definitions constants ───────────────────────────────────────

def test_import_region_definitions():
    from egfr_pipeline.region_definitions import (
        get_region,
        ATP_SITE_RESIDUES,
        N_LOBE,
        ATP_SITE,
        C_LOBE_SURFACE,
        C_LOBE_CORE,
        REGION_LOOKUP,
        KINASE_DOMAIN_START,
        KINASE_DOMAIN_END,
    )
    assert KINASE_DOMAIN_START == 699
    assert KINASE_DOMAIN_END == 1007


def test_atp_site_residues_count():
    from egfr_pipeline.region_definitions import ATP_SITE_RESIDUES

    assert len(ATP_SITE_RESIDUES) == 37


def test_atp_site_contains_core_residues():
    from egfr_pipeline.region_definitions import ATP_SITE_RESIDUES

    core_five = {745, 790, 791, 793, 855}
    assert core_five.issubset(ATP_SITE_RESIDUES)


# ── 3. get_region() ──────────────────────────────────────────────────────

def test_get_region_atp_site():
    from egfr_pipeline.region_definitions import get_region

    assert get_region(745) == "atp_site"
    assert get_region(790) == "atp_site"
    assert get_region(855) == "atp_site"


def test_get_region_n_lobe():
    from egfr_pipeline.region_definitions import get_region

    assert get_region(700) == "n_lobe"
    assert get_region(750) == "n_lobe"


def test_get_region_c_lobe_surface():
    from egfr_pipeline.region_definitions import get_region

    # 860 is in C-lobe, not ATP site, not core
    result = get_region(860)
    assert result in ("c_lobe_surface", "c_lobe_core")


def test_get_region_outside_domain():
    from egfr_pipeline.region_definitions import get_region

    assert get_region(500) is None
    assert get_region(1100) is None


# ── 4. Mutual exclusivity and completeness ────────────────────────────────

def test_regions_mutually_exclusive():
    from egfr_pipeline.region_definitions import (
        N_LOBE, ATP_SITE, C_LOBE_SURFACE, C_LOBE_CORE,
    )

    assert not (N_LOBE & ATP_SITE)
    assert not (N_LOBE & C_LOBE_SURFACE)
    assert not (N_LOBE & C_LOBE_CORE)
    assert not (ATP_SITE & C_LOBE_SURFACE)
    assert not (ATP_SITE & C_LOBE_CORE)
    assert not (C_LOBE_SURFACE & C_LOBE_CORE)


def test_regions_cover_full_domain():
    from egfr_pipeline.region_definitions import (
        N_LOBE, ATP_SITE, C_LOBE_SURFACE, C_LOBE_CORE,
        KINASE_DOMAIN_START, KINASE_DOMAIN_END,
    )

    all_domain = frozenset(range(KINASE_DOMAIN_START, KINASE_DOMAIN_END + 1))
    all_regions = N_LOBE | ATP_SITE | C_LOBE_SURFACE | C_LOBE_CORE
    assert all_regions == all_domain
    assert len(all_regions) == 309


# ── 5. Ko et al. sheets ──────────────────────────────────────────────────

def test_ko_sheet_ranges():
    from egfr_pipeline.region_definitions import (
        KO_SHEET_8, KO_SHEET_9, KO_SHEET_10, KO_SHEET_11, KO_SHEET_12,
    )

    assert KO_SHEET_8 == frozenset({961, 962, 963, 964})
    assert KO_SHEET_9 == frozenset({968, 969, 970, 971, 972})
    assert KO_SHEET_10 == frozenset({977, 978, 979, 980})
    assert KO_SHEET_11 == frozenset({985, 986, 987, 988})
    assert KO_SHEET_12 == frozenset({993, 994, 995, 996, 997})


def test_ko_active_face():
    from egfr_pipeline.region_definitions import KO_SHEET_8_9

    assert len(KO_SHEET_8_9) == 9
    assert 962 in KO_SHEET_8_9
    assert 971 in KO_SHEET_8_9


# ── 6. SDF ligand files exist ────────────────────────────────────────────

_LIGAND_FILES = [
    "input/ligands/173940_ligand.sdf",
    "input/ligands/97806_ligand.sdf",
    "input/ligands/VAX-C12_0_ligand.sdf",
]


@pytest.mark.parametrize("relpath", _LIGAND_FILES)
def test_ligand_sdf_exists(relpath):
    assert (PROJECT_ROOT / relpath).is_file(), f"Missing: {relpath}"
