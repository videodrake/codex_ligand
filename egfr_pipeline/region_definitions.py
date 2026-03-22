"""Central residue classification module for EGFR kinase domain.

Provides the canonical definitions for:
- EGFR 5-region classification (n_lobe, atp_site, c_lobe_surface, c_lobe_core)
- ATP binding site residues (37 residues)
- Ko et al. MYO1D beta-meander sheet definitions

**Provenance:** ATP_SITE_RESIDUES and the 4-region classification are derived
from published EGFR kinase crystal structures (PDB 2ITY, 1M17, 3GT8) and
standard kinase domain secondary structure assignments.  The exact surface/core
split is an approximation — B-1.1 (FreeSASA 3-state union, SASA > 10.0 A^2)
has not yet been performed.  Update this file when FreeSASA results are
available.

**Numbering:** All EGFR residue numbers use PDB numbering (3GT8-consistent).
MYO1D residues also use their source PDB numbering.
"""

from typing import Optional

# ============================================================
# Kinase domain boundaries
# ============================================================

KINASE_DOMAIN_START = 699
KINASE_DOMAIN_END = 1007  # inclusive; 699-1007 = 309 residues
NLOBE_CLOBE_BOUNDARY = 838  # residues < 838 are N-lobe, >= 838 are C-lobe

# ============================================================
# ATP binding site (37 residues)
# ============================================================
# Literature: PDB 2ITY (EGFR-erlotinib), 1M17 (EGFR-lapatinib),
# 3GT8 (EGFR kinase domain).  Includes residues within ~4 A of
# ATP/inhibitor plus key catalytic residues.

ATP_SITE_RESIDUES: frozenset = frozenset({
    # P-loop / glycine-rich loop (cap + loop)
    718, 719, 720, 721, 722, 723,
    # beta-3 strand (pocket floor)
    726,
    # alphaC-helix (pocket wall, catalytic lysine)
    743, 744, 745,
    # alphaC-helix distal (E762, K745 salt bridge partner)
    762,
    # beta-4 strand (back wall)
    766,
    # beta-4/beta-5 linker
    777,
    # Hinge region
    788, 789, 790, 791, 792, 793, 794, 795, 796,
    # Post-hinge (C797 = osimertinib covalent target)
    797,
    # alphaD-helix (selectivity pocket entrance)
    800,
    # Catalytic loop (includes HRD motif: H835, R836, D837)
    831, 832, 833, 834, 835, 836, 837,
    # Pre-DFG
    844,
    # Activation loop / DFG motif (D855, F856, G857)
    854, 855, 856, 857, 858,
})

assert len(ATP_SITE_RESIDUES) == 37, (
    f"ATP_SITE_RESIDUES should have 37 members, got {len(ATP_SITE_RESIDUES)}"
)

# ============================================================
# 4-region classification (309 residues total)
# ============================================================
# Approximation without SASA data.  Strategy:
#   1. ATP site residues assigned first (takes priority).
#   2. Remaining N-lobe residues (699-837 minus ATP) → n_lobe.
#   3. C-lobe core: deeply buried hydrophobic residues (838-1007)
#      identified from secondary structure (interior of helices,
#      beta-sheet core, far from solvent).
#   4. C-lobe surface: all remaining C-lobe residues.

# C-lobe core: residues buried in the hydrophobic interior.
# Estimated from helical/strand interiors in EGFR C-lobe (PDB 3GT8).
_C_LOBE_CORE_RESIDUES = frozenset({
    # alphaE-helix interior
    845, 846, 848, 850, 852,
    # alphaF-helix core (deeply buried)
    862, 866, 869, 870,
    # alphaG-helix core
    878, 882, 885, 886,
    # alphaH-helix core
    896, 899, 900, 903,
    # alphaI-helix core
    924, 927, 928, 931,
    # beta-7/8 strand core
    938, 941, 944, 947,
    # C-terminal buried
    955, 960, 965, 975,
})

# Build region sets
_all_kinase = frozenset(range(KINASE_DOMAIN_START, KINASE_DOMAIN_END + 1))
assert len(_all_kinase) == 309

_n_lobe_range = frozenset(r for r in _all_kinase if r < NLOBE_CLOBE_BOUNDARY)
_c_lobe_range = frozenset(r for r in _all_kinase if r >= NLOBE_CLOBE_BOUNDARY)

N_LOBE: frozenset = _n_lobe_range - ATP_SITE_RESIDUES
ATP_SITE: frozenset = ATP_SITE_RESIDUES  # alias for consistency
C_LOBE_CORE: frozenset = _C_LOBE_CORE_RESIDUES - ATP_SITE_RESIDUES
C_LOBE_SURFACE: frozenset = _c_lobe_range - ATP_SITE_RESIDUES - C_LOBE_CORE

# Validate mutual exclusivity and completeness
assert N_LOBE | ATP_SITE | C_LOBE_SURFACE | C_LOBE_CORE == _all_kinase
assert not (N_LOBE & ATP_SITE)
assert not (N_LOBE & C_LOBE_SURFACE)
assert not (N_LOBE & C_LOBE_CORE)
assert not (ATP_SITE & C_LOBE_SURFACE)
assert not (ATP_SITE & C_LOBE_CORE)
assert not (C_LOBE_SURFACE & C_LOBE_CORE)

# ============================================================
# Region lookup
# ============================================================

REGION_LOOKUP: dict = {}
for _r in N_LOBE:
    REGION_LOOKUP[_r] = "n_lobe"
for _r in ATP_SITE:
    REGION_LOOKUP[_r] = "atp_site"
for _r in C_LOBE_SURFACE:
    REGION_LOOKUP[_r] = "c_lobe_surface"
for _r in C_LOBE_CORE:
    REGION_LOOKUP[_r] = "c_lobe_core"


def get_region(residue_number: int) -> Optional[str]:
    """Return the region name for an EGFR residue number.

    Returns one of ``"n_lobe"``, ``"atp_site"``, ``"c_lobe_surface"``,
    ``"c_lobe_core"``, or ``None`` if the residue is outside the kinase
    domain (699-1007).
    """
    return REGION_LOOKUP.get(residue_number)


# ============================================================
# Ko et al. MYO1D beta-meander sheet definitions
# ============================================================
# Values match orientation_filter.py and prepare_inputs.py.
# These are beta-strand core residues (MYO1D PDB numbering).

KO_SHEET_8: frozenset = frozenset({961, 962, 963, 964})
KO_SHEET_9: frozenset = frozenset({968, 969, 970, 971, 972})
KO_SHEET_10: frozenset = frozenset({977, 978, 979, 980})
KO_SHEET_11: frozenset = frozenset({985, 986, 987, 988})
KO_SHEET_12: frozenset = frozenset({993, 994, 995, 996, 997})

# Convenience unions
KO_SHEET_8_9: frozenset = KO_SHEET_8 | KO_SHEET_9       # Active face
KO_SHEET_10_11: frozenset = KO_SHEET_10 | KO_SHEET_11   # Neutral
ACTIVE_FACE_RESIDUES: frozenset = KO_SHEET_8_9
BACK_FACE_RESIDUES: frozenset = KO_SHEET_10 | KO_SHEET_11 | KO_SHEET_12
