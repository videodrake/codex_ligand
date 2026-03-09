# -*- coding: utf-8 -*-
import logging
import math
import os
import traceback
import re
from typing import TYPE_CHECKING, Any, Dict, Optional, Set, Tuple, Union

import pyrosetta.rosetta.core.pose
from pyrosetta import pose_from_pdb, create_score_function, get_score_function
from pyrosetta.rosetta.core.pose import getPoseExtraScore
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol, fa_elec, calpha_superimpose_pose
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

import common

if TYPE_CHECKING:
    from pyrosetta.rosetta.core.pose import Pose

# ---- Module Constants ---- #
DEFAULT_CONTACT_DISTANCE = 6.0
DEFAULT_EXCLUSION_DISTANCE = 10.0
BFACTOR_CLAMP_MIN = -5.0
BFACTOR_CLAMP_MAX = 10.0
RMSD_ERROR_VALUE = 999.9
SCORE_FUNCTION_NAME = "ref2015"
MIN_CHAINS = 2

internal_logger = logging.getLogger('python_internal')


# ---- Private Helpers ---- #

def _extract_iam_metric(pose: "Pose", iam: Any, score_key: str, method_name: str, default: Union[int, float]) -> Union[int, float]:
    """Extract IAM metric with 2-level fallback.

    1st: getPoseExtraScore (server-tested OK for all metrics)
    2nd: iam method (hasattr guard)
    3rd: safe default (conservative: fails filter direction)
    """
    val = default
    try:
        raw = getPoseExtraScore(pose, score_key)
        val = type(default)(raw) if isinstance(default, int) else float(raw)
    except Exception:
        if hasattr(iam, method_name):
            try:
                raw = getattr(iam, method_name)()
                val = type(default)(raw) if isinstance(default, int) else float(raw)
            except Exception:
                pass
    return val


def _load_pose(target: str) -> Tuple[Optional["Pose"], str]:
    """Load Pose from file path or PDB string."""
    if os.path.exists(str(target)) and len(str(target)) < 255:
        return pose_from_pdb(target), target
    else:
        return common.string_to_pose(target), "memory_pose"


def _get_interface_def(pose: "Pose") -> str:
    """Build interface definition string from pose chains (e.g. 'A_B')."""
    conf = pose.conformation()
    c1 = pose.pdb_info().chain(conf.chain_begin(1))
    c2 = pose.pdb_info().chain(conf.chain_begin(2))
    return f"{c1}_{c2}"


def inject_bfactors(pose: "Pose") -> None:
    """ 
    Injects residue energies into B-factor column.
    Iterates through atoms to support all PyRosetta versions.
    """
    if not pose.pdb_info():
        pose.pdb_info(pyrosetta.rosetta.core.pose.PDBInfo(pose))
    scorefxn = get_score_function() 
    scorefxn(pose)
    
    pdb_info = pose.pdb_info()
    
    for i in range(1, pose.total_residue() + 1):
        try:
            e = pose.energies().residue_total_energy(i)
            e_clamped = max(BFACTOR_CLAMP_MIN, min(BFACTOR_CLAMP_MAX, e))
            
            res = pose.residue(i)
            for atom_idx in range(1, res.natoms() + 1):
                pdb_info.bfactor(i, atom_idx, e_clamped)
        except Exception as e:
            internal_logger.warning(f"Error injecting B-factor for residue {i}: {e}")
            internal_logger.debug(traceback.format_exc())

def get_residue_energy_csv_string(pose: "Pose", scorefxn: Any) -> str:
    """
    Generates per-residue energy breakdown.
    Uses residue_total_energies() for per-term EMapVector access.
    """
    weights = scorefxn.weights()
    lines = ["Residue_ID,Residue_Name,Total_Energy,fa_atr,fa_rep,fa_sol,fa_elec"]

    has_details = hasattr(pose.energies(), 'residue_total_energies')

    for i in range(1, pose.total_residue() + 1):
        name = pose.residue(i).name3()
        pdb_num = pose.pdb_info().number(i)
        chain = pose.pdb_info().chain(i)
        total = pose.energies().residue_total_energy(i)

        atr, rep, sol, elec = 0.0, 0.0, 0.0, 0.0

        if has_details:
            try:
                emap = pose.energies().residue_total_energies(i)
                atr = emap[fa_atr] * weights[fa_atr]
                rep = emap[fa_rep] * weights[fa_rep]
                sol = emap[fa_sol] * weights[fa_sol]
                elec = emap[fa_elec] * weights[fa_elec]
            except Exception as e:
                internal_logger.warning(f"Error getting residue energies for residue {i}: {e}")
                internal_logger.debug(traceback.format_exc())

        lines.append(f"{chain}{pdb_num},{name},{total:.4f},{atr:.4f},{rep:.4f},{sol:.4f},{elec:.4f}")
    return "\n".join(lines)

def get_interface_energy_csv_string(pose: "Pose", scorefxn: Any, threshold: float = 0.5) -> str:
    """
    Generates per-residue interface energy breakdown (DeltaE).
    DeltaE(i) = E_complex(i) - E_separated(i).
    Only residues with |DeltaE| > threshold are included.
    Uses split_by_chain() to obtain separated chains.
    Returns CSV string or empty string on failure.
    """
    try:
        # 1. Score complex
        scorefxn(pose)
        weights = scorefxn.weights()

        total_res = pose.total_residue()
        complex_energies = {}
        has_details = hasattr(pose.energies(), 'residue_total_energies')

        for i in range(1, total_res + 1):
            total_e = pose.energies().residue_total_energy(i)
            atr, rep, sol, elec = 0.0, 0.0, 0.0, 0.0
            if has_details:
                try:
                    emap = pose.energies().residue_total_energies(i)
                    atr = emap[fa_atr] * weights[fa_atr]
                    rep = emap[fa_rep] * weights[fa_rep]
                    sol = emap[fa_sol] * weights[fa_sol]
                    elec = emap[fa_elec] * weights[fa_elec]
                except Exception:
                    pass
            chain = pose.pdb_info().chain(i)
            pdb_num = pose.pdb_info().number(i)
            res_name = pose.residue(i).name3()
            complex_energies[i] = {
                'chain': chain, 'pdb_num': pdb_num, 'name': res_name,
                'total': total_e, 'fa_atr': atr, 'fa_rep': rep,
                'fa_sol': sol, 'fa_elec': elec,
            }

        # 2. Split by chain and score individually
        try:
            chains = pose.split_by_chain()
        except Exception:
            return ""

        if len(chains) < 2:
            return ""

        sep_energies = {}
        rosetta_offset = 0

        for chain_idx in range(1, len(chains) + 1):
            chain_pose = chains[chain_idx]
            scorefxn(chain_pose)
            chain_has_details = hasattr(chain_pose.energies(), 'residue_total_energies')

            for j in range(1, chain_pose.total_residue() + 1):
                orig_i = rosetta_offset + j
                total_e = chain_pose.energies().residue_total_energy(j)
                atr, rep, sol, elec = 0.0, 0.0, 0.0, 0.0
                if chain_has_details:
                    try:
                        emap = chain_pose.energies().residue_total_energies(j)
                        atr = emap[fa_atr] * weights[fa_atr]
                        rep = emap[fa_rep] * weights[fa_rep]
                        sol = emap[fa_sol] * weights[fa_sol]
                        elec = emap[fa_elec] * weights[fa_elec]
                    except Exception:
                        pass
                sep_energies[orig_i] = {
                    'total': total_e, 'fa_atr': atr, 'fa_rep': rep,
                    'fa_sol': sol, 'fa_elec': elec,
                }

            rosetta_offset += chain_pose.total_residue()

        # 3. Calculate DeltaE and filter
        lines = ["Residue_ID,Residue_Name,Chain,DeltaE_total,DeltaE_fa_atr,DeltaE_fa_rep,DeltaE_fa_sol,DeltaE_fa_elec"]

        for i in range(1, total_res + 1):
            if i not in complex_energies or i not in sep_energies:
                continue
            ce = complex_energies[i]
            se = sep_energies[i]
            d_total = ce['total'] - se['total']

            if abs(d_total) < threshold:
                continue

            d_atr = ce['fa_atr'] - se['fa_atr']
            d_rep = ce['fa_rep'] - se['fa_rep']
            d_sol = ce['fa_sol'] - se['fa_sol']
            d_elec = ce['fa_elec'] - se['fa_elec']

            res_id = f"{ce['chain']}{ce['pdb_num']}"
            lines.append(f"{res_id},{ce['name']},{ce['chain']},"
                         f"{d_total:.4f},{d_atr:.4f},{d_rep:.4f},"
                         f"{d_sol:.4f},{d_elec:.4f}")

        if len(lines) <= 1:
            return ""
        return "\n".join(lines)
    except Exception as e:
        internal_logger.warning(f"Error generating interface energy CSV: {e}")
        internal_logger.debug(traceback.format_exc())
        return ""


def analyze_interface_contacts(pose: "Pose", contact_distance: float = DEFAULT_CONTACT_DISTANCE) -> Dict[str, str]:
    """Returns per-chain binding residues with chain ID prefix.

    Returns {'A': 'A:LEU819,A:ALA822', 'B': 'B:VAL962,B:ASN963'}
    Chain A = receptor-side contacts, Chain B = ligand-side contacts.
    """
    try:
        if pose.num_chains() < MIN_CHAINS:
            return {"A": "No_Chain_2", "B": "No_Chain_2"}

        conf = pose.conformation()
        chain1_id = pose.pdb_info().chain(conf.chain_begin(1))
        chain2_id = pose.pdb_info().chain(conf.chain_begin(2))

        # Forward: neighbors of Chain B within Chain A (receptor contacts)
        ligand_sel = ChainSelector(chain2_id)
        fwd_mask = NeighborhoodResidueSelector(ligand_sel, contact_distance, False).apply(pose)
        contacts_a = []
        for i in range(1, len(fwd_mask) + 1):
            if fwd_mask[i] and pose.pdb_info().chain(i) == chain1_id:
                pdb_num = pose.pdb_info().number(i)
                res_name = pose.residue(i).name3()
                contacts_a.append(f"{chain1_id}:{res_name}{pdb_num}")

        # Reverse: neighbors of Chain A within Chain B (ligand contacts)
        receptor_sel = ChainSelector(chain1_id)
        rev_mask = NeighborhoodResidueSelector(receptor_sel, contact_distance, False).apply(pose)
        contacts_b = []
        for i in range(1, len(rev_mask) + 1):
            if rev_mask[i] and pose.pdb_info().chain(i) == chain2_id:
                pdb_num = pose.pdb_info().number(i)
                res_name = pose.residue(i).name3()
                contacts_b.append(f"{chain2_id}:{res_name}{pdb_num}")

        return {
            chain1_id: ",".join(contacts_a) if contacts_a else "None",
            chain2_id: ",".join(contacts_b) if contacts_b else "None",
        }
    except Exception as e:
        internal_logger.error(f"Error analyzing interface contacts: {e}")
        internal_logger.debug(traceback.format_exc())
        return {"A": "Analysis_Failed", "B": "Analysis_Failed"}

count_excluded_contacts = common.count_excluded_contacts


def _derive_key_contact_from_binding(
    binding_res_b_str: str,
    key_residues_B: Set[int],
    key_weights: Optional[Dict[int, float]] = None,
) -> Tuple[float, int, int]:
    """Derive key_contact_ratio directly from Binding_Residues_B string.

    Parses PDB numbers from binding string (e.g. "B:VAL962,B:ASN963")
    and intersects with key_residues_B. Supports optional per-residue weights.
    Returns (ratio, n_contacted, n_present) where n_present = len(key_residues_B).
    """
    if not key_residues_B:
        return (0.0, 0, 0)

    # Parse PDB residue numbers from binding string
    contacted_pdb_nums: Set[int] = set()
    if binding_res_b_str and binding_res_b_str not in ("None", "Analysis_Failed", "No_Chain_2"):
        for token in binding_res_b_str.split(','):
            token = token.strip()
            # Format: "B:VAL962" or "B:ALA1004" — extract trailing digits
            m = re.search(r'(-?\d+)$', token)
            if m:
                contacted_pdb_nums.add(int(m.group(1)))

    n_present = len(key_residues_B)
    if n_present == 0:
        return (0.0, 0, 0)

    if key_weights:
        w_contacted = sum(key_weights.get(r, 1.0) for r in key_residues_B if r in contacted_pdb_nums)
        w_total = sum(key_weights.get(r, 1.0) for r in key_residues_B)
        ratio = w_contacted / w_total if w_total > 0 else 0.0
        n_contacted = sum(1 for r in key_residues_B if r in contacted_pdb_nums)
    else:
        n_contacted = sum(1 for r in key_residues_B if r in contacted_pdb_nums)
        ratio = n_contacted / n_present

    return (ratio, n_contacted, n_present)


def calculate_key_contact_ratio(pose: "Pose", key_residues_B: Set[int], contact_dist: float = DEFAULT_CONTACT_DISTANCE, key_weights: Optional[Dict[int, float]] = None) -> Tuple[float, int, int]:
    """
    Calculate the fraction of key Chain B residues that participate in
    the binding interface using NeighborhoodResidueSelector (all-atom,
    consistent with analyze_interface_contacts).
    Residues not present in the PDB are excluded from the denominator.
    Returns (ratio, n_contacted, n_present).
    """
    if not key_residues_B or pose.num_chains() < MIN_CHAINS:
        return (0.0, 0, 0)

    try:
        conf = pose.conformation()
        chain1_id = pose.pdb_info().chain(conf.chain_begin(1))

        # All-atom neighbor detection: Chain B residues near Chain A
        receptor_sel = ChainSelector(chain1_id)
        rev_mask = NeighborhoodResidueSelector(receptor_sel, contact_dist, False).apply(pose)

        chain_b_start = conf.chain_begin(2)
        chain_b_end = conf.chain_end(2)

        n_present = 0
        n_contacted = 0
        w_present = 0.0
        w_contacted = 0.0

        for i in range(chain_b_start, chain_b_end + 1):
            pdb_num = pose.pdb_info().number(i)
            if pdb_num not in key_residues_B:
                continue
            n_present += 1
            w = key_weights.get(pdb_num, 1.0) if key_weights else 1.0
            w_present += w
            if rev_mask[i]:
                n_contacted += 1
                w_contacted += w

        if key_weights and w_present > 0:
            ratio = w_contacted / w_present
        else:
            ratio = n_contacted / n_present if n_present > 0 else 0.0
        return (ratio, n_contacted, n_present)
    except Exception:
        internal_logger.warning("NeighborhoodResidueSelector failed in key_contact_ratio, using legacy CA-CA")
        return _calculate_key_contact_ratio_legacy(pose, key_residues_B, contact_dist)


def _calculate_key_contact_ratio_legacy(pose: "Pose", key_residues_B: Set[int], contact_dist: float) -> Tuple[float, int, int]:
    """Legacy CA-CA distance fallback for calculate_key_contact_ratio."""
    conf = pose.conformation()
    chain_a_end = conf.chain_end(1)
    chain_b_start = conf.chain_begin(2)
    chain_b_end = conf.chain_end(2)

    a_cas = []
    for i in range(1, chain_a_end + 1):
        if pose.residue(i).has("CA"):
            a_cas.append(pose.residue(i).xyz("CA"))

    if not a_cas:
        return (0.0, 0, 0)

    dist_sq = contact_dist * contact_dist
    n_present = 0
    n_contacted = 0

    for i in range(chain_b_start, chain_b_end + 1):
        pdb_num = pose.pdb_info().number(i)
        if pdb_num not in key_residues_B:
            continue
        if not pose.residue(i).has("CA"):
            continue

        n_present += 1
        b_ca = pose.residue(i).xyz("CA")
        for a_ca in a_cas:
            dx = a_ca.x - b_ca.x
            dy = a_ca.y - b_ca.y
            dz = a_ca.z - b_ca.z
            if dx * dx + dy * dy + dz * dz <= dist_sq:
                n_contacted += 1
                break

    ratio = n_contacted / n_present if n_present > 0 else 0.0
    return (ratio, n_contacted, n_present)


def get_chain_b_center(pose: "Pose") -> Tuple[float, float, float]:
    """ Calculates geometric center of Chain B safely. """
    if pose.num_chains() < MIN_CHAINS: return (0,0,0)

    conf = pose.conformation()
    start = conf.chain_begin(2)
    end = conf.chain_end(2)
    
    sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
    count = 0

    for i in range(start, end+1):
        if pose.residue(i).has("CA"):
            coords = pose.residue(i).xyz("CA")
            sum_x += coords.x
            sum_y += coords.y
            sum_z += coords.z
            count += 1

    if count == 0: return (0,0,0)
    return (sum_x/count, sum_y/count, sum_z/count)

def calculate_l_rmsd(mobile: "Pose", ref: "Pose") -> float:
    """ Calculates Ligand RMSD relative to reference. """
    try:
        calpha_superimpose_pose(mobile, ref)
        
        chain_a_end = mobile.conformation().chain_end(1)
        chain_b_start = chain_a_end + 1
        chain_b_end = mobile.total_residue()
        
        rmsd_sum = 0.0
        count = 0

        for i in range(chain_b_start, chain_b_end + 1):
            if mobile.residue(i).has("CA"):
                d2 = mobile.residue(i).xyz("CA").distance_squared(ref.residue(i).xyz("CA"))
                rmsd_sum += d2
                count += 1

        if count == 0: return RMSD_ERROR_VALUE
        return math.sqrt(rmsd_sum / count)
    except Exception as e:
        internal_logger.error(f"Error calculating L_RMSD: {e}")
        internal_logger.debug(traceback.format_exc())
        return RMSD_ERROR_VALUE

def run_lrmsd_dual_task(args: tuple) -> Dict[str, Any]:
    """Compute L_RMSD against both relaxed reference and best-dG reference.
    args: (pdb_data_string, relaxed_ref_path, best_ref_path)
    """
    common.init_rosetta()
    pdb_data, relaxed_ref_path, best_ref_path = args
    try:
        pose = common.string_to_pose(pdb_data)
        if pose is None or pose.num_chains() < MIN_CHAINS:
            return {"status": "error", "L_RMSD": RMSD_ERROR_VALUE,
                    "L_RMSD_best": RMSD_ERROR_VALUE}

        l_rmsd = RMSD_ERROR_VALUE
        if relaxed_ref_path and os.path.exists(relaxed_ref_path):
            try:
                ref_pose = pose_from_pdb(relaxed_ref_path)
                l_rmsd = calculate_l_rmsd(pose.clone(), ref_pose)
            except Exception:
                pass

        l_rmsd_best = RMSD_ERROR_VALUE
        if best_ref_path and os.path.exists(best_ref_path):
            try:
                best_pose = pose_from_pdb(best_ref_path)
                l_rmsd_best = calculate_l_rmsd(pose.clone(), best_pose)
            except Exception:
                l_rmsd_best = RMSD_ERROR_VALUE

        return {"status": "success", "L_RMSD": l_rmsd, "L_RMSD_best": l_rmsd_best}
    except Exception:
        return {"status": "error", "L_RMSD": RMSD_ERROR_VALUE,
                "L_RMSD_best": RMSD_ERROR_VALUE}


def run_fast_scoring_task(args: tuple) -> Dict[str, Any]:
    """
    Fast scoring for filtering stage.
    Computes only: dG, dSASA, sc_value, Chain B center, constraint metrics.
    Skips: packstat, L_RMSD, B-factors, per-residue CSV, binding residues.
    args: (target, ref_path, contact_dist)
       or (target, ref_path, contact_dist, constraints)
    constraints dict: {excluded_residues_A, key_residues_B, exclusion_contact_dist}
    """
    common.init_rosetta()

    constraints = {}
    if len(args) >= 4:
        target, ref_path, contact_dist, constraints = args[:4]
    elif len(args) == 3:
        target, ref_path, contact_dist = args
    else:
        target, ref_path = args
        contact_dist = DEFAULT_CONTACT_DISTANCE

    excluded_A = constraints.get('excluded_residues_A', set())
    key_B = constraints.get('key_residues_B', set())
    excl_dist = constraints.get('exclusion_contact_dist', DEFAULT_EXCLUSION_DISTANCE)
    step_status = "Init"

    try:
        step_status = "Loading Pose"
        pose, _ = _load_pose(target)

        if pose is None or pose.num_chains() < MIN_CHAINS:
            return {"status": "error", "error": "Invalid Pose or Chain count < 2"}

        step_status = "Detecting Chains"
        interface_def = _get_interface_def(pose)

        step_status = "Scoring"
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)
        total_score = scorefxn(pose)

        step_status = "InterfaceAnalyzer"
        iam = InterfaceAnalyzerMover()
        iam.set_interface(interface_def)
        iam.set_scorefunction(scorefxn)
        iam.set_pack_separated(True)
        iam.set_compute_packstat(False)
        compute_sc = constraints.get('compute_sc', True)  # default True (backward compat)
        iam.set_compute_interface_sc(compute_sc)
        iam.apply(pose)

        sc_val = _extract_iam_metric(pose, iam, "sc_value", "get_sc_value", 0.0)

        step_status = "Center"
        center = get_chain_b_center(pose)

        # Constraint metrics
        step_status = "Constraint Metrics"
        excl_count = 0
        key_ratio, key_contacted, key_present = 0.0, 0, 0

        if excluded_A:
            excl_count = count_excluded_contacts(pose, excluded_A, excl_dist)
        if key_B:
            key_weights = constraints.get('key_residue_weights')
            key_ratio, key_contacted, key_present = calculate_key_contact_ratio(
                pose, key_B, contact_dist, key_weights)

        return {
            "status": "success",
            "dG_separated": iam.get_interface_dG(),
            "dSASA": iam.get_interface_delta_sasa(),
            "sc_value": sc_val,
            "total_score": total_score,
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "excluded_contacts": excl_count,
            "key_contact_ratio": key_ratio,
            "key_contacted": key_contacted,
            "key_present": key_present,
        }
    except Exception as e:
        return {
            "status": "error",
            "error": str(e),
            "step": step_status,
            "trace": traceback.format_exc()
        }


def run_intermediate_scoring_task(args: tuple) -> Dict[str, Any]:
    """
    Intermediate scoring for v2.0 Stage 2 expensive metrics.
    Computes: packstat, delta_unsatHbonds, nres_int.
    Called only on Stage 1 survivors to save compute time.
    args: (pdb_data_string,)
    """
    common.init_rosetta()

    if isinstance(args, tuple):
        pdb_data_string = args[0]
    else:
        pdb_data_string = args

    step_status = "Init"

    try:
        step_status = "Loading Pose"
        pose = common.string_to_pose(pdb_data_string)

        if pose is None or pose.num_chains() < MIN_CHAINS:
            return {"status": "error", "error": "Invalid Pose or Chain count < 2"}

        step_status = "Detecting Chains"
        interface_def = _get_interface_def(pose)

        step_status = "Scoring"
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)
        scorefxn(pose)

        step_status = "InterfaceAnalyzer (expensive)"
        iam = InterfaceAnalyzerMover()
        iam.set_interface(interface_def)
        iam.set_scorefunction(scorefxn)
        iam.set_pack_separated(True)
        iam.set_compute_packstat(True)
        iam.set_compute_interface_sc(False)
        iam.apply(pose)

        packstat_val = _extract_iam_metric(
            pose, iam, "packstat", "get_interface_packstat", 0.0)
        delta_unsatHbonds = _extract_iam_metric(
            pose, iam, "delta_unsatHbonds", "get_interface_delta_hbond_unsat", 99)
        nres_int = _extract_iam_metric(
            pose, iam, "nres_int", "get_num_interface_residues", 0)
        hbonds_int = _extract_iam_metric(
            pose, iam, "hbonds_int", "get_interface_hbonds", 0)

        return {
            "status": "success",
            "packstat": packstat_val,
            "delta_unsatHbonds": delta_unsatHbonds,
            "nres_int": nres_int,
            "hbonds_int": hbonds_int,
        }
    except Exception as e:
        return {
            "status": "error",
            "error": str(e),
            "step": step_status,
            "trace": traceback.format_exc()
        }


def run_scoring_task(args: tuple) -> Dict[str, Any]:
    """
    Full scoring task. Calculates Score, RMSD, SC, Contacts, B-factors, CSV.
    args: (target, ref_path, contact_dist)
       or (target, ref_path, contact_dist, constraints)
    """
    common.init_rosetta()

    constraints = {}
    if len(args) >= 4:
        target, ref_path, contact_dist, constraints = args[:4]
    elif len(args) == 3:
        target, ref_path, contact_dist = args
    else:
        target, ref_path = args
        contact_dist = DEFAULT_CONTACT_DISTANCE

    key_B = constraints.get('key_residues_B', set())
    step_status = "Init"
    interface_def = "Unknown"

    try:
        step_status = "Loading Pose"
        pose, source_name = _load_pose(target)

        if pose is None or pose.num_chains() < MIN_CHAINS:
            return {"error": "Invalid Pose or Chain count < 2", "step": step_status, "status": "error"}

        step_status = "Detecting Chains"
        interface_def = _get_interface_def(pose)

        step_status = "Calculating RMSD"
        l_rmsd = 0.0
        if ref_path and os.path.exists(ref_path):
            try:
                ref_pose = pose_from_pdb(ref_path)
                l_rmsd = calculate_l_rmsd(pose, ref_pose)
            except Exception as e:
                internal_logger.warning(f"Error calculating L_RMSD for pose {source_name} against ref {ref_path}: {e}")
                internal_logger.debug(traceback.format_exc())

        step_status = "Calculating ScoreFunction"
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)
        total_score = scorefxn(pose)

        step_status = f"Running InterfaceAnalyzer ({interface_def})"
        iam = InterfaceAnalyzerMover()
        iam.set_interface(interface_def)
        iam.set_scorefunction(scorefxn)
        iam.set_pack_separated(True)
        iam.set_compute_packstat(True)
        iam.set_compute_interface_sc(True)
        iam.apply(pose)

        sc_val = _extract_iam_metric(pose, iam, "sc_value", "get_sc_value", 0.0)
        packstat_val = _extract_iam_metric(
            pose, iam, "packstat", "get_interface_packstat", 0.0)
        delta_unsatHbonds = _extract_iam_metric(
            pose, iam, "delta_unsatHbonds", "get_interface_delta_hbond_unsat", 99)
        nres_int = _extract_iam_metric(
            pose, iam, "nres_int", "get_num_interface_residues", 0)
        hbonds_int = _extract_iam_metric(
            pose, iam, "hbonds_int", "get_interface_hbonds", 0)

        step_status = "Injecting B-factors"
        inject_bfactors(pose)

        step_status = "Generating CSV Data"
        csv_data = get_residue_energy_csv_string(pose, scorefxn)

        step_status = "Analyzing Contacts"
        binding_res_dict = analyze_interface_contacts(pose, contact_dist)

        step_status = "Generating Interface Energy CSV"
        interface_csv_data = get_interface_energy_csv_string(pose, scorefxn)

        step_status = "Calculating Center"
        center = get_chain_b_center(pose)

        # Key residue contact metric (soft) — derive from Binding_Residues_B
        step_status = "Key Residue Metric"
        key_ratio, key_contacted, key_present = 0.0, 0, 0
        key_weights = constraints.get('key_residue_weights')
        if key_B:
            br_b = binding_res_dict.get("B", "")
            if br_b and br_b not in ("None", "Analysis_Failed", "No_Chain_2"):
                key_ratio, key_contacted, key_present = _derive_key_contact_from_binding(
                    br_b, key_B, key_weights)
            else:
                # Fallback: direct pose-based calculation
                key_ratio, key_contacted, key_present = calculate_key_contact_ratio(
                    pose, key_B, contact_dist, key_weights)

        step_status = "Saving String"
        aligned_pdb_string = common.pose_to_string(pose)

        return {
            "status": "success",
            "File": source_name,
            "dG_separated": iam.get_interface_dG(),
            "Total_Score": total_score,
            "dSASA": iam.get_interface_delta_sasa(),
            "sc_value": sc_val,
            "packstat": packstat_val,
            "delta_unsatHbonds": delta_unsatHbonds,
            "nres_int": nres_int,
            "hbonds_int": hbonds_int,
            "L_RMSD": l_rmsd,
            "Binding_Residues": binding_res_dict.get("A", ""),
            "Binding_Residues_A": binding_res_dict.get("A", ""),
            "Binding_Residues_B": binding_res_dict.get("B", ""),
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "key_contact_ratio": key_ratio,
            "key_contacted": key_contacted,
            "key_present": key_present,
            "pdb_data": aligned_pdb_string,
            "csv_data": csv_data,
            "interface_csv_data": interface_csv_data,
        }
    except Exception as e:
        return {
            "error": str(e),
            "status": "error",
            "step": step_status,
            "chains": interface_def,
            "trace": traceback.format_exc()
        }