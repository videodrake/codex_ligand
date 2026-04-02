# -*- coding: utf-8 -*-
import logging
import traceback
from typing import TYPE_CHECKING, Any, Dict, Set, Tuple

from pyrosetta import pose_from_pdb, create_score_function
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.docking import DockMCMProtocol, DockingSlideIntoContact, setup_foldtree
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover
from pyrosetta.rosetta.utility import vector1_int

from . import pyrosetta_init

if TYPE_CHECKING:
    from pyrosetta.rosetta.core.pose import Pose

# ---- Module Constants ---- #
DOCK_JUMP_INDEX = 1
GLOBAL_ROT_MAG = 360
GLOBAL_TRANS_MAG = 100
DEFAULT_EXCLUSION_DIST = 10.0
MIN_CHAINS = 2
SCORE_FUNCTION_NAME = "ref2015"
DOCK_PARTNERS = "A_B"

from .logging_config import get_worker_logger
internal_logger = get_worker_logger()


def _setup_docking(pose: "Pose") -> None:
    """Setup FoldTree and docking jump for 2-chain pose."""
    movable_jumps = vector1_int()
    movable_jumps.append(DOCK_JUMP_INDEX)
    setup_foldtree(pose, DOCK_PARTNERS, movable_jumps)


def run_relax_task(pdb_path: str) -> Dict[str, Any]:
    """
    Runs FastRelax on the input structure.
    Returns a dictionary containing status and PDB data.
    """
    pyrosetta_init.init_rosetta()

    try:
        internal_logger.info(f"[Relax] Starting FastRelax on {pdb_path}")
        pose = pose_from_pdb(pdb_path)
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)

        relax = FastRelax()
        relax.set_scorefxn(scorefxn)
        relax.apply(pose)

        internal_logger.info(f"[Relax] FastRelax completed successfully.")
        return {
            "status": "success",
            "pdb_data": pyrosetta_init.pose_to_string(pose),
            "error": None
        }
    except Exception as e:
        internal_logger.error(f"[Relax] FastRelax failed: {e}")
        internal_logger.debug(traceback.format_exc())
        return {
            "status": "error",
            "pdb_data": None,
            "error": str(e),
            "trace": traceback.format_exc()
        }

def run_global_docking_task(args: tuple) -> Dict[str, Any]:
    """
    Performs global blind docking (Randomize -> SlideIntoContact -> Docking).
    args: (index, relaxed_pdb_string)
       or (index, relaxed_pdb_string, excluded_residues_A, exclusion_contact_dist, max_excluded_contacts)
    """
    pyrosetta_init.init_rosetta()

    if len(args) >= 5:
        idx, input_data, excluded_residues_A, exclusion_contact_dist, max_excluded_contacts = args[:5]
    elif len(args) == 4:
        idx, input_data, excluded_residues_A, exclusion_contact_dist = args
        max_excluded_contacts = 0
    elif len(args) == 3:
        idx, input_data, excluded_residues_A = args
        exclusion_contact_dist = DEFAULT_EXCLUSION_DIST
        max_excluded_contacts = 0
    else:
        idx, input_data = args
        excluded_residues_A = set()
        exclusion_contact_dist = DEFAULT_EXCLUSION_DIST
        max_excluded_contacts = 0

    step = "Init"

    try:
        pose = pyrosetta_init.string_to_pose(input_data)
        if pose is None: raise ValueError("Invalid PDB String")

        # 0. Setup Fold Tree for Docking
        step = "FoldTree Setup"
        _setup_docking(pose)

        # 1. Randomize Position (Global Search)
        step = "Randomization"
        pert_mover = RigidBodyPerturbMover(DOCK_JUMP_INDEX, GLOBAL_ROT_MAG, GLOBAL_TRANS_MAG)
        pert_mover.apply(pose)

        # 2. Slide Into Contact (bring chains together)
        step = "Slide Into Contact"
        slide = DockingSlideIntoContact(DOCK_JUMP_INDEX)
        slide.apply(pose)

        # 2.1. Early rejection: check if Chain B contacts excluded residues
        # This runs BEFORE the expensive DockMCMProtocol to save compute
        if excluded_residues_A:
            step = "Exclusion Check"
            n_contacts = pyrosetta_init.count_excluded_contacts(
                pose, excluded_residues_A, exclusion_contact_dist)
            if n_contacts > max_excluded_contacts:
                return {
                    "status": "rejected",
                    "id": idx,
                    "reason": "excluded_zone_contact",
                    "excluded_contacts": n_contacts,
                    "error": None
                }

        # 3. Docking Protocol (High-Resolution)
        step = "Docking Protocol"
        docking = DockMCMProtocol()
        docking.set_scorefxn(create_score_function(SCORE_FUNCTION_NAME))
        docking.apply(pose)

        # 4. Get Center of Mass
        step = "Center Calculation"
        center = (0, 0, 0)

        if pose.num_chains() >= MIN_CHAINS:
            conf = pose.conformation()
            chain_b_start = conf.chain_begin(2)
            chain_b_end = conf.chain_end(2)

            xyz = xyzVector_double_t(0, 0, 0)
            count = 0

            for i in range(chain_b_start, chain_b_end + 1):
                if pose.residue(i).has("CA"):
                    xyz += pose.residue(i).xyz("CA")
                    count += 1

            if count > 0:
                xyz /= count
                center = (xyz.x, xyz.y, xyz.z)

        return {
            "status": "success",
            "id": idx,
            "pdb_data": pyrosetta_init.pose_to_string(pose),
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "error": None
        }

    except Exception as e:
        internal_logger.error(f"[Docking] Task {idx} failed at step '{step}': {e}")
        internal_logger.debug(traceback.format_exc())
        return {
            "status": "error",
            "id": idx,
            "error": str(e),
            "step": step,
            "trace": traceback.format_exc()
        }

def run_refinement_task(args: tuple) -> Dict[str, Any]:
    """
    Performs local refinement using DockMCMProtocol (Monte Carlo Minimization).
    Includes rigid-body perturbation, sidechain repacking, and minimization
    with Metropolis acceptance criterion — far more effective than simple
    perturb+minimize.
    args: (index, pdb_string, parent_name, trans_mag, rot_mag, sub_rank)
    """
    pyrosetta_init.init_rosetta()

    idx, pdb_data, parent_name, t_mag, r_mag, sub_rank = args
    step = "Init"

    try:
        step = "Loading Pose"
        pose = pyrosetta_init.string_to_pose(pdb_data)
        if pose is None: raise ValueError("Invalid PDB String")

        # Setup FoldTree (lost during deserialization)
        step = "FoldTree Setup"
        _setup_docking(pose)

        # DockMCMProtocol: iterative rigid-body MCM with repacking
        step = "DockMCMProtocol"
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)
        docking = DockMCMProtocol()
        docking.set_scorefxn(scorefxn)
        docking.apply(pose)

        return {
            "status": "success",
            "id": idx,
            "parent": parent_name,
            "sub_rank": sub_rank,
            "pdb_data": pyrosetta_init.pose_to_string(pose),
            "error": None
        }
    except Exception as e:
        internal_logger.error(f"[Refinement] Task {idx} ({parent_name}) failed at step '{step}': {e}")
        internal_logger.debug(traceback.format_exc())
        return {
            "status": "error",
            "id": idx,
            "error": str(e),
            "step": step,
            "trace": traceback.format_exc()
        }


def run_mini_refinement_task(args: tuple) -> Dict[str, Any]:
    """
    Lightweight refinement for Stage 1 survivors.
    Repacks (and optionally minimizes) interface sidechains to improve
    packing quality and shape complementarity before Stage 2 filtering.

    args: (idx, pdb_data, mode)            — legacy 3-element
       or (idx, pdb_data, mode, n_rounds)  — enhanced 4-element
        idx:      decoy index
        pdb_data: PDB string
        mode:     'pack_only' = interface repacking only
                  'pack_min'  = repacking + interface minimization
        n_rounds: number of repacking rounds (default 3)
    """
    pyrosetta_init.init_rosetta()

    if len(args) >= 4:
        idx, pdb_data, mode, n_rounds = args[:4]
    else:
        idx, pdb_data, mode = args
        n_rounds = 3
    step = "Init"

    try:
        step = "Loading Pose"
        pose = pyrosetta_init.string_to_pose(pdb_data)
        if pose is None:
            raise ValueError("Invalid PDB String")

        step = "FoldTree Setup"
        _setup_docking(pose)

        step = "Score Function"
        scorefxn = create_score_function(SCORE_FUNCTION_NAME)

        # --- Select interface residues ---
        step = "Interface Selection"
        try:
            from pyrosetta.rosetta.core.select.residue_selector import (
                InterGroupInterfaceByVectorSelector,
                ChainSelector,
            )
            chain_a_sel = ChainSelector("A")
            chain_b_sel = ChainSelector("B")
            iface_sel = InterGroupInterfaceByVectorSelector(chain_a_sel, chain_b_sel)
            iface_mask = iface_sel.apply(pose)
        except Exception as e_sel:
            # Fallback: select residues within 10A of the other chain
            internal_logger.debug(f"[MiniRefine] InterGroupInterfaceByVectorSelector unavailable: {e_sel}")
            iface_mask = _fallback_interface_mask(pose)

        # --- Repacking ---
        step = "Repacking"
        try:
            from pyrosetta.rosetta.core.pack.task import TaskFactory
            from pyrosetta.rosetta.core.pack.task.operation import (
                RestrictToRepacking,
                PreventRepacking,
            )
            from pyrosetta.rosetta.protocols.minimization_packing import (
                PackRotamersMover,
            )

            tf = TaskFactory()
            tf.push_back(RestrictToRepacking())

            # Include current rotamers as candidates (preserve good ones)
            try:
                from pyrosetta.rosetta.core.pack.task.operation import IncludeCurrent
                tf.push_back(IncludeCurrent())
            except Exception as e_ic:
                internal_logger.debug(f"[MiniRefine] IncludeCurrent unavailable: {e_ic}")

            # Expand chi sampling with extra rotamers (ex1 + ex2)
            try:
                from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamersGeneric
                ex = ExtraRotamersGeneric()
                ex.ex1(True)
                ex.ex2(True)
                tf.push_back(ex)
            except Exception as e_er:
                internal_logger.debug(f"[MiniRefine] ExtraRotamersGeneric unavailable: {e_er}")

            # Prevent repacking of non-interface residues
            prevent = PreventRepacking()
            for i in range(1, pose.total_residue() + 1):
                if not iface_mask[i]:
                    prevent.include_residue(i)
            tf.push_back(prevent)

            packer = PackRotamersMover(scorefxn)
            packer.task_factory(tf)
            for _ in range(n_rounds):
                packer.apply(pose)
        except Exception as e_pack:
            internal_logger.warning(
                f"[MiniRefine] Task {idx} repacking failed, returning original: {e_pack}")
            return {
                "status": "success",
                "id": idx,
                "pdb_data": pdb_data,
                "refined": False,
                "error": None,
            }

        # --- Optional Minimization ---
        if mode == "pack_min":
            step = "Minimization"
            try:
                from pyrosetta.rosetta.core.kinematics import MoveMap
                from pyrosetta.rosetta.protocols.minimization_packing import (
                    MinMover,
                )

                mm = MoveMap()
                mm.set_bb(False)
                mm.set_chi(False)
                mm.set_jump(False)  # keep rigid-body fixed
                for i in range(1, pose.total_residue() + 1):
                    if iface_mask[i]:
                        mm.set_chi(i, True)
                        # backbone stays fixed to prevent rigid-body drift

                min_mover = MinMover()
                min_mover.movemap(mm)
                min_mover.score_function(scorefxn)
                min_mover.min_type("lbfgs_armijo_nonmonotone")
                min_mover.tolerance(0.01)
                min_mover.apply(pose)
            except Exception as e_min:
                internal_logger.warning(
                    f"[MiniRefine] Task {idx} minimization failed (pack OK): {e_min}")

        return {
            "status": "success",
            "id": idx,
            "pdb_data": pyrosetta_init.pose_to_string(pose),
            "refined": True,
            "error": None,
        }

    except Exception as e:
        internal_logger.error(
            f"[MiniRefine] Task {idx} failed at step '{step}': {e}")
        internal_logger.debug(traceback.format_exc())
        return {
            "status": "error",
            "id": idx,
            "pdb_data": pdb_data,
            "refined": False,
            "error": str(e),
            "step": step,
            "trace": traceback.format_exc(),
        }


def _fallback_interface_mask(pose: "Pose", dist_cutoff: float = 10.0):
    """Fallback interface detection when InterGroupInterfaceByVectorSelector unavailable."""
    from pyrosetta.rosetta.utility import vector1_bool
    n = pose.total_residue()
    mask = vector1_bool(n, False)

    conf = pose.conformation()
    chain_a_end = conf.chain_end(1)
    chain_b_start = conf.chain_begin(2)
    chain_b_end = conf.chain_end(2)

    for i in range(1, chain_a_end + 1):
        if not pose.residue(i).has("CA"):
            continue
        ca_i = pose.residue(i).xyz("CA")
        for j in range(chain_b_start, chain_b_end + 1):
            if not pose.residue(j).has("CA"):
                continue
            if ca_i.distance(pose.residue(j).xyz("CA")) < dist_cutoff:
                mask[i] = True
                mask[j] = True
                break

    for j in range(chain_b_start, chain_b_end + 1):
        if mask[j]:
            continue
        if not pose.residue(j).has("CA"):
            continue
        ca_j = pose.residue(j).xyz("CA")
        for i in range(1, chain_a_end + 1):
            if not pose.residue(i).has("CA"):
                continue
            if ca_j.distance(pose.residue(i).xyz("CA")) < dist_cutoff:
                mask[j] = True
                break

    return mask
