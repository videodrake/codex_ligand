#!/usr/bin/env python
"""Restore V924R mutation to wild-type VAL in PPI dimer PDBs.

3GT8 crystal structure contains V924R (PDB resnum 948) engineering mutation
(Jura et al. 2009 Cell) for crystallization. This must be reverted to
wild-type VAL before PPI docking, since MYO1D binds wild-type EGFR.

Usage (on server with PyRosetta):
    conda activate pyrosetta
    python scripts/fix_v924r_mutation.py

This script:
  1. Loads each dimer PDB
  2. Mutates chain A residue 948 ARG → VAL using PyRosetta MutateResidue
  3. Repacks neighbors (8Å shell) for proper sidechain placement
  4. Saves to *_wt.pdb (original files preserved)
  5. Verifies the mutation was applied
"""
import sys
import os
from pathlib import Path


def init_pyrosetta():
    """Initialize PyRosetta with minimal output."""
    import pyrosetta
    pyrosetta.init(
        "-mute all -ignore_unrecognized_res true",
        set_logging_handler=None,
    )


def mutate_and_repack(pose, chain_id, resnum, new_aa, repack_radius=8.0):
    """Mutate a residue and repack its neighborhood.

    Args:
        pose: PyRosetta Pose object
        chain_id: Chain letter (e.g., 'A')
        resnum: PDB residue number
        new_aa: One-letter code of target amino acid (e.g., 'V' for VAL)
        repack_radius: Radius for neighbor repacking (Angstroms)
    """
    from pyrosetta.rosetta.core.pose import PDBInfo
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking,
        PreventRepacking,
    )
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.core.select.residue_selector import (
        NeighborhoodResidueSelector,
        ResidueIndexSelector,
    )

    # Find Rosetta internal residue index for this PDB residue
    pdb_info = pose.pdb_info()
    pose_idx = pdb_info.pdb2pose(chain_id, resnum)
    if pose_idx == 0:
        raise ValueError(f"Residue {chain_id}:{resnum} not found in pose")

    old_name = pose.residue(pose_idx).name3()
    print(f"  Residue {chain_id}:{resnum} = {old_name} (pose index {pose_idx})")

    # Step 1: Mutate
    mutator = MutateResidue(pose_idx, new_aa)
    mutator.apply(pose)
    new_name = pose.residue(pose_idx).name3()
    print(f"  Mutated: {old_name} → {new_name}")

    # Step 2: Repack neighborhood
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")

    # Select the mutated residue and its neighbors
    mut_selector = ResidueIndexSelector(str(pose_idx))
    nbr_selector = NeighborhoodResidueSelector(mut_selector, repack_radius, True)

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())

    # Only repack residues near the mutation site
    packer = PackRotamersMover(sfxn)
    packer.task_factory(tf)

    # Restrict to neighborhood only
    from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset
    from pyrosetta.rosetta.core.pack.task.operation import PreventRepackingRLT

    prevent_rlt = PreventRepackingRLT()
    from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
    not_nbr = NotResidueSelector(nbr_selector)
    restrict_op = OperateOnResidueSubset(prevent_rlt, not_nbr)
    tf.push_back(restrict_op)

    packer = PackRotamersMover(sfxn)
    packer.task_factory(tf)
    packer.apply(pose)
    print(f"  Repacked {repack_radius}Å neighborhood")

    return pose


def process_pdb(input_path, output_path, chain_id="A", resnum=948, new_aa="V"):
    """Load PDB, apply mutation, save."""
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.protocols.docking import setup_foldtree

    print(f"\nProcessing: {input_path}")

    pose = pose_from_pdb(str(input_path))
    print(f"  Loaded: {pose.size()} residues")

    # Setup foldtree for 2-chain system
    try:
        movable_jumps = pyrosetta.rosetta.utility.vector1_int()
        movable_jumps.append(1)
        setup_foldtree(pose, "A_B", movable_jumps)
    except Exception:
        pass  # Single chain or already set up

    # Apply mutation
    pose = mutate_and_repack(pose, chain_id, resnum, new_aa)

    # Verify
    pdb_info = pose.pdb_info()
    pose_idx = pdb_info.pdb2pose(chain_id, resnum)
    final_name = pose.residue(pose_idx).name3()
    assert final_name == "VAL", f"Mutation failed! Got {final_name}, expected VAL"

    # Save
    pose.dump_pdb(str(output_path))
    print(f"  Saved: {output_path}")
    print(f"  Verified: residue {chain_id}:{resnum} = {final_name}")

    return output_path


def main():
    init_pyrosetta()

    base = Path("input/PPI/prepared")
    targets = [
        ("EGFR_dimer_beta_meander.pdb", "EGFR_dimer_beta_meander_wt.pdb"),
        ("EGFR_dimer_TH1.pdb", "EGFR_dimer_TH1_wt.pdb"),
    ]

    print("=" * 60)
    print("  V924R → VAL Restoration (PDB resnum 948)")
    print("  Jura et al. 2009 Cell crystallization artifact")
    print("=" * 60)

    results = []
    for src, dst in targets:
        src_path = base / src
        dst_path = base / dst

        if not src_path.exists():
            print(f"\n[SKIP] {src_path} not found")
            continue

        if dst_path.exists():
            print(f"\n[SKIP] {dst_path} already exists (delete to regenerate)")
            continue

        out = process_pdb(src_path, dst_path)
        results.append(out)

    if results:
        print(f"\n{'=' * 60}")
        print(f"  Done! {len(results)} PDB(s) restored to wild-type.")
        print(f"  Original files preserved (not modified).")
        print(f"\n  Next steps:")
        print(f"  1. Update config INI files to point to *_wt.pdb")
        print(f"  2. Delete relaxed_cache/ (force re-relaxation with new PDB)")
        print(f"  3. Run PPI production docking")
        print(f"{'=' * 60}")
    else:
        print("\nNo files processed.")


if __name__ == "__main__":
    import pyrosetta
    main()
