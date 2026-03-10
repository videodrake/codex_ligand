#!/usr/bin/env python
"""Restore V924R mutation to wild-type VAL in PPI dimer PDBs.

3GT8 crystal structure contains V924R (PDB resnum 948) engineering mutation
(Jura et al. 2009 Cell) for crystallization. This must be reverted to
wild-type VAL before PPI docking, since MYO1D binds wild-type EGFR.

Usage:
    python scripts/fix_v924r_mutation.py              # PDB text method (no PyRosetta)
    python scripts/fix_v924r_mutation.py --pyrosetta   # PyRosetta method (repack)

This script:
  1. Loads each dimer PDB
  2. Replaces ARG 948 chain A with VAL (keeping backbone + CB)
  3. Saves to *_wt.pdb (original files preserved)
  4. Verifies the mutation was applied
"""
import sys
import os
from pathlib import Path


# ── ARG sidechain atoms beyond CB (these get removed for VAL) ─────────
_ARG_SIDECHAIN_PAST_CB = {"CG", "CD", "NE", "CZ", "NH1", "NH2",
                           "HG2", "HG3", "HD2", "HD3", "HE",
                           "HH11", "HH12", "HH21", "HH22",
                           "1HG", "2HG", "1HD", "2HD", "1HH1",
                           "2HH1", "1HH2", "2HH2"}

# VAL sidechain atoms beyond CB that we need to add: CG1, CG2
# We place them using ideal geometry relative to CA-CB bond.
# For PDB text method, we skip adding CG1/CG2 — PyRosetta's Relax
# step (first step of the PPI pipeline) will rebuild them automatically.


def _pdb_resnum_chain(line: str):
    """Extract (chain, resnum) from an ATOM/HETATM PDB line."""
    chain = line[21]
    try:
        resnum = int(line[22:26].strip())
    except ValueError:
        return None, None
    return chain, resnum


def _pdb_resname(line: str) -> str:
    return line[17:20].strip()


def _pdb_atom_name(line: str) -> str:
    return line[12:16].strip()


def mutate_arg_to_val_pdb_text(input_path: str, output_path: str,
                                chain_id: str = "A", resnum: int = 948) -> bool:
    """Pure PDB text manipulation: ARG→VAL at specified position.

    Strategy:
      - Keep all backbone atoms (N, CA, C, O, H, HA) and CB, HB
      - Remove ARG-specific sidechain atoms (CG, CD, NE, CZ, NH1, NH2 + H's)
      - Rename residue from ARG to VAL
      - PyRosetta's FastRelax (pipeline step 1) will add CG1/CG2 automatically

    Returns True if mutation applied, False if residue was not ARG.
    """
    lines = []
    found = False
    is_arg = False
    removed_count = 0

    with open(input_path, "r") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                lines.append(line)
                continue

            ch, rn = _pdb_resnum_chain(line)
            if ch == chain_id and rn == resnum:
                found = True
                resname = _pdb_resname(line)
                atom_name = _pdb_atom_name(line)

                if resname != "ARG":
                    if resname == "VAL":
                        print(f"  Residue {chain_id}:{resnum} is already VAL — no change needed")
                        return False
                    print(f"  WARNING: Residue {chain_id}:{resnum} is {resname}, not ARG")
                    lines.append(line)
                    continue

                is_arg = True

                # Remove sidechain atoms past CB
                if atom_name in _ARG_SIDECHAIN_PAST_CB:
                    removed_count += 1
                    continue  # skip this atom line

                # Rename ARG → VAL in the line
                line = line[:17] + "VAL" + line[20:]
                lines.append(line)
            else:
                lines.append(line)

    if not found:
        raise ValueError(f"Residue {chain_id}:{resnum} not found in {input_path}")

    if not is_arg:
        return False

    with open(output_path, "w") as f:
        f.writelines(lines)

    print(f"  Mutated: ARG → VAL (removed {removed_count} sidechain atoms)")
    print(f"  Note: CG1/CG2 will be built by FastRelax (pipeline step 1)")
    return True


def verify_pdb(output_path: str, chain_id: str = "A", resnum: int = 948) -> bool:
    """Verify the output PDB has VAL at the target position."""
    with open(output_path, "r") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            ch, rn = _pdb_resnum_chain(line)
            if ch == chain_id and rn == resnum:
                resname = _pdb_resname(line)
                if resname == "VAL":
                    return True
                else:
                    print(f"  ERROR: Expected VAL at {chain_id}:{resnum}, got {resname}")
                    return False
    print(f"  ERROR: Residue {chain_id}:{resnum} not found in output")
    return False


# ── PyRosetta method (fallback, better sidechain placement) ───────────

def process_pdb_pyrosetta(input_path, output_path, chain_id="A", resnum=948):
    """PyRosetta method: mutate + repack. May segfault on some versions."""
    import pyrosetta
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
    from pyrosetta.rosetta.core.select.residue_selector import (
        NeighborhoodResidueSelector,
        ResidueIndexSelector,
        NotResidueSelector,
    )
    from pyrosetta.rosetta.core.pack.task.operation import (
        OperateOnResidueSubset,
        PreventRepackingRLT,
    )

    pose = pose_from_pdb(str(input_path))
    print(f"  Loaded: {pose.size()} residues")

    pdb_info = pose.pdb_info()
    pose_idx = pdb_info.pdb2pose(chain_id, resnum)
    if pose_idx == 0:
        raise ValueError(f"Residue {chain_id}:{resnum} not found")

    old_name = pose.residue(pose_idx).name3()
    print(f"  Residue {chain_id}:{resnum} = {old_name} (pose index {pose_idx})")

    # Mutate (no foldtree setup — not needed for single-residue mutation)
    mutator = MutateResidue(pose_idx, "V")
    mutator.apply(pose)
    print(f"  Mutated: {old_name} → {pose.residue(pose_idx).name3()}")

    # Repack neighborhood
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    mut_selector = ResidueIndexSelector(str(pose_idx))
    nbr_selector = NeighborhoodResidueSelector(mut_selector, 8.0, True)
    not_nbr = NotResidueSelector(nbr_selector)

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), not_nbr))

    packer = PackRotamersMover(sfxn)
    packer.task_factory(tf)
    packer.apply(pose)
    print(f"  Repacked 8Å neighborhood")

    # Verify & save
    final_name = pose.residue(pose_idx).name3()
    assert final_name == "VAL", f"Mutation failed: got {final_name}"
    pose.dump_pdb(str(output_path))
    print(f"  Saved: {output_path}")


def main():
    use_pyrosetta = "--pyrosetta" in sys.argv

    base = Path("input/PPI/prepared")
    targets = [
        ("EGFR_dimer_beta_meander.pdb", "EGFR_dimer_beta_meander_wt.pdb"),
        ("EGFR_dimer_TH1.pdb", "EGFR_dimer_TH1_wt.pdb"),
    ]

    method = "PyRosetta" if use_pyrosetta else "PDB text"
    print("=" * 60)
    print("  V924R → VAL Restoration (PDB resnum 948)")
    print("  Jura et al. 2009 Cell crystallization artifact")
    print(f"  Method: {method}")
    print("=" * 60)

    if use_pyrosetta:
        import pyrosetta
        pyrosetta.init(
            "-mute all -ignore_unrecognized_res true",
            set_logging_handler=None,
        )

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

        print(f"\nProcessing: {src_path}")

        try:
            if use_pyrosetta:
                process_pdb_pyrosetta(src_path, dst_path)
            else:
                ok = mutate_arg_to_val_pdb_text(str(src_path), str(dst_path))
                if not ok:
                    continue

            # Verify
            if verify_pdb(str(dst_path)):
                print(f"  Verified: {dst_path.name} has VAL at A:948")
                results.append(dst_path)
            else:
                print(f"  FAILED verification — removing {dst_path}")
                dst_path.unlink(missing_ok=True)
        except Exception as e:
            print(f"  ERROR: {e}")
            dst_path.unlink(missing_ok=True)

    if results:
        print(f"\n{'=' * 60}")
        print(f"  Done! {len(results)} PDB(s) restored to wild-type.")
        print(f"  Original files preserved (not modified).")
        if not use_pyrosetta:
            print(f"\n  VAL sidechain (CG1/CG2) will be built automatically")
            print(f"  by FastRelax (pipeline step 1).")
        print(f"\n  Next: qsub -v RUN_MODE=both config/run_ppi_test.pbs")
        print(f"{'=' * 60}")
    else:
        print("\nNo files processed.")


if __name__ == "__main__":
    main()
