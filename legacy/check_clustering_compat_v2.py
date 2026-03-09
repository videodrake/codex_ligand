# -*- coding: utf-8 -*-
"""
L_RMSD 클러스터링 호환성 테스트 V2
실제 PDB 파일 사용 (polyalanine 대신)
서버에서 실행: python check_clustering_compat_v2.py
"""
import sys, os, time, math, glob

# PyRosetta init
original_stdout, original_stderr = sys.stdout, sys.stderr
sys.stdout = open(os.devnull, 'w')
sys.stderr = open(os.devnull, 'w')
try:
    import pyrosetta
    pyrosetta.init("-mute all -multithreading:total_threads 1")
except RuntimeError:
    pass
finally:
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = original_stdout
    sys.stderr = original_stderr

from pyrosetta import pose_from_pdb, create_score_function
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover
from pyrosetta.rosetta.protocols.docking import DockingSlideIntoContact, setup_foldtree
from pyrosetta.rosetta.utility import vector1_int
import pyrosetta.rosetta.core.import_pose
import pyrosetta.rosetta.std
import numpy as np

TOP_DIR = os.path.dirname(os.path.abspath(__file__))

print("=" * 60)
print("Clustering Compatibility Test V2 (Real PDB)")
print("=" * 60)

# Find a real PDB file
pdb_files = glob.glob(os.path.join(TOP_DIR, "relaxed_cache", "*.pdb"))
if not pdb_files:
    pdb_files = glob.glob(os.path.join(TOP_DIR, "input_PDB", "*.pdb"))
if not pdb_files:
    print("[FAIL] No PDB files found in relaxed_cache/ or input_PDB/")
    sys.exit(1)

test_pdb = pdb_files[0]
print(f"\nUsing: {os.path.basename(test_pdb)}")

# Load reference pose
ref_pose = pose_from_pdb(test_pdb)
print(f"Loaded: {ref_pose.total_residue()} residues, {ref_pose.num_chains()} chains")

if ref_pose.num_chains() < 2:
    print("[FAIL] Need 2+ chains for docking test")
    sys.exit(1)

# Convert to PDB string (simulate pipeline's in-memory workflow)
s = pyrosetta.rosetta.std.ostringstream()
ref_pose.dump_pdb(s)
ref_pdb_string = s.str()


# --- Test 1: clone + calpha_superimpose on REAL structure ---
print("\n[1] calpha_superimpose_pose on real structure...")
try:
    mobile = ref_pose.clone()
    rmsd = calpha_superimpose_pose(mobile, ref_pose)
    print(f"    [OK] Self-RMSD: {rmsd:.6f}")
except Exception as e:
    # Self-superimpose might fail due to zero offset
    # Try with a slightly perturbed structure instead
    print(f"    [WARN] Self-superimpose failed (expected for zero RMSD): {e}")
    print("    Retrying with perturbed structure...")
    try:
        mobile = ref_pose.clone()
        dock_jump = 1
        movable_jumps = vector1_int()
        movable_jumps.append(dock_jump)
        setup_foldtree(mobile, "A_B", movable_jumps)
        pert = RigidBodyPerturbMover(dock_jump, 5.0, 2.0)  # small perturbation
        pert.apply(mobile)
        rmsd = calpha_superimpose_pose(mobile, ref_pose)
        print(f"    [OK] Perturbed RMSD: {rmsd:.6f}")
    except Exception as e2:
        print(f"    [FAIL] {e2}")
        sys.exit(1)


# --- Test 2: Generate 2 docked structures and compute L_RMSD between them ---
print("\n[2] L_RMSD between two different docked structures...")
try:
    dock_jump = 1
    movable_jumps = vector1_int()
    movable_jumps.append(dock_jump)

    # Docked structure A
    pose_a = ref_pose.clone()
    setup_foldtree(pose_a, "A_B", movable_jumps)
    RigidBodyPerturbMover(dock_jump, 45.0, 20.0).apply(pose_a)
    DockingSlideIntoContact(dock_jump).apply(pose_a)

    # Docked structure B
    pose_b = ref_pose.clone()
    setup_foldtree(pose_b, "A_B", movable_jumps)
    RigidBodyPerturbMover(dock_jump, 90.0, 30.0).apply(pose_b)
    DockingSlideIntoContact(dock_jump).apply(pose_b)

    # Compute L_RMSD: superimpose A onto B, then measure chain B deviation
    mobile = pose_a.clone()
    calpha_superimpose_pose(mobile, pose_b)

    chain_b_start = mobile.conformation().chain_begin(2)
    chain_b_end = mobile.conformation().chain_end(2)
    rmsd_sum = 0.0
    cnt = 0
    for i in range(chain_b_start, chain_b_end + 1):
        if mobile.residue(i).has("CA"):
            d2 = mobile.residue(i).xyz("CA").distance_squared(pose_b.residue(i).xyz("CA"))
            rmsd_sum += d2
            cnt += 1
    l_rmsd = math.sqrt(rmsd_sum / cnt) if cnt > 0 else 999.9
    print(f"    [OK] L_RMSD(A vs B): {l_rmsd:.2f} A ({cnt} CA atoms)")
except Exception as e:
    print(f"    [FAIL] {e}")


# --- Test 3: Performance benchmark with real structures ---
print("\n[3] Performance: clone + superimpose (real structure)...")
try:
    N = 100
    # Create a reference docked structure
    ref_docked = ref_pose.clone()
    setup_foldtree(ref_docked, "A_B", movable_jumps)
    RigidBodyPerturbMover(dock_jump, 30.0, 15.0).apply(ref_docked)
    DockingSlideIntoContact(dock_jump).apply(ref_docked)

    t0 = time.time()
    for _ in range(N):
        mob = ref_pose.clone()
        calpha_superimpose_pose(mob, ref_docked)
    elapsed = time.time() - t0

    per_op = elapsed / N * 1000
    print(f"    [OK] {N} ops in {elapsed:.2f}s ({per_op:.1f} ms/op)")

    # Estimates
    est_10k = per_op * 10000 / 1000 / 60
    est_30k = per_op * 30000 / 1000 / 60
    est_60k = per_op * 60000 / 1000 / 60
    print(f"    [EST] 10,000 ops: ~{est_10k:.1f} min")
    print(f"    [EST] 30,000 ops: ~{est_30k:.1f} min")
    print(f"    [EST] 60,000 ops: ~{est_60k:.1f} min")
except Exception as e:
    print(f"    [FAIL] {e}")


# --- Test 4: string_to_pose performance ---
print("\n[4] Performance: string_to_pose (real PDB string)...")
try:
    N = 30
    t0 = time.time()
    for _ in range(N):
        p = pyrosetta.Pose()
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(p, ref_pdb_string)
    elapsed = time.time() - t0
    per_op = elapsed / N * 1000
    print(f"    [OK] {per_op:.1f} ms/op ({ref_pose.total_residue()} residues)")

    est_3k = per_op * 3000 / 1000 / 60
    print(f"    [EST] 3,000 string_to_pose: ~{est_3k:.1f} min")
except Exception as e:
    print(f"    [FAIL] {e}")


# --- Test 5: Full clustering simulation ---
print("\n[5] Full clustering simulation (10 structures, greedy L_RMSD)...")
try:
    # Generate 10 docked structures
    structures = []
    for i in range(10):
        p = ref_pose.clone()
        setup_foldtree(p, "A_B", movable_jumps)
        RigidBodyPerturbMover(dock_jump, 360.0, 100.0).apply(p)
        DockingSlideIntoContact(dock_jump).apply(p)
        structures.append(p)
    print(f"    Generated 10 docked structures")

    # Greedy L_RMSD clustering
    THRESHOLD = 4.0
    leaders = [structures[0]]
    cluster_ids = [0]

    t0 = time.time()
    for i in range(1, len(structures)):
        cand = structures[i]
        assigned = False
        for lid, leader in enumerate(leaders):
            mob = cand.clone()
            calpha_superimpose_pose(mob, leader)
            chain_b_start = mob.conformation().chain_begin(2)
            chain_b_end = mob.conformation().chain_end(2)
            rmsd_sum, cnt = 0.0, 0
            for r in range(chain_b_start, chain_b_end + 1):
                if mob.residue(r).has("CA"):
                    d2 = mob.residue(r).xyz("CA").distance_squared(leader.residue(r).xyz("CA"))
                    rmsd_sum += d2
                    cnt += 1
            l_rmsd = math.sqrt(rmsd_sum / cnt) if cnt > 0 else 999.9
            if l_rmsd <= THRESHOLD:
                cluster_ids.append(lid)
                assigned = True
                break
        if not assigned:
            cluster_ids.append(len(leaders))
            leaders.append(cand)
    elapsed = time.time() - t0

    print(f"    [OK] {len(leaders)} clusters found in {elapsed:.2f}s")
    for i, lid in enumerate(cluster_ids):
        print(f"         Structure {i} → Cluster {lid}")
except Exception as e:
    print(f"    [FAIL] {e}")


print("\n" + "=" * 60)
print("All checks complete.")
print("=" * 60)
