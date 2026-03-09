# -*- coding: utf-8 -*-
"""
L_RMSD 기반 클러스터링 전환을 위한 호환성 + 성능 테스트
서버에서 실행: python check_clustering_compat.py
"""
import sys, os, time

# PyRosetta init (suppress banner)
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

from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
import numpy as np

print("=" * 60)
print("Clustering Compatibility & Performance Test")
print("=" * 60)

# 1. Create a 2-chain test pose
print("\n[1] Creating 2-chain test pose...")
try:
    pose = pyrosetta.pose_from_sequence("A" * 50 + "/" + "A" * 30)
    print(f"    [OK] Test pose: {pose.total_residue()} residues, {pose.num_chains()} chains")
except Exception as e:
    print(f"    [FAIL] {e}")
    sys.exit(1)

# 2. pose.clone()
print("\n[2] Testing pose.clone()...")
try:
    clone = pose.clone()
    print(f"    [OK] pose.clone() works ({clone.total_residue()} residues)")
except Exception as e:
    print(f"    [FAIL] {e}")

# 3. calpha_superimpose_pose
print("\n[3] Testing calpha_superimpose_pose...")
try:
    mobile = pose.clone()
    ref = pose.clone()
    rmsd = calpha_superimpose_pose(mobile, ref)
    print(f"    [OK] calpha_superimpose_pose works (self-RMSD: {rmsd:.4f})")
except Exception as e:
    print(f"    [FAIL] {e}")

# 4. Performance: clone + superimpose speed
print("\n[4] Performance benchmark (clone + superimpose)...")
try:
    ref_pose = pose.clone()
    N_ITERS = 100

    t0 = time.time()
    for _ in range(N_ITERS):
        mob = pose.clone()
        calpha_superimpose_pose(mob, ref_pose)
    elapsed = time.time() - t0

    per_op = elapsed / N_ITERS * 1000  # ms
    print(f"    [OK] {N_ITERS} iterations in {elapsed:.2f}s ({per_op:.1f} ms/op)")

    # Estimate for real workload:
    # 3000 candidates x avg 10 leader comparisons = 30,000 ops
    est_30k = per_op * 30000 / 1000 / 60  # minutes
    print(f"    [EST] ~30,000 ops (3000 cands x 10 leaders): ~{est_30k:.1f} min")
except Exception as e:
    print(f"    [FAIL] {e}")

# 5. string_to_pose + clone speed (realistic scenario)
print("\n[5] Performance: string_to_pose + clone...")
try:
    import pyrosetta.rosetta.std
    import pyrosetta.rosetta.core.import_pose

    # Generate PDB string
    s = pyrosetta.rosetta.std.ostringstream()
    pose.dump_pdb(s)
    pdb_string = s.str()

    N_ITERS = 50
    t0 = time.time()
    for _ in range(N_ITERS):
        p = pyrosetta.Pose()
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(p, pdb_string)
    elapsed = time.time() - t0

    per_op = elapsed / N_ITERS * 1000  # ms
    print(f"    [OK] string_to_pose: {per_op:.1f} ms/op")

    # For clustering: create each candidate pose once, then clone for each leader
    # 3000 string_to_pose + 30000 clones
    est_total = (per_op * 3000 / 1000 / 60)
    print(f"    [EST] 3000 string_to_pose calls: ~{est_total:.1f} min")
except Exception as e:
    print(f"    [FAIL] {e}")

# 6. numpy percentile
print("\n[6] Testing numpy.percentile...")
try:
    data = np.random.randn(10000)
    p10 = np.percentile(data, 10)
    print(f"    [OK] np.percentile works (10th percentile of random data: {p10:.2f})")
except Exception as e:
    print(f"    [FAIL] {e}")

# 7. L_RMSD computation (chain B only)
print("\n[7] Testing L_RMSD computation (chain B)...")
try:
    import math
    mobile = pose.clone()
    ref = pose.clone()
    calpha_superimpose_pose(mobile, ref)

    chain_b_start = mobile.conformation().chain_begin(2)
    chain_b_end = mobile.conformation().chain_end(2)

    rmsd_sum = 0.0
    cnt = 0
    for i in range(chain_b_start, chain_b_end + 1):
        if mobile.residue(i).has("CA"):
            d2 = mobile.residue(i).xyz("CA").distance_squared(ref.residue(i).xyz("CA"))
            rmsd_sum += d2
            cnt += 1

    l_rmsd = math.sqrt(rmsd_sum / cnt) if cnt > 0 else 999.9
    print(f"    [OK] L_RMSD computation works (self-RMSD: {l_rmsd:.4f}, {cnt} CA atoms)")
except Exception as e:
    print(f"    [FAIL] {e}")

print("\n" + "=" * 60)
print("All checks complete.")
print("=" * 60)
