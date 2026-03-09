# -*- coding: utf-8 -*-
"""
Pipeline Improvement Compatibility Check
=========================================
Tier 1~3 개선사항에 필요한 PyRosetta 함수/클래스 호환성 점검.
서버에서 실행: python check_improvements_compat.py
"""
import sys
import os
import time
import traceback

# ================================================================
#  PyRosetta 초기화 (배너 억제)
# ================================================================
original_stdout = sys.stdout
original_stderr = sys.stderr
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

print("=" * 70)
print("  Pipeline Improvement Compatibility Check")
print("=" * 70)

# ================================================================
#  Helper: 2-chain test pose
# ================================================================
print("\n[0] Preparing test pose...")
try:
    test_pose = pyrosetta.pose_from_sequence("ACDEFGHIKLMN/ACDEFGHIKLMN")
    from pyrosetta.rosetta.protocols.docking import setup_foldtree
    from pyrosetta.rosetta.utility import vector1_int
    movable = vector1_int()
    movable.append(1)
    setup_foldtree(test_pose, "A_B", movable)
    print(f"    [OK] Test pose: {test_pose.total_residue()} residues, "
          f"{test_pose.num_chains()} chains")
except Exception as e:
    print(f"    [FAIL] Cannot create test pose: {e}")
    sys.exit(1)

passed = 0
failed = 0
warnings = 0


def report(status, msg):
    global passed, failed, warnings
    if status == "OK":
        passed += 1
        print(f"    [OK]      {msg}")
    elif status == "WARN":
        warnings += 1
        print(f"    [WARNING] {msg}")
    else:
        failed += 1
        print(f"    [FAIL]    {msg}")


# ================================================================
#  TIER 1-1: Fast Scoring (dG + dSASA only, skip packstat/sc)
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 1-1] Fast Scoring: InterfaceAnalyzerMover minimal mode")
print("=" * 70)

try:
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta import create_score_function

    scorefxn = create_score_function("ref2015")

    # Test 1: Full scoring (current)
    iam_full = InterfaceAnalyzerMover()
    iam_full.set_interface("A_B")
    iam_full.set_scorefunction(scorefxn)
    iam_full.set_pack_separated(True)
    iam_full.set_compute_packstat(True)
    iam_full.set_compute_interface_sc(True)

    pose_full = test_pose.clone()
    scorefxn(pose_full)
    t0 = time.time()
    iam_full.apply(pose_full)
    t_full = time.time() - t0

    dG_full = iam_full.get_interface_dG()
    dSASA_full = iam_full.get_interface_delta_sasa()
    report("OK", f"Full scoring: dG={dG_full:.2f}, dSASA={dSASA_full:.2f} ({t_full*1000:.0f}ms)")

    # Test 2: Minimal scoring (skip packstat & sc)
    iam_fast = InterfaceAnalyzerMover()
    iam_fast.set_interface("A_B")
    iam_fast.set_scorefunction(scorefxn)
    iam_fast.set_pack_separated(True)
    iam_fast.set_compute_packstat(False)
    iam_fast.set_compute_interface_sc(False)

    pose_fast = test_pose.clone()
    scorefxn(pose_fast)
    t0 = time.time()
    iam_fast.apply(pose_fast)
    t_fast = time.time() - t0

    dG_fast = iam_fast.get_interface_dG()
    dSASA_fast = iam_fast.get_interface_delta_sasa()
    report("OK", f"Fast scoring: dG={dG_fast:.2f}, dSASA={dSASA_fast:.2f} ({t_fast*1000:.0f}ms)")

    # Compare
    if abs(dG_full - dG_fast) < 0.01 and abs(dSASA_full - dSASA_fast) < 0.01:
        report("OK", f"dG/dSASA identical between full and fast mode")
    else:
        report("WARN", f"dG/dSASA differ: full({dG_full:.4f},{dSASA_full:.4f}) "
               f"vs fast({dG_fast:.4f},{dSASA_fast:.4f})")

    speedup = t_full / t_fast if t_fast > 0 else 0
    report("OK", f"Speed: full={t_full*1000:.0f}ms, fast={t_fast*1000:.0f}ms "
           f"(speedup: {speedup:.1f}x)")

except Exception as e:
    report("FAIL", f"InterfaceAnalyzerMover fast mode: {e}")
    traceback.print_exc()

# ================================================================
#  TIER 1-2: DockMCMProtocol for refinement
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 1-2] DockMCMProtocol for local refinement")
print("=" * 70)

try:
    from pyrosetta.rosetta.protocols.docking import DockMCMProtocol

    # Test: create with default constructor
    dmc = DockMCMProtocol()
    report("OK", "DockMCMProtocol() default constructor")

    # Test: set_scorefxn
    dmc.set_scorefxn(create_score_function("ref2015"))
    report("OK", "set_scorefxn(ref2015)")

    # Check available methods for controlling refinement
    mcm_methods = [
        'set_scorefxn', 'set_scorefxn_pack',
        'set_move_since_last_cycle', 'set_filter_score',
        'set_partners', 'set_native_pose',
        'set_ignore_default_docking_task',
    ]
    for m in mcm_methods:
        if hasattr(dmc, m):
            report("OK", f"Method: {m}")
        else:
            report("WARN", f"Method not found: {m}")

    # Functional test: apply to pre-docked pose
    print("\n    --- Functional Test: DockMCMProtocol.apply() ---")
    from pyrosetta.rosetta.protocols.docking import DockingSlideIntoContact
    from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover

    ref_pose = test_pose.clone()

    # Small perturbation (simulating refinement, not global search)
    pert = RigidBodyPerturbMover(1, 1.0, 0.1)
    pert.apply(ref_pose)

    slide = DockingSlideIntoContact(1)
    slide.apply(ref_pose)

    scorefxn_ref = create_score_function("ref2015")
    score_before = scorefxn_ref(ref_pose)

    dmc_test = DockMCMProtocol()
    dmc_test.set_scorefxn(scorefxn_ref)

    t0 = time.time()
    dmc_test.apply(ref_pose)
    t_mcm = time.time() - t0

    score_after = scorefxn_ref(ref_pose)
    report("OK", f"DockMCMProtocol.apply() completed in {t_mcm:.2f}s")
    report("OK", f"Score: before={score_before:.2f} -> after={score_after:.2f} "
           f"(delta={score_after-score_before:.2f})")

except Exception as e:
    report("FAIL", f"DockMCMProtocol: {e}")
    traceback.print_exc()

# ================================================================
#  TIER 1-3: DockMCMProtocol with scorefxn_pack (separate packing sfxn)
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 1-3] DockMCMProtocol with separate packing scorefunction")
print("=" * 70)

try:
    dmc2 = DockMCMProtocol()
    sfxn_dock = create_score_function("ref2015")
    sfxn_pack = create_score_function("ref2015")

    dmc2.set_scorefxn(sfxn_dock)

    if hasattr(dmc2, 'set_scorefxn_pack'):
        dmc2.set_scorefxn_pack(sfxn_pack)
        report("OK", "Separate packing scorefunction supported")
    else:
        report("WARN", "set_scorefxn_pack not available (single sfxn will be used)")

except Exception as e:
    report("FAIL", f"Packing scorefunction test: {e}")

# ================================================================
#  TIER 1-4: Center-of-Mass distance pre-filter (pure Python)
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 1-4] Center-of-Mass pre-filter (numpy)")
print("=" * 70)

try:
    import numpy as np

    # Simulate center coordinates for 15000 models
    n_test = 15000
    centers = np.random.randn(n_test, 3) * 30  # spread across 30A
    threshold = 4.0

    # Test: pairwise distance check against a leader
    leader = centers[0]
    t0 = time.time()
    dists = np.sqrt(np.sum((centers - leader) ** 2, axis=1))
    within = np.sum(dists <= threshold)
    t_dist = time.time() - t0

    report("OK", f"CoM distance filter for {n_test} models: {t_dist*1000:.1f}ms "
           f"({within} within {threshold}A)")

    # Compare: 20 leaders vs 15000 candidates (batch)
    n_leaders = 20
    leaders_arr = centers[:n_leaders]
    t0 = time.time()
    # Broadcast: (15000, 1, 3) - (1, 20, 3) -> (15000, 20) distances
    all_dists = np.sqrt(np.sum(
        (centers[:, np.newaxis, :] - leaders_arr[np.newaxis, :, :]) ** 2,
        axis=2
    ))
    any_within = np.any(all_dists <= threshold, axis=1)
    t_batch = time.time() - t0
    report("OK", f"Batch CoM: {n_test} vs {n_leaders} leaders: {t_batch*1000:.1f}ms")

except Exception as e:
    report("FAIL", f"NumPy CoM pre-filter: {e}")

# ================================================================
#  TIER 2-1: Pareto ranking (pure Python)
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 2-1] Pareto ranking (numpy)")
print("=" * 70)

try:
    import numpy as np

    # Pareto front: no model dominated in ALL objectives
    def find_pareto_front(costs):
        """costs: (N, M) array where lower is better for all objectives."""
        is_pareto = np.ones(costs.shape[0], dtype=bool)
        for i in range(costs.shape[0]):
            if not is_pareto[i]:
                continue
            # A point dominates i if it's <= in all objectives and < in at least one
            dominated = np.all(costs <= costs[i], axis=1) & np.any(costs < costs[i], axis=1)
            is_pareto[dominated] = False
            is_pareto[i] = True  # keep self
        return is_pareto

    # Test with synthetic data (dG: lower better, -sc: lower better, -dSASA: lower better)
    n = 100
    test_data = np.column_stack([
        np.random.randn(n) * 5 - 15,  # dG
        -np.random.rand(n),            # -sc_value (lower = better sc)
        -np.random.rand(n) * 1000,     # -dSASA (lower = more buried)
    ])

    t0 = time.time()
    pareto = find_pareto_front(test_data)
    t_pareto = time.time() - t0
    report("OK", f"Pareto front: {np.sum(pareto)}/{n} models on front ({t_pareto*1000:.1f}ms)")

except Exception as e:
    report("FAIL", f"Pareto ranking: {e}")

# ================================================================
#  TIER 2-2: Boltzmann-weighted cluster scoring
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 2-2] Boltzmann-weighted scoring (numpy)")
print("=" * 70)

try:
    import numpy as np

    # Boltzmann probability: P(i) = exp(-dG_i / kT) / Z
    dGs = np.array([-23.3, -20.9, -20.6, -20.4, -19.4, -19.3, -19.2])
    clusters = ['C03', 'C03', 'C01', 'C03', 'C06', 'C06', 'C10']
    kT = 1.0  # Rosetta energy units

    # Shift to prevent overflow
    dGs_shifted = dGs - np.min(dGs)
    boltz = np.exp(-dGs_shifted / kT)
    Z = np.sum(boltz)
    probs = boltz / Z

    # Accumulate per cluster
    cluster_probs = {}
    for c, p in zip(clusters, probs):
        cluster_probs[c] = cluster_probs.get(c, 0) + p

    report("OK", "Boltzmann weighting computed successfully")
    for c, p in sorted(cluster_probs.items(), key=lambda x: -x[1]):
        print(f"              {c}: P={p:.4f} ({p*100:.1f}%)")

except Exception as e:
    report("FAIL", f"Boltzmann scoring: {e}")

# ================================================================
#  TIER 3-1: sc_value from getPoseExtraScore (for filtering)
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 3-1] sc_value accessibility for filtering")
print("=" * 70)

try:
    from pyrosetta.rosetta.core.pose import getPoseExtraScore

    # Apply IAM first to populate extra scores
    iam_sc = InterfaceAnalyzerMover()
    iam_sc.set_interface("A_B")
    iam_sc.set_scorefunction(create_score_function("ref2015"))
    iam_sc.set_pack_separated(True)
    iam_sc.set_compute_interface_sc(True)
    iam_sc.set_compute_packstat(False)

    pose_sc = test_pose.clone()
    create_score_function("ref2015")(pose_sc)
    iam_sc.apply(pose_sc)

    sc_val = getPoseExtraScore(pose_sc, "sc_value")
    report("OK", f"sc_value via getPoseExtraScore: {sc_val:.4f}")

    # Test without compute_interface_sc
    iam_no_sc = InterfaceAnalyzerMover()
    iam_no_sc.set_interface("A_B")
    iam_no_sc.set_scorefunction(create_score_function("ref2015"))
    iam_no_sc.set_pack_separated(True)
    iam_no_sc.set_compute_interface_sc(False)
    iam_no_sc.set_compute_packstat(False)

    pose_no_sc = test_pose.clone()
    create_score_function("ref2015")(pose_no_sc)
    iam_no_sc.apply(pose_no_sc)

    try:
        sc_val2 = getPoseExtraScore(pose_no_sc, "sc_value")
        report("WARN", f"sc_value available even without compute_interface_sc: {sc_val2:.4f}")
    except Exception:
        report("OK", "sc_value correctly unavailable when compute_interface_sc=False")

except Exception as e:
    report("FAIL", f"sc_value test: {e}")

# ================================================================
#  TIER 3-2: Memory estimation for 100K models
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 3-2] Memory estimation for 100K scale")
print("=" * 70)

try:
    import sys

    # Estimate PDB string size
    pdb_str = pyrosetta.rosetta.core.io.pdb.dump_pdb(test_pose)
    pdb_size = sys.getsizeof(pdb_str)
    residue_count = test_pose.total_residue()

    # Extrapolate for a typical 500-residue complex
    factor = 500 / residue_count if residue_count > 0 else 20
    est_pdb_size = pdb_size * factor

    est_100k_gb = (est_pdb_size * 100000) / (1024 ** 3)
    est_15k_gb = (est_pdb_size * 15000) / (1024 ** 3)

    report("OK", f"Test pose PDB string: {pdb_size:,} bytes ({residue_count} residues)")
    report("OK", f"Estimated for 500-res complex: ~{est_pdb_size/1024:.0f} KB per model")
    report("OK", f"100K models in memory: ~{est_100k_gb:.1f} GB")
    report("OK", f"15K filtered models: ~{est_15k_gb:.1f} GB")

    # Check available memory
    try:
        with open('/proc/meminfo', 'r') as f:
            for line in f:
                if line.startswith('MemTotal'):
                    total_mem_kb = int(line.split()[1])
                    total_mem_gb = total_mem_kb / (1024 ** 2)
                    report("OK", f"System total memory: {total_mem_gb:.1f} GB")
                    if total_mem_gb < est_100k_gb * 2:
                        report("WARN", f"100K may exceed memory "
                               f"(need ~{est_100k_gb*2:.0f} GB, have {total_mem_gb:.0f} GB)")
                    else:
                        report("OK", f"Memory sufficient for 100K models")
                    break
    except Exception:
        report("WARN", "Cannot read /proc/meminfo")

except Exception as e:
    report("FAIL", f"Memory estimation: {e}")

# ================================================================
#  TIER 3-3: Pose serialization speed benchmark
# ================================================================
print("\n" + "=" * 70)
print("  [Tier 3-3] Pose serialization/deserialization benchmark")
print("=" * 70)

try:
    import common

    pdb_string = common.pose_to_string(test_pose)

    # Benchmark string_to_pose
    n_iter = 10
    t0 = time.time()
    for _ in range(n_iter):
        p = common.string_to_pose(pdb_string)
    t_deser = (time.time() - t0) / n_iter

    report("OK", f"string_to_pose: {t_deser*1000:.1f}ms per call")
    report("OK", f"Estimated for 15K clustering candidates: {t_deser*15000:.0f}s "
           f"({t_deser*15000/60:.1f} min)")

    # Benchmark pose_to_string
    t0 = time.time()
    for _ in range(n_iter):
        s = common.pose_to_string(test_pose)
    t_ser = (time.time() - t0) / n_iter

    report("OK", f"pose_to_string: {t_ser*1000:.1f}ms per call")

except Exception as e:
    report("FAIL", f"Serialization benchmark: {e}")

# ================================================================
#  SUMMARY
# ================================================================
print("\n" + "=" * 70)
print(f"  SUMMARY: {passed} OK / {warnings} WARNING / {failed} FAIL")
print("=" * 70)

if failed == 0:
    print("\n  All critical checks passed. Improvements can proceed.")
else:
    print(f"\n  {failed} critical failure(s) detected. Review before proceeding.")

if warnings > 0:
    print(f"  {warnings} warning(s) - non-blocking but may need workarounds.")

print()
