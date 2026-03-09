# -*- coding: utf-8 -*-
"""
Filter v2.0 Compatibility Check
================================
v2.0 필터링 파이프라인에 필요한 새 메트릭/함수 호환성 점검.
서버에서 실행: python check_filter_v2_compat.py
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
print("  Filter v2.0 Compatibility Check")
print("=" * 70)

# ================================================================
#  Helper: 2-chain test pose
# ================================================================
print("\n[0] Preparing test pose...")
try:
    test_pose = pyrosetta.pose_from_sequence("ACDEFGHIKLMN/ACDEFGHIKLMN")
    from pyrosetta.rosetta.protocols.docking import setup_foldtree, DockingSlideIntoContact
    from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover
    from pyrosetta.rosetta.utility import vector1_int
    movable = vector1_int()
    movable.append(1)
    setup_foldtree(test_pose, "A_B", movable)

    # Separate chains to avoid zero-length vector errors in IAM
    pert = RigidBodyPerturbMover(1, 8.0, 30.0)
    pert.apply(test_pose)
    slide = DockingSlideIntoContact(1)
    slide.apply(test_pose)

    from pyrosetta import create_score_function as _csf
    _csf("ref2015")(test_pose)

    print(f"    [OK] Test pose: {test_pose.total_residue()} residues, "
          f"{test_pose.num_chains()} chains (chains separated)")
except Exception as e:
    print(f"    [FAIL] Cannot create test pose: {e}")
    traceback.print_exc()
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
#  Prepare scored pose for metric tests
# ================================================================
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta import create_score_function
from pyrosetta.rosetta.core.pose import getPoseExtraScore

scorefxn = create_score_function("ref2015")

# ================================================================
#  Test 1: delta_unsatHbonds (get_interface_delta_hbond_unsat)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 1] delta_unsatHbonds: iam.get_interface_delta_hbond_unsat()")
print("=" * 70)

try:
    iam = InterfaceAnalyzerMover()
    iam.set_interface("A_B")
    iam.set_scorefunction(scorefxn)
    iam.set_pack_separated(True)
    iam.set_compute_packstat(False)
    iam.set_compute_interface_sc(False)

    pose1 = test_pose.clone()
    scorefxn(pose1)
    iam.apply(pose1)

    if hasattr(iam, 'get_interface_delta_hbond_unsat'):
        val = iam.get_interface_delta_hbond_unsat()
        report("OK", f"get_interface_delta_hbond_unsat() = {val}")
    else:
        report("WARN", "get_interface_delta_hbond_unsat() not available "
               "(will use fallback default=99)")

    # Also check getPoseExtraScore key
    try:
        val2 = getPoseExtraScore(pose1, "delta_unsatHbonds")
        report("OK", f"getPoseExtraScore 'delta_unsatHbonds' = {val2}")
    except Exception:
        report("WARN", "getPoseExtraScore 'delta_unsatHbonds' not available "
               "(will use method call)")

except Exception as e:
    report("FAIL", f"delta_unsatHbonds test: {e}")
    traceback.print_exc()

# ================================================================
#  Test 2: nres_int (get_num_interface_residues)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 2] nres_int: iam.get_num_interface_residues()")
print("=" * 70)

try:
    iam2 = InterfaceAnalyzerMover()
    iam2.set_interface("A_B")
    iam2.set_scorefunction(scorefxn)
    iam2.set_pack_separated(True)
    iam2.set_compute_packstat(False)
    iam2.set_compute_interface_sc(False)

    pose2 = test_pose.clone()
    scorefxn(pose2)
    iam2.apply(pose2)

    if hasattr(iam2, 'get_num_interface_residues'):
        val = iam2.get_num_interface_residues()
        report("OK", f"get_num_interface_residues() = {val}")
    else:
        report("WARN", "get_num_interface_residues() not available "
               "(will use fallback default=0)")

    # Also check getPoseExtraScore key
    try:
        val2 = getPoseExtraScore(pose2, "nres_int")
        report("OK", f"getPoseExtraScore 'nres_int' = {val2}")
    except Exception:
        report("WARN", "getPoseExtraScore 'nres_int' not available "
               "(will use method call)")

except Exception as e:
    report("FAIL", f"nres_int test: {e}")
    traceback.print_exc()

# ================================================================
#  Test 3: packstat method name (get_packstat vs get_interface_packstat)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 3] packstat method name variants")
print("=" * 70)

try:
    iam3 = InterfaceAnalyzerMover()
    iam3.set_interface("A_B")
    iam3.set_scorefunction(scorefxn)
    iam3.set_pack_separated(True)
    iam3.set_compute_packstat(True)
    iam3.set_compute_interface_sc(False)

    pose3 = test_pose.clone()
    scorefxn(pose3)
    iam3.apply(pose3)

    # Check getPoseExtraScore
    try:
        ps_extra = getPoseExtraScore(pose3, "packstat")
        report("OK", f"getPoseExtraScore 'packstat' = {ps_extra:.4f}")
    except Exception:
        report("WARN", "getPoseExtraScore 'packstat' not available")

    # Check get_packstat
    if hasattr(iam3, 'get_packstat'):
        try:
            ps_val = iam3.get_packstat()
            report("OK", f"iam.get_packstat() = {ps_val:.4f}")
        except Exception as e:
            report("WARN", f"iam.get_packstat() exists but failed: {e}")
    else:
        report("WARN", "iam.get_packstat() not available")

    # Check get_interface_packstat
    if hasattr(iam3, 'get_interface_packstat'):
        try:
            ps_val2 = iam3.get_interface_packstat()
            report("OK", f"iam.get_interface_packstat() = {ps_val2:.4f}")
        except Exception as e:
            report("WARN", f"iam.get_interface_packstat() exists but failed: {e}")
    else:
        report("WARN", "iam.get_interface_packstat() not available")

except Exception as e:
    report("FAIL", f"packstat method test: {e}")
    traceback.print_exc()

# ================================================================
#  Test 3.5: hbonds_int (interface hydrogen bonds)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 3.5] hbonds_int: interface hydrogen bond count")
print("=" * 70)

try:
    iam3h = InterfaceAnalyzerMover()
    iam3h.set_interface("A_B")
    iam3h.set_scorefunction(scorefxn)
    iam3h.set_pack_separated(True)
    iam3h.set_compute_packstat(False)
    iam3h.set_compute_interface_sc(False)

    pose3h = test_pose.clone()
    scorefxn(pose3h)
    iam3h.apply(pose3h)

    # Check getPoseExtraScore
    try:
        hb_val = getPoseExtraScore(pose3h, "hbonds_int")
        report("OK", f"getPoseExtraScore 'hbonds_int' = {hb_val:.4f}")
    except Exception:
        report("WARN", "getPoseExtraScore 'hbonds_int' not available")

    # Check method
    if hasattr(iam3h, 'get_interface_hbonds'):
        try:
            hb_val2 = iam3h.get_interface_hbonds()
            report("OK", f"iam.get_interface_hbonds() = {hb_val2}")
        except Exception as e:
            report("WARN", f"iam.get_interface_hbonds() exists but failed: {e}")
    else:
        report("WARN", "iam.get_interface_hbonds() not available")

except Exception as e:
    report("FAIL", f"hbonds_int test: {e}")
    traceback.print_exc()

# ================================================================
#  Test 4: total_score capture from scorefxn(pose)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 4] total_score = scorefxn(pose) return value")
print("=" * 70)

try:
    pose4 = test_pose.clone()
    total_score = scorefxn(pose4)
    if isinstance(total_score, float):
        report("OK", f"scorefxn(pose) returns float: {total_score:.4f}")
    else:
        report("WARN", f"scorefxn(pose) returns {type(total_score)}: {total_score}")
except Exception as e:
    report("FAIL", f"total_score capture: {e}")

# ================================================================
#  Test 5: dG_density derived metric
# ================================================================
print("\n" + "=" * 70)
print("  [Test 5] dG_density = dG_separated / dSASA * 100")
print("=" * 70)

try:
    iam5 = InterfaceAnalyzerMover()
    iam5.set_interface("A_B")
    iam5.set_scorefunction(scorefxn)
    iam5.set_pack_separated(True)
    iam5.set_compute_packstat(False)
    iam5.set_compute_interface_sc(False)

    pose5 = test_pose.clone()
    scorefxn(pose5)
    iam5.apply(pose5)

    dG = iam5.get_interface_dG()
    dSASA = iam5.get_interface_delta_sasa()

    if abs(dSASA) > 1e-6:
        dG_density = dG / dSASA * 100
        report("OK", f"dG={dG:.2f}, dSASA={dSASA:.2f}, "
               f"dG_density={dG_density:.4f}")
    else:
        report("WARN", f"dSASA ~ 0 ({dSASA:.6f}), dG_density undefined "
               "(will use safe division)")

except Exception as e:
    report("FAIL", f"dG_density calculation: {e}")

# ================================================================
#  Test 6: packstat performance benchmark (fast vs expensive)
# ================================================================
print("\n" + "=" * 70)
print("  [Test 6] IAM performance: packstat ON vs OFF")
print("=" * 70)

try:
    # Without packstat
    iam_fast = InterfaceAnalyzerMover()
    iam_fast.set_interface("A_B")
    iam_fast.set_scorefunction(scorefxn)
    iam_fast.set_pack_separated(True)
    iam_fast.set_compute_packstat(False)
    iam_fast.set_compute_interface_sc(True)

    pose_fast = test_pose.clone()
    scorefxn(pose_fast)
    t0 = time.time()
    iam_fast.apply(pose_fast)
    t_fast = time.time() - t0

    # With packstat
    iam_exp = InterfaceAnalyzerMover()
    iam_exp.set_interface("A_B")
    iam_exp.set_scorefunction(scorefxn)
    iam_exp.set_pack_separated(True)
    iam_exp.set_compute_packstat(True)
    iam_exp.set_compute_interface_sc(False)

    pose_exp = test_pose.clone()
    scorefxn(pose_exp)
    t0 = time.time()
    iam_exp.apply(pose_exp)
    t_exp = time.time() - t0

    overhead = t_exp / t_fast if t_fast > 0 else float('inf')
    report("OK", f"Without packstat: {t_fast*1000:.1f}ms")
    report("OK", f"With packstat:    {t_exp*1000:.1f}ms "
           f"(overhead: {overhead:.1f}x)")

    if overhead > 5:
        report("WARN", f"Packstat is {overhead:.1f}x slower - "
               "consider only for Stage 2 survivors")
    else:
        report("OK", f"Packstat overhead acceptable ({overhead:.1f}x)")

except Exception as e:
    report("FAIL", f"Packstat benchmark: {e}")

# ================================================================
#  Test 7: getPoseExtraScore keys for new metrics
# ================================================================
print("\n" + "=" * 70)
print("  [Test 7] getPoseExtraScore keys for v2.0 metrics")
print("=" * 70)

try:
    # Full IAM with everything enabled
    iam7 = InterfaceAnalyzerMover()
    iam7.set_interface("A_B")
    iam7.set_scorefunction(scorefxn)
    iam7.set_pack_separated(True)
    iam7.set_compute_packstat(True)
    iam7.set_compute_interface_sc(True)

    pose7 = test_pose.clone()
    scorefxn(pose7)
    iam7.apply(pose7)

    keys_to_check = [
        "dG_separated", "dSASA_int", "delta_unsatHbonds",
        "packstat", "sc_value", "nres_int",
        "dG_separated/dSASAx100", "hbonds_int",
        "per_residue_energy_int", "side1_score", "side2_score",
    ]

    for key in keys_to_check:
        try:
            val = getPoseExtraScore(pose7, key)
            report("OK", f"'{key}' = {val:.4f}")
        except Exception:
            report("WARN", f"'{key}' not available in getPoseExtraScore")

except Exception as e:
    report("FAIL", f"getPoseExtraScore key test: {e}")

# ================================================================
#  Test 8: sklearn availability
# ================================================================
print("\n" + "=" * 70)
print("  [Test 8] sklearn availability")
print("=" * 70)

try:
    import sklearn
    report("OK", f"sklearn version: {sklearn.__version__}")
except ImportError:
    report("WARN", "sklearn not available (expected on network-isolated HPC)")

# ================================================================
#  Test 9: scipy availability
# ================================================================
print("\n" + "=" * 70)
print("  [Test 9] scipy availability")
print("=" * 70)

try:
    import scipy
    report("OK", f"scipy version: {scipy.__version__}")
    from scipy.spatial.distance import cdist
    report("OK", "scipy.spatial.distance.cdist available")
except ImportError:
    report("WARN", "scipy not available (clustering will use pure numpy)")

# ================================================================
#  SUMMARY
# ================================================================
print("\n" + "=" * 70)
print(f"  SUMMARY: {passed} OK / {warnings} WARNING / {failed} FAIL")
print("=" * 70)

if failed == 0:
    print("\n  All critical checks passed. v2.0 filtering can proceed.")
    print("  WARNINGs indicate optional features that will use fallback defaults.")
else:
    print(f"\n  {failed} critical failure(s) detected. Review before proceeding.")

if warnings > 0:
    print(f"  {warnings} warning(s) - non-blocking, fallback strategies available.")

print()
