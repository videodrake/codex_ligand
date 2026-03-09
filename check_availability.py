# -*- coding: utf-8 -*-
"""
PyRosetta 도킹 관련 클래스/함수 가용성 확인 스크립트
서버에서 실행: python check_availability.py
"""
import sys
import os

# PyRosetta 초기화 (배너 억제)
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

print("=" * 60)
print("PyRosetta Availability Check")
print("=" * 60)

# 1. PyRosetta 버전 정보
print("\n[1] PyRosetta Version Info")
try:
    print(f"    Version: {pyrosetta.rosetta.utility.Version.version()}")
except:
    pass
try:
    print(f"    Commit:  {pyrosetta.rosetta.utility.Version.commit()}")
except:
    pass
try:
    print(f"    Date:    {pyrosetta.rosetta.utility.Version.date()}")
except:
    pass

# 2. 글로벌 도킹 핵심 클래스
print("\n[2] Global Docking - Slide Into Contact")
checks = {
    "DockingSlideIntoContact": "pyrosetta.rosetta.protocols.docking.DockingSlideIntoContact",
    "FaDock": "pyrosetta.rosetta.protocols.docking.FaDock",
    "DockingProtocol": "pyrosetta.rosetta.protocols.docking.DockingProtocol",
    "DockingLowRes": "pyrosetta.rosetta.protocols.docking.DockingLowRes",
    "DockingHighRes": "pyrosetta.rosetta.protocols.docking.DockingHighRes",
    "DockMCMProtocol": "pyrosetta.rosetta.protocols.docking.DockMCMProtocol",
    "setup_foldtree": "pyrosetta.rosetta.protocols.docking.setup_foldtree",
}

for name, path in checks.items():
    try:
        parts = path.rsplit('.', 1)
        module = __import__(parts[0], fromlist=[parts[1]])
        obj = getattr(module, parts[1])
        print(f"    [OK] {name}")
        # 생성자 시그니처 확인
        try:
            import inspect
            sig = inspect.signature(obj.__init__)
            print(f"         Signature: {name}{sig}")
        except:
            pass
    except (ImportError, AttributeError) as e:
        print(f"    [MISSING] {name} - {e}")

# 3. RigidBody 관련
print("\n[3] RigidBody Movers")
rb_checks = {
    "RigidBodyPerturbMover": "pyrosetta.rosetta.protocols.rigid.RigidBodyPerturbMover",
    "RigidBodyRandomizeMover": "pyrosetta.rosetta.protocols.rigid.RigidBodyRandomizeMover",
    "RigidBodyTransMover": "pyrosetta.rosetta.protocols.rigid.RigidBodyTransMover",
    "RigidBodySpinMover": "pyrosetta.rosetta.protocols.rigid.RigidBodySpinMover",
}

for name, path in rb_checks.items():
    try:
        parts = path.rsplit('.', 1)
        module = __import__(parts[0], fromlist=[parts[1]])
        obj = getattr(module, parts[1])
        print(f"    [OK] {name}")
    except (ImportError, AttributeError) as e:
        print(f"    [MISSING] {name} - {e}")

# 4. InterfaceAnalyzerMover 메서드
print("\n[4] InterfaceAnalyzerMover Methods")
try:
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    iam = InterfaceAnalyzerMover()
    methods_to_check = [
        'set_interface', 'set_pack_separated', 'set_compute_packstat',
        'set_scorefunction', 'set_compute_interface_sc',
        'get_interface_dG', 'get_interface_delta_sasa',
        'get_sc_value', 'get_packstat',
        'get_separated_interface_energy', 'get_complexed_sasa',
        'get_interface_delta_hbond_unsat',
    ]
    for m in methods_to_check:
        if hasattr(iam, m):
            print(f"    [OK] {m}")
        else:
            print(f"    [MISSING] {m}")
except Exception as e:
    print(f"    [ERROR] Cannot create InterfaceAnalyzerMover: {e}")

# 5. ScoreFunction 관련
print("\n[5] Score Functions")
for sfxn_name in ["ref2015", "interchain_cen", "docking", "docking_low_res"]:
    try:
        sfxn = pyrosetta.create_score_function(sfxn_name)
        print(f"    [OK] {sfxn_name}")
    except Exception as e:
        print(f"    [MISSING] {sfxn_name} - {e}")

# 6. 실제 SlideIntoContact 동작 테스트
print("\n[6] Functional Test - SlideIntoContact")
try:
    from pyrosetta.rosetta.protocols.docking import DockingSlideIntoContact
    # 간단한 2-chain 테스트 포즈 생성
    test_pose = pyrosetta.pose_from_sequence("AAAAA/AAAAA")
    if test_pose.num_chains() >= 2:
        slide = DockingSlideIntoContact(1)
        slide.apply(test_pose)
        print(f"    [OK] DockingSlideIntoContact works (test pose: {test_pose.num_chains()} chains)")
    else:
        print(f"    [WARN] Test pose has only {test_pose.num_chains()} chain(s)")
except Exception as e:
    print(f"    [FAIL] {e}")

# 7. DockingProtocol 동작 테스트
print("\n[7] Functional Test - DockingProtocol")
try:
    from pyrosetta.rosetta.protocols.docking import DockingProtocol
    # DockingProtocol(jump, low_res, high_res)
    dp = DockingProtocol(1, False, True)
    print(f"    [OK] DockingProtocol(1, False, True) created")
except Exception as e:
    print(f"    [FAIL] {e}")

try:
    from pyrosetta.rosetta.protocols.docking import DockingProtocol
    dp = DockingProtocol(1, True, True)
    print(f"    [OK] DockingProtocol(1, True, True) created")
except Exception as e:
    print(f"    [FAIL] {e}")

# 8. Vector1 유틸리티
print("\n[8] Utility Classes")
try:
    from pyrosetta.rosetta.utility import vector1_int
    v = vector1_int()
    v.append(1)
    print(f"    [OK] vector1_int")
except (ImportError, AttributeError) as e:
    print(f"    [MISSING] vector1_int - {e}")

try:
    from pyrosetta.rosetta.std import vector_unsigned_long
    print(f"    [OK] vector_unsigned_long")
except (ImportError, AttributeError) as e:
    print(f"    [MISSING] vector_unsigned_long - {e}")

print("\n" + "=" * 60)
print("Check Complete.")
print("=" * 60)
