#!/usr/bin/env python
"""구문검증 스크립트 - 서버에서 실행: python check_syntax.py"""
import ast
import sys

files = [
    "common.py",
    "docking.py",
    "analysis.py",
    "pipeline_manager.py",
]

errors = 0
for f in files:
    try:
        with open(f) as fh:
            ast.parse(fh.read())
        print(f"  [OK] {f}")
    except SyntaxError as e:
        print(f"  [FAIL] {f} : line {e.lineno} - {e.msg}")
        errors += 1

print()
if errors == 0:
    print("All files passed syntax check.")
else:
    print(f"{errors} file(s) failed.")
    sys.exit(1)

# --- 추가: 리팩토링 검증 ---
print("\n--- Refactoring Sanity Check ---")
with open("pipeline_manager.py") as fh:
    src = fh.read()

checks = [
    ("_step3_full_scoring",   "NEW Step 3 method exists"),
    ("_step4_clustering",     "Renamed Step 4 method exists"),
    ("_step5_selection_and_save", "Renamed Step 5 method exists"),
    ("_step6_visualization",  "Renamed Step 6 method exists"),
]
removed = [
    ("_step4_refinement",     "Old refinement method removed"),
    ("_step5_final_scoring",  "Old final scoring method removed"),
    ("_step7_visualization",  "Old step7 method removed"),
    ("_step3_clustering",     "Old step3_clustering method removed"),
    ("_step6_selection_and_save", "Old step6 method removed"),
]

ok = True
for name, desc in checks:
    if f"def {name}" in src:
        print(f"  [OK] {desc}")
    else:
        print(f"  [FAIL] {desc} - 'def {name}' not found!")
        ok = False

for name, desc in removed:
    if f"def {name}" in src:
        print(f"  [FAIL] {desc} - 'def {name}' still exists!")
        ok = False
    else:
        print(f"  [OK] {desc}")

# fallback 확인
if "fallback=0)" in src and "'Refinement'" in src:
    print("  [OK] Refinement config fallback present")
else:
    print("  [WARN] Refinement fallback not detected")

# 삭제된 stats 키 참조 확인
for key in ["refine_success", "total_refine_jobs", "total_final_scored"]:
    if key in src:
        print(f"  [WARN] Removed stats key '{key}' still referenced")
        ok = False

if "total_full_scored" in src:
    print("  [OK] New stats key 'total_full_scored' present")

print()
if ok:
    print("All refactoring checks passed!")
else:
    print("Some refactoring checks failed - review above.")
