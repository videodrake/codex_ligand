#!/usr/bin/env python3
"""
================================================================
  MD 자동 분석 + 그래프 생성
================================================================
  사용법:
    python analyze.py              # 대화형 선택
    python analyze.py --prefix md  # 직접 지정
    python analyze.py --plot-only  # 그래프만 다시 생성
================================================================
  자동 감지 시스템:
    - Protein only (수용액)
    - Protein + Ligand
    - Protein + Membrane
    - Protein + Ligand + Membrane
  시스템 구성에 따라 분석 항목이 자동으로 결정됩니다.
================================================================
"""

import subprocess
import os
import sys
import re
import argparse
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime


# ================================================================
#  분석 레지스트리 - 모든 분석 항목을 여기서 정의
# ================================================================
# 각 항목은 다음 필드를 가짐:
#   id        : 고유 식별자 (출력 파일명에 사용)
#   label     : 표시 이름
#   category  : "common" | "protein" | "ligand" | "membrane"
#   requires  : 실행에 필요한 조건 키 목록 (모두 충족해야 실행)
#   source    : "energy" (EDR) | "trajectory" (XTC)
#   gmx_tool  : GROMACS 도구명
#   gmx_args  : 추가 인자 (lambda로 동적 생성)
#   input_fn  : echo 입력 생성 함수
#   output_flag : 출력 파일 플래그 (-o, -od, -num 등)
#   plot      : 플롯 설정 dict
#   stat      : 통계 표시용 라벨 (None이면 통계 생략)
#   summary   : True이면 종합 요약 그래프에 포함

def build_analysis_registry():
    """시스템 구성과 무관하게, 가능한 모든 분석 항목을 정의"""
    return [
        # ==================== 공통 (에너지) ====================
        {
            'id': 'pot_energy',
            'label': 'Potential Energy',
            'category': 'common',
            'requires': ['pot_num'],
            'source': 'energy',
            'gmx_tool': 'energy',
            'gmx_args': lambda ctx: ["-f", ctx['EDR']],
            'input_fn': lambda ctx: str(ctx['pot_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Energy (kJ/mol)'},
            'stat': 'Potential (kJ/mol)',
            'summary': True,
            'description': (
                "시스템의 퍼텐셜 에너지(kJ/mol)를 시간에 따라 추적합니다.\n"
                "시뮬레이션이 안정적으로 평형에 도달했는지 확인하는 기본 지표입니다.\n"
                "값이 일정하게 수렴하면 시스템이 평형 상태임을 의미합니다.\n"
                "급격한 변동이나 지속적인 드리프트는 시스템 불안정을 나타낼 수 있습니다."
            ),
        },
        {
            'id': 'temperature',
            'label': 'Temperature',
            'category': 'common',
            'requires': ['temp_num'],
            'source': 'energy',
            'gmx_tool': 'energy',
            'gmx_args': lambda ctx: ["-f", ctx['EDR']],
            'input_fn': lambda ctx: str(ctx['temp_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Temperature (K)'},
            'stat': 'Temperature (K)',
            'summary': True,
            'description': (
                "시스템 온도(K)의 시간 변화입니다.\n"
                "thermostat(온도 조절기)이 정상 작동하는지 검증합니다.\n"
                "목표 온도(보통 300K 또는 310K) 주위로 안정적으로 요동해야 합니다.\n"
                "큰 편차는 thermostat 설정 오류나 시스템 문제를 시사합니다."
            ),
        },
        {
            'id': 'pressure',
            'label': 'Pressure',
            'category': 'common',
            'requires': ['press_num'],
            'source': 'energy',
            'gmx_tool': 'energy',
            'gmx_args': lambda ctx: ["-f", ctx['EDR']],
            'input_fn': lambda ctx: str(ctx['press_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Pressure (bar)'},
            'stat': 'Pressure (bar)',
            'summary': False,
            'description': (
                "시스템 압력(bar)의 시간 변화입니다.\n"
                "NPT 앙상블에서 barostat(압력 조절기)의 정상 작동을 확인합니다.\n"
                "압력은 자연적으로 큰 요동을 보이므로 (수백 bar) 이는 정상입니다.\n"
                "평균값이 목표 압력(보통 1 bar)에 가까운지가 중요합니다."
            ),
        },
        {
            'id': 'density',
            'label': 'Density',
            'category': 'common',
            'requires': ['dens_num'],
            'source': 'energy',
            'gmx_tool': 'energy',
            'gmx_args': lambda ctx: ["-f", ctx['EDR']],
            'input_fn': lambda ctx: str(ctx['dens_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Density (kg/m\u00b3)'},
            'stat': 'Density (kg/m\u00b3)',
            'summary': False,
            'description': (
                "시스템 밀도(kg/m3)의 시간 변화입니다.\n"
                "수용액 시스템에서 물의 실험적 밀도(~1000 kg/m3, 300K)와\n"
                "비교하여 force field와 시뮬레이션 조건의 타당성을 검증합니다.\n"
                "안정적인 밀도는 시스템이 올바르게 평형화되었음을 나타냅니다."
            ),
        },
        # ==================== 단백질 구조 ====================
        {
            'id': 'rmsd',
            'label': 'RMSD (Backbone)',
            'category': 'protein',
            'requires': ['backbone_num'],
            'source': 'trajectory',
            'gmx_tool': 'rms',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: f"{ctx['backbone_num']}\n{ctx['backbone_num']}",
            'output_flag': '-o',
            'plot': {'ylabel': 'RMSD (nm)'},
            'stat': 'RMSD Backbone (nm)',
            'summary': True,
            'description': (
                "단백질 백본(N, CA, C) 원자의 RMSD (Root Mean Square Deviation)입니다.\n"
                "초기 구조 대비 단백질의 전체적인 구조 변화를 측정합니다.\n"
                "값이 수렴하면 단백질이 안정적인 conformation에 도달했음을 의미합니다.\n"
                "일반적으로 0.1-0.3 nm이면 안정적, 0.3 nm 이상이면 큰 구조 변화가 있습니다.\n"
                "fit 그룹과 계산 그룹 모두 Backbone을 사용합니다."
            ),
        },
        {
            'id': 'rmsf',
            'label': 'RMSF per Residue',
            'category': 'protein',
            'requires': ['calpha_num'],
            'source': 'trajectory',
            'gmx_tool': 'rmsf',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'],
                                     "-res", "-n", ctx['NDX']],
            'input_fn': lambda ctx: str(ctx['calpha_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'RMSF (nm)', 'xlabel': 'Residue',
                     'time_convert': False, 'show_stats': False},
            'stat': None,
            'summary': True,
            'description': (
                "잔기별 RMSF (Root Mean Square Fluctuation)입니다. C-alpha 원자 기준.\n"
                "각 잔기의 유연성(flexibility)을 나타냅니다.\n"
                "높은 RMSF: 루프, 말단 등 유연한 영역 (주로 용매에 노출된 부분)\n"
                "낮은 RMSF: 안정적인 2차 구조 영역 (alpha-helix, beta-sheet 등)\n"
                "리간드 결합 부위 근처의 RMSF 변화는 결합에 의한 구조적 영향을 보여줍니다.\n"
                "X축은 잔기 번호, Y축은 요동 크기(nm)입니다."
            ),
        },
        {
            'id': 'gyrate',
            'label': 'Radius of Gyration',
            'category': 'protein',
            'requires': ['protein_num'],
            'source': 'trajectory',
            'gmx_tool': 'gyrate',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: str(ctx['protein_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Rg (nm)'},
            'stat': 'Rg (nm)',
            'summary': True,
            'description': (
                "단백질의 회전반경 (Radius of Gyration, Rg)입니다.\n"
                "단백질의 전체적인 크기/컴팩트함을 측정합니다.\n"
                "Rg가 안정적이면 단백질이 접힌 상태를 유지하고 있음을 의미합니다.\n"
                "Rg 증가는 단백질 펼침(unfolding), 감소는 더 컴팩트한 접힘을 시사합니다.\n"
                "Protein 그룹 전체 원자 기준으로 계산됩니다."
            ),
        },
        # ==================== 리간드 ====================
        {
            'id': 'rmsd_ligand',
            'label': 'RMSD (Ligand)',
            'category': 'ligand',
            'requires': ['ligand_num', 'backbone_num'],
            'source': 'trajectory',
            'gmx_tool': 'rms',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: f"{ctx['backbone_num']}\n{ctx['ligand_num']}",
            'output_flag': '-o',
            'plot': {'ylabel': 'RMSD (nm)'},
            'stat': 'RMSD Ligand (nm)',
            'summary': True,
            'description': (
                "리간드의 RMSD입니다. 단백질 Backbone에 fit한 후 리간드 원자의 편차를 계산합니다.\n"
                "리간드가 결합 부위에서 얼마나 움직이는지를 보여줍니다.\n"
                "낮고 안정적인 값: 리간드가 결합 포켓에 안정적으로 위치\n"
                "높거나 증가하는 값: 리간드가 결합 부위에서 이탈하거나 크게 움직임\n"
                "일반적으로 0.1-0.2 nm이면 안정적 결합, 0.3 nm 이상이면 불안정합니다."
            ),
        },
        # Ligand RMSF 제거: 소수 원자(~44개)로 least-squares fitting 시 GROMACS 크래시 발생
        {
            'id': 'mindist_prot_lig',
            'label': 'Protein-Ligand Min Distance',
            'category': 'ligand',
            'requires': ['protein_num', 'ligand_num'],
            'source': 'trajectory',
            'gmx_tool': 'mindist',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: f"{ctx['protein_num']}\n{ctx['ligand_num']}",
            'output_flag': '-od',
            'plot': {'ylabel': 'Distance (nm)'},
            'stat': 'Prot-Lig MinDist (nm)',
            'summary': True,
            'description': (
                "단백질과 리간드 사이의 최소 거리(nm)를 시간에 따라 추적합니다.\n"
                "리간드가 단백질에 결합된 상태를 유지하는지 직접적으로 확인할 수 있습니다.\n"
                "~0.2-0.4 nm: 안정적 접촉/결합 상태\n"
                "갑자기 증가: 리간드 이탈 (dissociation) 이벤트\n"
                "이 값이 꾸준히 낮게 유지되면 리간드가 결합 상태를 유지하고 있음을 의미합니다."
            ),
        },
        {
            'id': 'hbond_prot_lig',
            'label': 'Protein-Ligand H-bonds',
            'category': 'ligand',
            'requires': ['protein_num', 'ligand_num'],
            'source': 'trajectory',
            'gmx_tool': 'hbond',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: f"{ctx['protein_num']}\n{ctx['ligand_num']}",
            'output_flag': '-num',
            'plot': {'ylabel': 'Number of H-bonds'},
            'stat': 'Prot-Lig H-bonds',
            'summary': True,
            'description': (
                "단백질-리간드 간 수소결합 개수의 시간 변화입니다.\n"
                "수소결합은 리간드 결합의 핵심 상호작용 중 하나입니다.\n"
                "안정적인 수소결합 수는 리간드가 특정 결합 모드를 유지함을 나타냅니다.\n"
                "수소결합이 사라지면 리간드 결합 모드 변화나 이탈을 시사합니다.\n"
                "기준: donor-acceptor 거리 0.35nm 이하, 각도 30도 이하."
            ),
        },
        # ==================== 막 (Membrane) ====================
        {
            'id': 'density_profile',
            'label': 'Density Profile (Z-axis)',
            'category': 'membrane',
            'requires': ['system_num'],
            'source': 'trajectory',
            'gmx_tool': 'density',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'],
                                     "-n", ctx['NDX'], "-d", "Z"],
            'input_fn': lambda ctx: str(ctx['system_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Density (kg/m\u00b3)', 'xlabel': 'Z (nm)',
                     'time_convert': False, 'show_stats': False},
            'stat': None,
            'summary': True,
            'description': (
                "Z축 방향(막 법선)의 밀도 프로파일입니다.\n"
                "막의 구조적 완전성을 시각적으로 확인할 수 있습니다.\n"
                "전형적 패턴: 막 headgroup 영역에서 높은 밀도 피크 2개,\n"
                "그 사이 tail 영역에서 낮은 밀도, 양쪽 물 영역에서 ~1000 kg/m3.\n"
                "두 headgroup 피크 사이 거리가 대략적인 막 두께를 나타냅니다."
            ),
        },
        {
            'id': 'rmsd_membrane',
            'label': 'RMSD (Membrane)',
            'category': 'membrane',
            'requires': ['membrane_num'],
            'source': 'trajectory',
            'gmx_tool': 'rms',
            'gmx_args': lambda ctx: ["-s", ctx['TPR'], "-f", ctx['ANALYSIS_XTC'], "-n", ctx['NDX']],
            'input_fn': lambda ctx: f"{ctx['membrane_num']}\n{ctx['membrane_num']}",
            'output_flag': '-o',
            'plot': {'ylabel': 'RMSD (nm)'},
            'stat': 'RMSD Membrane (nm)',
            'summary': False,
            'description': (
                "막(lipid) 원자의 RMSD입니다.\n"
                "막 전체의 구조적 변화를 추적합니다.\n"
                "지질은 자연적으로 유동적이므로 단백질보다 높은 RMSD를 보이는 것이 정상입니다.\n"
                "급격한 증가는 막 구조의 붕괴나 비정상적인 변형을 나타낼 수 있습니다."
            ),
        },
        {
            'id': 'membrane_thickness',
            'label': 'Membrane Thickness (Box Z)',
            'category': 'membrane',
            'requires': ['box_z_num'],
            'source': 'energy',
            'gmx_tool': 'energy',
            'gmx_args': lambda ctx: ["-f", ctx['EDR']],
            'input_fn': lambda ctx: str(ctx['box_z_num']),
            'output_flag': '-o',
            'plot': {'ylabel': 'Box-Z (nm)'},
            'stat': 'Box-Z (nm)',
            'summary': True,
            'description': (
                "시뮬레이션 박스의 Z축 크기(nm) 변화입니다.\n"
                "막 시스템에서 Z축 방향은 막 법선과 일치하므로,\n"
                "Box-Z 변화는 막 두께 및 물 층 두께 변화를 반영합니다.\n"
                "NPT 앙상블에서 semi-isotropic 압력 커플링 시\n"
                "Z축이 독립적으로 변할 수 있으므로 모니터링이 중요합니다."
            ),
        },
    ]


# ================================================================
#  0. 분석 대상 탐색 & 선택
# ================================================================

def find_all_targets():
    """현재 디렉토리에서 분석 가능한 모든 tpr 파일을 찾아 정보와 함께 반환"""
    tpr_files = sorted(glob.glob("*.tpr"))
    if not tpr_files:
        return []

    targets = []
    for tpr in tpr_files:
        prefix = tpr.replace(".tpr", "")

        has_xtc = os.path.exists(f"{prefix}.xtc")
        has_edr = os.path.exists(f"{prefix}.edr")
        has_log = os.path.exists(f"{prefix}.log")
        has_cpt = os.path.exists(f"{prefix}.cpt")

        if not (has_xtc or has_edr):
            continue

        xtc_size = ""
        if has_xtc:
            size = os.path.getsize(f"{prefix}.xtc")
            if size > 1024**3:
                xtc_size = f"{size/(1024**3):.1f} GB"
            elif size > 1024**2:
                xtc_size = f"{size/(1024**2):.0f} MB"
            else:
                xtc_size = f"{size/1024:.0f} KB"

        sim_time = ""
        ns_day = ""
        actual_ns = 0
        if has_xtc:
            try:
                result = subprocess.run(
                    ["gmx_mpi", "check", "-f", f"{prefix}.xtc"],
                    capture_output=True, text=True, timeout=60
                )
                output = result.stdout + result.stderr
                time_match = re.search(r'Last frame\s+\d+\s+time\s+([\d.]+)', output)
                if time_match:
                    actual_ns = float(time_match.group(1)) / 1000
            except Exception:
                pass

        total_ns = 0
        finished = False
        if has_log:
            try:
                with open(f"{prefix}.log", 'r') as f:
                    content = f.read()
                dt_match = re.search(r'dt\s*=\s*([\d.e-]+)', content)
                nsteps_match = re.search(r'nsteps\s*=\s*(\d+)', content)
                if dt_match and nsteps_match:
                    total_ns = float(dt_match.group(1)) * int(nsteps_match.group(1)) / 1000
                finished = bool(re.search(r'Performance:', content))
                perf_match = re.search(r'Performance:\s+([\d.]+)', content)
                if perf_match:
                    ns_day = f"{float(perf_match.group(1)):.1f} ns/day"
            except Exception:
                pass

        if actual_ns > 0:
            if finished or abs(actual_ns - total_ns) < 1:
                sim_time = f"{actual_ns:.0f} ns \u2713"
            else:
                sim_time = f"{actual_ns:.0f}/{total_ns:.0f} ns"
        elif total_ns > 0:
            sim_time = f"{total_ns:.0f} ns"

        mtimes = []
        for ext in ['.xtc', '.edr', '.log', '.cpt']:
            path = f"{prefix}{ext}"
            if os.path.exists(path):
                mtimes.append(os.path.getmtime(path))
        mtime = max(mtimes) if mtimes else os.path.getmtime(tpr)
        date_str = datetime.fromtimestamp(mtime).strftime("%m/%d")

        targets.append({
            'prefix': prefix, 'has_xtc': has_xtc, 'has_edr': has_edr,
            'has_log': has_log, 'has_cpt': has_cpt,
            'xtc_size': xtc_size, 'sim_time': sim_time,
            'ns_day': ns_day, 'date': date_str,
        })

    def sort_key(t):
        match = re.search(r'step(\d+\.?\d*)', t['prefix'])
        return float(match.group(1)) if match else 0
    targets.sort(key=sort_key)
    return targets


def select_target(targets):
    """대화형으로 분석 대상 선택"""
    print("\n" + "=" * 70)
    print("  분석 가능한 시뮬레이션 목록")
    print("=" * 70)
    print(f"  {'#':>3}  {'Prefix':<30} {'시간':>15} {'XTC':>8} {'날짜':>6}")
    print(f"  {'-'*67}")

    for i, t in enumerate(targets, 1):
        xtc_info = t['xtc_size'] if t['has_xtc'] else '-'
        time_info = t['sim_time'] if t['sim_time'] else '-'
        print(f"  {i:>3}) {t['prefix']:<30} {time_info:>15} {xtc_info:>8} {t['date']:>6}")

    print(f"  {'-'*67}")
    print(f"  총 {len(targets)}개 발견\n")

    while True:
        try:
            choice = input("  분석할 번호를 선택하세요 (q=종료): ").strip()
            if choice.lower() == 'q':
                sys.exit(0)
            num = int(choice)
            if 1 <= num <= len(targets):
                print(f"\n  \u2192 {targets[num-1]['prefix']} 선택됨")
                return targets[num - 1]['prefix']
            print(f"  1~{len(targets)} 사이 숫자를 입력하세요.")
        except ValueError:
            print("  숫자를 입력하세요.")
        except (EOFError, KeyboardInterrupt):
            print()
            sys.exit(0)


# ================================================================
#  1. 시스템 타입 자동 감지
# ================================================================

MEMBRANE_KEYWORDS = [
    "MEMB", "DPPC", "POPC", "POPE", "DOPC", "SOLU_MEMB",
    "Lipid", "lipid", "DMPC", "DPPE", "CHOL", "DSPC", "DLPC",
    "POPS", "DPPS", "SM", "CER",
]

LIGAND_KEYWORDS = [
    "UNL", "LIG", "MOL", "Drug", "Ligand",
    "Protein_Ligand", "Protein_UNL", "Prot_Lig",
]


def detect_system_type(groups):
    """인덱스 그룹 dict에서 시스템 구성요소를 파악"""
    group_names = set(groups.keys())

    has_protein = any(g in group_names for g in ["Protein", "SOLU"])
    has_membrane = any(g in group_names for g in MEMBRANE_KEYWORDS)
    has_ligand = any(g in group_names for g in LIGAND_KEYWORDS)

    parts = []
    if has_protein:
        parts.append("Protein")
    if has_ligand:
        parts.append("Ligand")
    if has_membrane:
        parts.append("Membrane")
    sys_type = " + ".join(parts) if parts else "Unknown"

    return sys_type, has_protein, has_membrane, has_ligand


# ================================================================
#  2. 인덱스 자동 생성
# ================================================================

def ensure_full_index(gmx, tpr_file, original_ndx):
    full_ndx = "index_analysis.ndx"
    required_groups = ["Backbone", "C-alpha", "Protein"]

    if os.path.exists(full_ndx):
        with open(full_ndx, 'r') as f:
            existing = set(re.findall(r'\[\s*(.+?)\s*\]', f.read()))
        missing = [g for g in required_groups if g not in existing]
        if not missing:
            print(f"  {full_ndx} 이미 존재 (필수 그룹 포함)")
            return full_ndx
        print(f"  {full_ndx}에 {', '.join(missing)} 없음 \u2192 재생성")

    print("  표준 인덱스 그룹 생성 중...")
    try:
        subprocess.run(
            [gmx, "make_ndx", "-f", tpr_file, "-o", full_ndx],
            input="q\n", capture_output=True, text=True, timeout=60
        )
    except Exception as e:
        print(f"  [WARNING] make_ndx 실패: {e}")
        return original_ndx

    if not os.path.exists(full_ndx):
        return original_ndx

    # 기존 ndx에서 커스텀 그룹 병합
    if os.path.exists(original_ndx):
        with open(full_ndx, 'r') as f:
            standard_groups = set(re.findall(r'\[\s*(.+?)\s*\]', f.read()))

        custom_block = ""
        copying = False
        with open(original_ndx, 'r') as f:
            for line in f:
                match = re.match(r'\[\s*(.+?)\s*\]', line)
                if match:
                    copying = match.group(1) not in standard_groups
                    if copying:
                        custom_block += "\n" + line
                elif copying:
                    custom_block += line

        if custom_block:
            with open(full_ndx, 'a') as f:
                f.write(custom_block)
            print(f"  커스텀 그룹 {custom_block.count('[')}개 병합 완료")

    with open(full_ndx, 'r') as f:
        all_groups = re.findall(r'\[\s*(.+?)\s*\]', f.read())
    print(f"  {full_ndx}: 총 {len(all_groups)}개 그룹")
    for g in required_groups:
        print(f"    {'✓' if g in all_groups else '✗'} {g}")

    return full_ndx


# ================================================================
#  3. 그룹/에너지 탐색
# ================================================================

def parse_ndx_groups(ndx_file):
    groups = {}
    idx = 0
    with open(ndx_file, 'r') as f:
        for line in f:
            match = re.match(r'\[\s*(.+?)\s*\]', line)
            if match:
                groups[match.group(1)] = idx
                idx += 1
    return groups


def find_group(groups, *keywords):
    """정확 매치 → 대소문자 무시 → 부분 매치 순으로 검색"""
    for kw in keywords:
        if kw in groups:
            return groups[kw], kw
        for name, num in groups.items():
            if name.lower() == kw.lower():
                return num, name
        for name, num in groups.items():
            if kw.lower() in name.lower():
                return num, name
    return None, None


def find_energy_terms(gmx, edr_file):
    """EDR 파일에서 사용 가능한 에너지 항목 목록을 파싱"""
    try:
        result = subprocess.run(
            [gmx, "energy", "-f", edr_file],
            input="0\n", capture_output=True, text=True, timeout=30
        )
        output = result.stdout + result.stderr
    except Exception as e:
        print(f"  [WARNING] gmx energy 파싱 실패: {e}")
        return {}

    terms = {}
    for m in re.finditer(r'(\d+)\s+([\w#()\-]+(?:\s[\w#()\-]+)*)', output):
        terms[m.group(2).strip()] = int(m.group(1))
    return terms


def find_energy_term(energy_terms, *keywords):
    for kw in keywords:
        for name, num in energy_terms.items():
            if kw.lower() == name.lower():
                return num, name
        for name, num in energy_terms.items():
            if kw.lower() in name.lower():
                return num, name
    return None, None


# ================================================================
#  4. GROMACS 명령 실행
# ================================================================

def run_gmx(gmx, args, echo_input=None, label=""):
    cmd = [gmx] + args
    cmd_str = ' '.join(cmd)
    if echo_input:
        cmd_str = f'echo "{echo_input.replace(chr(10), " ")}" | {cmd_str}'

    print(f"\n  [{label}]")
    print(f"  $ {cmd_str}")

    try:
        result = subprocess.run(
            cmd,
            input=echo_input + "\n" if echo_input else None,
            capture_output=True, text=True, timeout=3600
        )
        if result.returncode != 0:
            stderr = result.stderr
            if "Fatal error" in stderr:
                fatal = re.search(r'Fatal error:(.+?)(?:\n\n|$)', stderr, re.DOTALL)
                print(f"  [ERROR] {fatal.group(1).strip() if fatal else stderr[-300:]}")
                return False
            elif result.returncode != 0:
                print(f"  [WARNING] returncode={result.returncode}")
                return False
        print("  \u2713 완료")
        return True
    except subprocess.TimeoutExpired:
        print("  [ERROR] 타임아웃 (1시간)")
        return False
    except FileNotFoundError:
        print(f"  [ERROR] '{gmx}' 명령을 찾을 수 없습니다. --gmx 옵션을 확인하세요.")
        return False
    except Exception as e:
        print(f"  [ERROR] {e}")
        return False


# ================================================================
#  5. XVG 파싱 & 플롯
# ================================================================

def parse_xvg(filename):
    data = []
    x_label, y_label, title = "X", "Y", os.path.basename(filename).replace('.xvg', '')
    legends = []

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            elif line.startswith('@'):
                if 'xaxis' in line and 'label' in line and '"' in line:
                    x_label = line.split('"')[1]
                elif 'yaxis' in line and 'label' in line and '"' in line:
                    y_label = line.split('"')[1]
                elif 'title' in line and 'legend' not in line and '"' in line:
                    title = line.split('"')[1]
                elif 'legend' in line and '"' in line:
                    legends.append(line.split('"')[1])
            else:
                vals = line.split()
                if vals:
                    try:
                        data.append([float(v) for v in vals])
                    except ValueError:
                        continue

    return np.array(data) if data else np.array([]), x_label, y_label, title, legends


def calc_stats(y):
    return {'mean': np.mean(y), 'std': np.std(y),
            'min': np.min(y), 'max': np.max(y)}


def plot_xvg(filename, title=None, xlabel=None, ylabel=None,
             time_convert=True, multiplot=False, show_stats=True,
             running_avg_ns=1.0):
    data, x_lab, y_lab, title_raw, legends = parse_xvg(filename)
    if data.size == 0:
        print(f"  [SKIP] {filename} - 데이터 없음")
        return

    title = title or title_raw
    xlabel = xlabel or x_lab
    ylabel = ylabel or y_lab

    x = data[:, 0]
    if time_convert and x.max() > 1000:
        x = x / 1000.0
        xlabel = "Time (ns)"

    n_cols = data.shape[1] - 1
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    if multiplot and n_cols > 1:
        fig, axes = plt.subplots(n_cols, 1, figsize=(10, 3 * n_cols), sharex=True)
        if n_cols == 1:
            axes = [axes]
        for i in range(n_cols):
            ax = axes[i]
            y = data[:, i + 1]
            label = legends[i] if i < len(legends) else f"Coord {i+1}"
            ax.plot(x, y, color=colors[i % len(colors)], linewidth=0.5, alpha=0.7)
            dt = x[1] - x[0] if len(x) > 1 else 1
            window = max(1, int(running_avg_ns / dt))
            if len(y) > window:
                avg = np.convolve(y, np.ones(window)/window, mode='valid')
                ax.plot(x[:len(avg)], avg, color='black', linewidth=1.5,
                        label=f'{running_avg_ns:.0f} ns avg')
            if show_stats:
                ax.axhline(calc_stats(y)['mean'], color='red',
                           linestyle='--', alpha=0.4, linewidth=0.8)
            ax.set_ylabel(ylabel)
            ax.set_title(label)
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(True, alpha=0.3)
        axes[-1].set_xlabel(xlabel)
    else:
        fig, ax = plt.subplots(figsize=(10, 5))
        for i in range(n_cols):
            y = data[:, i + 1]
            lbl = legends[i] if i < len(legends) else (None if n_cols == 1 else f"Col {i+1}")
            ax.plot(x, y, linewidth=0.8, label=lbl, alpha=0.7,
                    color=colors[i % len(colors)])

            if time_convert and len(x) > 10:
                dt = x[1] - x[0] if len(x) > 1 else 1
                window = max(1, int(running_avg_ns / dt))
                if len(y) > window:
                    avg = np.convolve(y, np.ones(window)/window, mode='valid')
                    ax.plot(x[:len(avg)], avg, color='black', linewidth=1.5,
                            alpha=0.8, label=f'{running_avg_ns:.0f} ns avg' if i == 0 else None)

        if show_stats and n_cols == 1:
            s = calc_stats(data[:, 1])
            ax.set_title(f"{title}\n(Mean: {s['mean']:.4f} | Std: {s['std']:.4f})", fontsize=11)
        else:
            ax.set_title(title)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        if legends or (time_convert and len(x) > 10):
            ax.legend(fontsize=8)

    plt.tight_layout()
    outname = filename.replace('.xvg', '.png')
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  \u2713 {outname}")


def plot_summary(prefix, items, out_dir="."):
    """실제 생성된 분석 결과들만 종합 그래프로"""
    available = [(a, p) for a, p in items if os.path.exists(p)]
    if len(available) < 2:
        return

    n = len(available)
    cols = min(n, 3)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes = np.array(axes).flatten() if n > 1 else [axes]

    for idx, (analysis, xvg_file) in enumerate(available):
        if idx >= len(axes):
            break
        ax = axes[idx]
        data, x_lab, y_lab, _, _ = parse_xvg(xvg_file)
        if data.size == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(analysis['label'])
            continue

        x = data[:, 0]
        is_timeseries = analysis['plot'].get('time_convert', True)
        if is_timeseries and x.max() > 1000:
            x = x / 1000.0
            x_lab = "Time (ns)"

        y = data[:, 1]
        ax.plot(x, y, linewidth=0.6, alpha=0.7, color='#1f77b4')

        if is_timeseries and len(x) > 10:
            dt = x[1] - x[0] if len(x) > 1 else 1
            window = max(1, int(1.0 / dt))
            if len(y) > window:
                avg = np.convolve(y, np.ones(window)/window, mode='valid')
                ax.plot(x[:len(avg)], avg, color='red', linewidth=1.2, alpha=0.8)

        s = calc_stats(y)
        ax.set_title(f"{analysis['label']}\n(Mean={s['mean']:.3f}, Std={s['std']:.3f})", fontsize=9)
        ax.set_xlabel(analysis['plot'].get('xlabel', x_lab), fontsize=8)
        ax.set_ylabel(analysis['plot'].get('ylabel', y_lab), fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

    for idx in range(len(available), len(axes)):
        fig.delaxes(axes[idx])

    fig.suptitle(f"{prefix} - MD Analysis Summary", fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    outname = os.path.join(out_dir, f"{prefix}_summary.png")
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  \u2713 {outname} (종합 요약)")


# ================================================================
#  7. 분석 리포트 생성
# ================================================================

def write_report(prefix, sys_type, generated, out_dir="."):
    """분석 결과에 대한 설명 리포트를 텍스트 파일로 생성"""
    report_path = os.path.join(out_dir, f"{prefix}_analysis_report.txt")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cat_labels = {
        'common': '시스템 검증 (System Validation)',
        'protein': '단백질 구조 분석 (Protein Structure)',
        'ligand': '리간드 상호작용 분석 (Ligand Interaction)',
        'membrane': '막 분석 (Membrane)',
    }

    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write(f"  MD Analysis Report\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Date       : {now}\n")
        f.write(f"  Prefix     : {prefix}\n")
        f.write(f"  System     : {sys_type}\n")
        f.write(f"  Analyses   : {len(generated)}개\n")
        f.write("=" * 70 + "\n\n")

        # 입력 파일 설명
        f.write("-" * 70 + "\n")
        f.write("  입력 파일 (Input Files)\n")
        f.write("-" * 70 + "\n\n")
        input_files = [
            (f"{prefix}.tpr", "GROMACS run input file. 시뮬레이션의 모든 매개변수, 토폴로지, 좌표를 포함하는 바이너리 파일."),
            (f"{prefix}.xtc", "압축된 trajectory 파일. 시뮬레이션 동안의 모든 원자 좌표를 시간순으로 저장."),
            (f"{prefix}.edr", "에너지 파일. 시뮬레이션 동안의 에너지, 온도, 압력 등 열역학적 데이터를 저장."),
            (f"{prefix}.gro", "GROMACS 좌표 파일. 시스템의 원자 좌표와 박스 크기를 텍스트 형식으로 저장."),
            ("index.ndx", "인덱스 파일. 원자 그룹 정의 (Protein, Ligand, Membrane, Water 등)."),
        ]
        for fname, desc in input_files:
            if os.path.exists(os.path.join(out_dir, fname)) or os.path.exists(fname):
                f.write(f"  {fname}\n")
                f.write(f"    {desc}\n\n")

        # ChimeraX 파일 설명
        nw_gro = os.path.join(out_dir, f"{prefix}_nw.gro")
        nw_xtc = os.path.join(out_dir, f"{prefix}_nw.xtc")
        if os.path.exists(nw_gro) or os.path.exists(nw_xtc):
            f.write("-" * 70 + "\n")
            f.write("  시각화 파일 (Visualization Files)\n")
            f.write("-" * 70 + "\n\n")
            f.write(f"  {prefix}_nw.gro\n")
            f.write(f"    물과 이온을 제거한 좌표 파일. ChimeraX/VMD에서 구조를 로드할 때 사용.\n\n")
            f.write(f"  {prefix}_nw.xtc\n")
            f.write(f"    물과 이온을 제거한 trajectory. ChimeraX/VMD에서 동적 시각화에 사용.\n")
            f.write(f"    단백질, 리간드, 막 등 관심 구성요소만 포함됩니다.\n\n")
            f.write(f"  사용법 (ChimeraX):\n")
            f.write(f"    open {prefix}_nw.gro\n")
            f.write(f"    open {prefix}_nw.xtc structureModel #1\n\n")

        # 카테고리별 분석 결과 설명
        current_cat = None
        for a, xvg_path in generated:
            cat = a.get('category', 'common')
            if cat != current_cat:
                current_cat = cat
                f.write("-" * 70 + "\n")
                f.write(f"  {cat_labels.get(cat, cat)}\n")
                f.write("-" * 70 + "\n\n")

            xvg_name = os.path.basename(xvg_path)
            png_name = xvg_name.replace('.xvg', '.png')

            f.write(f"  [{a['label']}]\n")
            f.write(f"  데이터 : {xvg_name}\n")
            f.write(f"  그래프 : {png_name}\n")

            # 통계
            if a.get('stat') and os.path.exists(xvg_path):
                data, _, _, _, _ = parse_xvg(xvg_path)
                if data.size > 0:
                    s = calc_stats(data[:, 1])
                    f.write(f"  통계   : Mean={s['mean']:.4f}, Std={s['std']:.4f}, "
                            f"Min={s['min']:.4f}, Max={s['max']:.4f}\n")

            # 설명
            desc = a.get('description', '')
            if desc:
                f.write("\n")
                for line in desc.strip().split('\n'):
                    f.write(f"    {line}\n")

            f.write("\n")

        # 종합 요약
        summary_path = os.path.join(out_dir, f"{prefix}_summary.png")
        if os.path.exists(summary_path):
            f.write("-" * 70 + "\n")
            f.write("  종합 요약 (Summary)\n")
            f.write("-" * 70 + "\n\n")
            f.write(f"  {prefix}_summary.png\n")
            f.write(f"    모든 주요 분석 결과를 한 눈에 볼 수 있는 종합 그래프입니다.\n")
            f.write(f"    빨간 선은 이동평균(running average)을 나타냅니다.\n\n")

        f.write("=" * 70 + "\n")
        f.write("  End of Report\n")
        f.write("=" * 70 + "\n")

    print(f"  \u2713 {report_path}")
    return report_path


# ================================================================
#  메인
# ================================================================

def main():
    parser = argparse.ArgumentParser(description="MD 자동 분석 + 그래프")
    parser.add_argument("--prefix", default=None, help="파일 접두사 (미지정 시 대화형 선택)")
    parser.add_argument("--ndx", default="index.ndx", help="원본 인덱스 파일")
    parser.add_argument("--gmx", default="gmx_mpi", help="GROMACS 실행 파일명")
    parser.add_argument("--plot-only", action="store_true", help="분석 건너뛰고 그래프만")
    parser.add_argument("--outdir", default=None, help="결과 출력 디렉토리")
    args = parser.parse_args()

    # --- 대상 선택 ---
    if args.prefix is None:
        targets = find_all_targets()
        if not targets:
            print("[ERROR] 분석 가능한 시뮬레이션이 없습니다. (.tpr + .xtc 또는 .edr 필요)")
            sys.exit(1)
        elif len(targets) == 1:
            args.prefix = targets[0]['prefix']
            print(f"  자동 선택: {args.prefix}")
        else:
            args.prefix = select_target(targets)

    PREFIX = args.prefix
    TPR = f"{PREFIX}.tpr"
    XTC = f"{PREFIX}.xtc"
    EDR = f"{PREFIX}.edr"
    GMX = args.gmx
    OUT = args.outdir or "."
    if args.outdir:
        os.makedirs(OUT, exist_ok=True)

    print("\n" + "=" * 60)
    print("  MD Analysis")
    print("=" * 60)
    print(f"  Prefix : {PREFIX}")
    for label, path in [("TPR", TPR), ("XTC", XTC), ("EDR", EDR)]:
        print(f"  {label:<7}: {path} {'✓' if os.path.exists(path) else '✗'}")
    nopbc = f"{PREFIX}_noPBC.xtc"
    if os.path.exists(nopbc):
        print(f"  noPBC  : {nopbc} ✓")
    if args.outdir:
        print(f"  Output : {OUT}")
    print("=" * 60)

    # ================================================================
    #  컨텍스트 구축 - 시스템 감지 & 그룹 탐색
    # ================================================================
    if not args.plot_only:
        for f in [TPR, XTC, EDR]:
            if not os.path.exists(f):
                print(f"[ERROR] 파일 없음: {f}")
                sys.exit(1)

        # 인덱스 준비
        print("\n[ 인덱스 준비 ]")
        NDX = ensure_full_index(GMX, TPR, args.ndx)

        # 그룹 파싱 & 시스템 감지
        groups = parse_ndx_groups(NDX)
        sys_type, has_protein, has_membrane, has_ligand = detect_system_type(groups)

        print(f"\n[ 시스템 감지: {sys_type} ]")

        # 모든 그룹 탐색 → ctx에 저장
        ctx = {
            'PREFIX': PREFIX, 'TPR': TPR, 'XTC': XTC, 'EDR': EDR,
            'GMX': GMX, 'NDX': NDX, 'OUT': OUT,
            'has_protein': has_protein,
            'has_membrane': has_membrane,
            'has_ligand': has_ligand,
        }

        # 단백질 그룹
        group_searches = [
            ('backbone_num', 'backbone_name', ["Backbone", "BB", "MainChain"]),
            ('calpha_num',   'calpha_name',   ["C-alpha", "CA", "Calpha"]),
            ('protein_num',  'protein_name',  ["Protein", "SOLU", "protein"]),
            ('system_num',   'system_name',   ["System"]),
        ]
        # 막 그룹
        if has_membrane:
            group_searches += [
                ('solumemb_num', 'solumemb_name', ["SOLU_MEMB", "Protein_MEMB"]),
                ('membrane_num', 'membrane_name', ["MEMB", "DPPC", "POPC", "POPE", "DOPC",
                                                    "DMPC", "DPPE", "Lipid"]),
            ]
        # 리간드 그룹
        if has_ligand:
            group_searches += [
                ('ligand_num',  'ligand_name',  ["UNL", "LIG", "MOL", "Drug", "Ligand"]),
                ('protlig_num', 'protlig_name', ["Protein_Ligand", "Protein_UNL", "Prot_Lig"]),
            ]

        print("\n  [ 그룹 ]")
        for num_key, name_key, keywords in group_searches:
            num, name = find_group(groups, *keywords)
            ctx[num_key] = num
            ctx[name_key] = name
            if num is not None:
                print(f"    [{num:>2}] {name}")

        # 에너지 항목 탐색
        energy_terms = find_energy_terms(GMX, EDR)
        energy_searches = [
            ('pot_num',   ["Potential"]),
            ('temp_num',  ["Temperature"]),
            ('press_num', ["Pressure"]),
            ('dens_num',  ["Density"]),
        ]
        if has_membrane:
            energy_searches.append(('box_z_num', ["Box-Z"]))

        print("\n  [ 에너지 항목 ]")
        for num_key, keywords in energy_searches:
            num, name = find_energy_term(energy_terms, *keywords)
            ctx[num_key] = num
            if num is not None:
                print(f"    [{num:>2}] {name}")

        # ============================================================
        #  ChimeraX 추출 그룹 결정 & trajectory 생성
        # ============================================================
        # 우선순위: SOLU_MEMB → Protein_Ligand → Protein
        chimera_grp, chimera_name = None, None
        if has_membrane and ctx.get('solumemb_num') is not None:
            chimera_grp, chimera_name = ctx['solumemb_num'], ctx['solumemb_name']
        elif has_ligand and ctx.get('protlig_num') is not None:
            chimera_grp, chimera_name = ctx['protlig_num'], ctx['protlig_name']
        elif ctx.get('protein_num') is not None:
            chimera_grp, chimera_name = ctx['protein_num'], ctx['protein_name']

        NW_XTC = os.path.join(OUT, f"{PREFIX}_nw.xtc")
        NW_GRO = os.path.join(OUT, f"{PREFIX}_nw.gro")

        # 분석용 trajectory: noPBC가 있으면 우선 사용, 없으면 원본 XTC
        NOPBC_XTC = f"{PREFIX}_noPBC.xtc"
        if os.path.exists(NOPBC_XTC):
            ctx['ANALYSIS_XTC'] = NOPBC_XTC
            print(f"\n  분석용 trajectory: {NOPBC_XTC} (noPBC)")
        else:
            ctx['ANALYSIS_XTC'] = XTC
            print(f"\n  분석용 trajectory: {XTC}")

        if chimera_grp is not None:
            print(f"\n  ChimeraX 추출 그룹: [{chimera_grp}] {chimera_name}")

            # fit 그룹: Backbone 우선, 없으면 Protein
            fit_grp = ctx.get('backbone_num')
            fit_name = ctx.get('backbone_name', 'Backbone')
            if fit_grp is None:
                fit_grp = ctx.get('protein_num')
                fit_name = ctx.get('protein_name', 'Protein')

            # Step 1: PBC 보정 + 단백질 중심 정렬 (전체 시스템, 임시 파일)
            TMP_CENTER = os.path.join(OUT, f"{PREFIX}_tmp_center.xtc")

            # Step 2: Backbone fit + 시각화 그룹 추출 → NW_XTC
            # Step 3: 첫 프레임 gro (fitted)

            regen_nw = not os.path.exists(NW_XTC) or \
                       os.path.getmtime(NW_XTC) < os.path.getmtime(XTC)
            regen_gro = not os.path.exists(NW_GRO) or \
                        os.path.getmtime(NW_GRO) < os.path.getmtime(XTC)

            if regen_nw or regen_gro:
                # Step 1: center + PBC
                # 입력: 중심 그룹 (Protein) + 출력 그룹 (System 전체)
                center_grp = ctx.get('protein_num', 0)
                system_grp = ctx.get('system_num', 0)
                run_gmx(GMX,
                    ["trjconv", "-s", TPR, "-f", XTC, "-o", TMP_CENTER,
                     "-pbc", "mol", "-center", "-n", NDX],
                    echo_input=f"{center_grp}\n{system_grp}",
                    label="PBC 보정 + 단백질 중심 정렬")

                if os.path.exists(TMP_CENTER):
                    # Step 2: fit to backbone + extract chimera group
                    if regen_nw:
                        run_gmx(GMX,
                            ["trjconv", "-s", TPR, "-f", TMP_CENTER, "-o", NW_XTC,
                             "-fit", "rot+trans", "-n", NDX],
                            echo_input=f"{fit_grp}\n{chimera_grp}",
                            label=f"Backbone fit + 추출 ({chimera_name})")

                    # Step 3: fitted gro (첫 프레임)
                    if regen_gro:
                        run_gmx(GMX,
                            ["trjconv", "-s", TPR, "-f", TMP_CENTER, "-o", NW_GRO,
                             "-fit", "rot+trans", "-dump", "0", "-n", NDX],
                            echo_input=f"{fit_grp}\n{chimera_grp}",
                            label=f"Backbone fit + gro 추출 ({chimera_name})")

                    # 임시 파일 정리
                    if os.path.exists(TMP_CENTER):
                        os.remove(TMP_CENTER)
                        print(f"  임시 파일 {TMP_CENTER} 삭제")
            else:
                print(f"\n  {NW_XTC}, {NW_GRO} 최신 상태, 건너뜀")

        # ============================================================
        #  분석 레지스트리 기반 실행
        # ============================================================
        registry = build_analysis_registry()

        # 카테고리 필터: 시스템 구성에 맞는 것만 선택
        active_categories = {"common"}
        if has_protein:
            active_categories.add("protein")
        if has_ligand:
            active_categories.add("ligand")
        if has_membrane:
            active_categories.add("membrane")

        active_analyses = [a for a in registry if a['category'] in active_categories]

        # requires 조건 체크 → 실행 가능한 것만
        runnable = []
        for a in active_analyses:
            missing = [r for r in a['requires'] if ctx.get(r) is None]
            if not missing:
                runnable.append(a)
            else:
                print(f"  [SKIP] {a['label']} - 필요 그룹 없음: {', '.join(missing)}")

        # 카테고리별 구분 출력 & 실행
        print(f"\n  실행할 분석: {len(runnable)}개")

        current_cat = None
        for a in runnable:
            if a['category'] != current_cat:
                current_cat = a['category']
                cat_labels = {
                    'common': '시스템 검증',
                    'protein': '단백질 구조',
                    'ligand': '리간드 상호작용',
                    'membrane': '막 분석',
                }
                print(f"\n{'=' * 40}")
                print(f"  {cat_labels.get(current_cat, current_cat)}")
                print(f"{'=' * 40}")

            outfile = os.path.join(OUT, f"{PREFIX}_{a['id']}.xvg")
            gmx_args = a['gmx_args'](ctx) + [a['output_flag'], outfile]

            run_gmx(GMX, [a['gmx_tool']] + gmx_args,
                    echo_input=a['input_fn'](ctx),
                    label=a['label'])

        print(f"\n{'=' * 60}")
        print("  분석 명령 완료!")
        print(f"{'=' * 60}")

    # ================================================================
    #  그래프 생성 (plot-only 모드에서도 동작)
    # ================================================================
    print("\n[ 그래프 생성 ]")

    # plot-only에서는 NDX가 없으므로, 존재하는 xvg 기반으로 동작
    registry = build_analysis_registry()

    generated = []
    for a in registry:
        xvg_path = os.path.join(OUT, f"{PREFIX}_{a['id']}.xvg")
        if os.path.exists(xvg_path):
            plot_kwargs = dict(a['plot'])
            plot_kwargs['title'] = f"{PREFIX} - {a['label']}"
            plot_xvg(xvg_path, **plot_kwargs)
            generated.append((a, xvg_path))

    # Pull 분석
    for pf in set(glob.glob(os.path.join(OUT, f"{PREFIX}*pullx.xvg"))):
        plot_xvg(pf, title=f"{PREFIX} - Pull Distance",
                 ylabel="Distance (nm)", multiplot=True)
        generated.append(({'label': 'Pull', 'plot': {}, 'stat': None, 'summary': False}, pf))

    print(f"\n  총 {len(generated)}개 그래프 생성")

    # --- 종합 요약 그래프 ---
    summary_items = [(a, p) for a, p in generated if a.get('summary', False)]
    if len(summary_items) >= 2:
        plot_summary(PREFIX, summary_items, OUT)

    # --- 통계 요약 ---
    stat_items = [(a, p) for a, p in generated if a.get('stat')]
    if stat_items:
        print(f"\n[ 통계 요약 ]")
        print(f"  {'항목':<30} {'Mean':>12} {'Std':>12} {'Min':>12} {'Max':>12}")
        print(f"  {'-'*78}")
        for a, path in stat_items:
            data, _, _, _, _ = parse_xvg(path)
            if data.size > 0:
                s = calc_stats(data[:, 1])
                print(f"  {a['stat']:<30} {s['mean']:>12.4f} {s['std']:>12.4f} "
                      f"{s['min']:>12.4f} {s['max']:>12.4f}")

    # --- ChimeraX 안내 ---
    NW_GRO = os.path.join(OUT, f"{PREFIX}_nw.gro")
    NW_XTC = os.path.join(OUT, f"{PREFIX}_nw.xtc")
    if os.path.exists(NW_GRO) and os.path.exists(NW_XTC):
        gro_mb = os.path.getsize(NW_GRO) / (1024**2)
        xtc_mb = os.path.getsize(NW_XTC) / (1024**2)
        print(f"\n  ChimeraX 트래젝토리 보기:")
        print(f"    {NW_GRO} ({gro_mb:.0f} MB)")
        print(f"    {NW_XTC} ({xtc_mb:.0f} MB)")
        print(f"    명령:")
        print(f"      open {NW_GRO}")
        print(f"      open {NW_XTC} structureModel #1")

    # --- 분석 리포트 생성 ---
    if generated:
        print("\n[ 분석 리포트 생성 ]")
        # plot-only 모드에서는 ndx에서 시스템 타입 재감지
        if args.plot_only:
            ndx_path = args.ndx
            if os.path.exists("index_analysis.ndx"):
                ndx_path = "index_analysis.ndx"
            if os.path.exists(ndx_path):
                _groups = parse_ndx_groups(ndx_path)
                sys_type, _, _, _ = detect_system_type(_groups)
            else:
                sys_type = "Unknown"
        write_report(PREFIX, sys_type, generated, OUT)

    print(f"\n{'=' * 60}")
    print("  모든 분석 완료!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
