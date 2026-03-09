#!/usr/bin/env python3
"""
================================================================
  Protein-Ligand 접촉 분석 (MDAnalysis)
================================================================
  - 리간드 cutoff 이내 잔기 목록 + 접촉 빈도 (occupancy)
  - 잔기별 최소 거리 시계열
  - 수소결합 분석

  사용법:
    python analyze_ligand_contacts.py                          # 자동 감지
    python analyze_ligand_contacts.py --prefix md              # prefix 지정
    python analyze_ligand_contacts.py -l "resname LIG"         # 리간드 지정
    python analyze_ligand_contacts.py --start 50 --stop 200    # 구간 지정
================================================================
"""

import argparse
import os
import glob
import numpy as np
import warnings
warnings.filterwarnings('ignore')

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
except ImportError:
    print("[ERROR] MDAnalysis 미설치. 다음 명령어로 설치:")
    print("  pip install MDAnalysis --break-system-packages")
    print("  # 또는 conda install -c conda-forge mdanalysis")
    exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError:
    HAS_PLT = False
    print("[WARNING] matplotlib 없음. 그래프 생성 생략.")


# ============================================================
# 0. 파일 자동 감지
# ============================================================
def find_prefix():
    """현재 디렉토리에서 tpr+xtc 조합을 찾아 prefix 반환"""
    tpr_files = sorted(glob.glob("*.tpr"))
    for tpr in tpr_files:
        prefix = tpr.replace(".tpr", "")
        if os.path.exists(f"{prefix}.xtc"):
            return prefix
    return None


def find_trajectory(prefix):
    """분석용 trajectory 선택: noPBC > 원본 XTC"""
    candidates = [
        f"{prefix}_noPBC.xtc",
        f"{prefix}.xtc",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def detect_ligand_resname(u):
    """Universe에서 리간드 resname 자동 감지"""
    # 단백질/물/이온이 아닌 잔기를 찾음
    known_solvent = {'SOL', 'HOH', 'WAT', 'TIP3', 'TIP4', 'SPC'}
    known_ions = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'NA+', 'CL-'}

    candidates = set()
    for res in u.residues:
        rn = res.resname
        if rn in known_solvent or rn in known_ions:
            continue
        # 단백질 잔기(표준 아미노산) 제외
        try:
            if res.atoms.select_atoms("protein"):
                continue
        except Exception:
            pass
        candidates.add(rn)

    # 일반적인 리간드 이름 우선
    priority = ['UNL', 'LIG', 'MOL', 'DRG', 'INH']
    for name in priority:
        if name in candidates:
            return name

    # 남은 후보 중 첫 번째
    if candidates:
        return sorted(candidates)[0]
    return None


# ============================================================
# 1. 접촉 빈도 (Occupancy) 분석
# ============================================================
def analyze_contacts(u, lig_sel, cutoff=5.0, start_ns=None, stop_ns=None):
    """
    매 프레임마다 리간드 cutoff A 이내 단백질 잔기를 카운트.
    Returns:
        contact_counts: {(resid, resname): frame_count}
        total_frames: int
        per_frame_data: {(resid, resname): [min_dist_per_frame]}
        times: np.array
    """
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(lig_sel)

    if len(ligand) == 0:
        print(f"[ERROR] 리간드 선택 '{lig_sel}'에 해당하는 원자 없음!")
        all_resnames = sorted(set(u.atoms.resnames))
        print(f"  사용 가능한 resname: {all_resnames}")
        exit(1)

    print(f"  Protein : {len(protein)} atoms ({len(protein.residues)} residues)")
    print(f"  Ligand  : {len(ligand)} atoms (resname: {sorted(set(ligand.resnames))})")

    # 프레임 범위 설정
    dt_ns = u.trajectory.dt / 1000.0  # ps -> ns
    start_frame = 0
    stop_frame = len(u.trajectory)

    if start_ns is not None:
        start_frame = int(start_ns / dt_ns)
    if stop_ns is not None:
        stop_frame = min(int(stop_ns / dt_ns), len(u.trajectory))

    total_frames = stop_frame - start_frame
    total_ns = total_frames * dt_ns
    print(f"  범위    : frame {start_frame}-{stop_frame} ({total_frames} frames, {total_ns:.1f} ns)")

    # 잔기별 접촉 카운트 + 최소 거리 시계열
    contact_counts = {}
    per_frame_mindist = {}
    times = []

    for i, ts in enumerate(u.trajectory[start_frame:stop_frame]):
        if i % 100 == 0:
            pct = i / total_frames * 100
            print(f"  진행: {i}/{total_frames} ({pct:.0f}%) - {i*dt_ns:.1f} ns", end='\r')

        times.append(ts.time / 1000.0)

        for res in protein.residues:
            key = (res.resid, res.resname)
            res_atoms = res.atoms

            dist_matrix = distances.distance_array(
                res_atoms.positions, ligand.positions, box=u.dimensions
            )
            min_dist = dist_matrix.min()

            if key not in per_frame_mindist:
                per_frame_mindist[key] = []
            per_frame_mindist[key].append(min_dist)

            if min_dist <= cutoff:
                contact_counts[key] = contact_counts.get(key, 0) + 1

    print(f"\n  완료: {len(contact_counts)}개 잔기 접촉 감지")

    return contact_counts, total_frames, per_frame_mindist, np.array(times)


# ============================================================
# 2. 수소결합 분석
# ============================================================
def analyze_hbonds(u, lig_sel, start_ns=None, stop_ns=None):
    """리간드-단백질 수소결합 분석"""

    dt_ns = u.trajectory.dt / 1000.0
    start_frame = int(start_ns / dt_ns) if start_ns else None
    stop_frame = int(stop_ns / dt_ns) if stop_ns else None

    # Protein -> Ligand (protein이 donor)
    hbonds_pl = HydrogenBondAnalysis(
        u,
        donors_sel="protein",
        acceptors_sel=lig_sel,
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=120,
    )

    # Ligand -> Protein (ligand가 donor)
    hbonds_lp = HydrogenBondAnalysis(
        u,
        donors_sel=lig_sel,
        acceptors_sel="protein",
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=120,
    )

    print("\n  수소결합 분석 (Protein -> Ligand)...")
    hbonds_pl.run(start=start_frame, stop=stop_frame, verbose=True)

    print("  수소결합 분석 (Ligand -> Protein)...")
    hbonds_lp.run(start=start_frame, stop=stop_frame, verbose=True)

    return hbonds_pl, hbonds_lp


# ============================================================
# 3. 결과 출력
# ============================================================
def write_results(outdir, prefix, contact_counts, total_frames,
                  per_frame_mindist, times, hbonds_pl, hbonds_lp,
                  u, lig_sel, cutoff):

    os.makedirs(outdir, exist_ok=True)

    # --- 3.1 접촉 빈도 테이블 ---
    occ_file = os.path.join(outdir, f"{prefix}_contact_occupancy.csv")
    sorted_contacts = sorted(contact_counts.items(), key=lambda x: -x[1])

    with open(occ_file, 'w') as f:
        f.write("resid,resname,contact_frames,occupancy(%),avg_min_dist(A),std_min_dist(A)\n")
        for (resid, resname), count in sorted_contacts:
            occ = count / total_frames * 100
            dists = per_frame_mindist[(resid, resname)]
            avg_d = np.mean(dists)
            std_d = np.std(dists)
            f.write(f"{resid},{resname},{count},{occ:.1f},{avg_d:.2f},{std_d:.2f}\n")

    print(f"\n  -> {occ_file}")

    # 상위 20개 콘솔 출력
    print(f"\n{'='*70}")
    print(f"  리간드 {cutoff}A 이내 잔기 접촉 빈도 (상위 20)")
    print(f"{'='*70}")
    print(f"  {'Resid':>6} {'Resname':>8} {'Occupancy':>10} {'Avg Dist':>10} {'Std':>6}")
    print(f"  {'-'*44}")
    for (resid, resname), count in sorted_contacts[:20]:
        occ = count / total_frames * 100
        dists = per_frame_mindist[(resid, resname)]
        print(f"  {resid:>6} {resname:>8} {occ:>9.1f}% {np.mean(dists):>9.2f} {np.std(dists):>6.2f}")

    # --- 3.2 잔기별 최소 거리 시계열 (상위 10 잔기) ---
    top10 = sorted_contacts[:10]

    dist_file = os.path.join(outdir, f"{prefix}_mindist_timeseries.csv")
    with open(dist_file, 'w') as f:
        header = "time(ns)," + ",".join([f"{rn}{ri}" for (ri, rn), _ in top10])
        f.write(header + "\n")
        for i, t in enumerate(times):
            vals = [f"{per_frame_mindist[(ri,rn)][i]:.3f}" for (ri, rn), _ in top10]
            f.write(f"{t:.3f}," + ",".join(vals) + "\n")

    print(f"  -> {dist_file}")

    # --- 3.3 수소결합 요약 ---
    hb_file = os.path.join(outdir, f"{prefix}_hbond_summary.csv")

    all_hbonds = []

    if hbonds_pl.results.hbonds is not None and len(hbonds_pl.results.hbonds) > 0:
        for hb in hbonds_pl.results.hbonds:
            frame, donor_idx, hydrogen_idx, acceptor_idx, dist, angle = hb
            donor = u.atoms[int(donor_idx)]
            acceptor = u.atoms[int(acceptor_idx)]
            all_hbonds.append({
                'frame': int(frame),
                'donor': f"{donor.resname}{donor.resid}:{donor.name}",
                'acceptor': f"{acceptor.resname}{acceptor.resid}:{acceptor.name}",
                'dist': dist,
                'angle': angle,
                'type': 'Prot->Lig'
            })

    if hbonds_lp.results.hbonds is not None and len(hbonds_lp.results.hbonds) > 0:
        for hb in hbonds_lp.results.hbonds:
            frame, donor_idx, hydrogen_idx, acceptor_idx, dist, angle = hb
            donor = u.atoms[int(donor_idx)]
            acceptor = u.atoms[int(acceptor_idx)]
            all_hbonds.append({
                'frame': int(frame),
                'donor': f"{donor.resname}{donor.resid}:{donor.name}",
                'acceptor': f"{acceptor.resname}{acceptor.resid}:{acceptor.name}",
                'dist': dist,
                'angle': angle,
                'type': 'Lig->Prot'
            })

    # 수소결합 쌍별 빈도
    hb_pairs = {}
    for hb in all_hbonds:
        key = (hb['donor'], hb['acceptor'], hb['type'])
        if key not in hb_pairs:
            hb_pairs[key] = {'count': 0, 'dists': [], 'angles': []}
        hb_pairs[key]['count'] += 1
        hb_pairs[key]['dists'].append(hb['dist'])
        hb_pairs[key]['angles'].append(hb['angle'])

    with open(hb_file, 'w') as f:
        f.write("donor,acceptor,type,count,occupancy(%),avg_dist(A),avg_angle(deg)\n")
        for (donor, acceptor, htype), data in sorted(hb_pairs.items(), key=lambda x: -x[1]['count']):
            occ = data['count'] / total_frames * 100
            f.write(f"{donor},{acceptor},{htype},{data['count']},{occ:.1f},"
                    f"{np.mean(data['dists']):.2f},{np.mean(data['angles']):.1f}\n")

    print(f"  -> {hb_file}")

    print(f"\n{'='*70}")
    print(f"  수소결합 요약 (상위 15)")
    print(f"{'='*70}")
    if hb_pairs:
        print(f"  {'Donor':>25} {'Acceptor':>25} {'Type':>10} {'Occ%':>6} {'Dist':>6} {'Angle':>6}")
        print(f"  {'-'*82}")
        for (donor, acceptor, htype), data in sorted(hb_pairs.items(), key=lambda x: -x[1]['count'])[:15]:
            occ = data['count'] / total_frames * 100
            print(f"  {donor:>25} {acceptor:>25} {htype:>10} {occ:>5.1f} "
                  f"{np.mean(data['dists']):>6.2f} {np.mean(data['angles']):>6.1f}")
    else:
        print("  수소결합 감지 안됨")

    # --- 3.4 그래프 ---
    if HAS_PLT:
        _plot_results(outdir, prefix, sorted_contacts, total_frames,
                      per_frame_mindist, times, all_hbonds, cutoff)


def _plot_results(outdir, prefix, sorted_contacts, total_frames,
                  per_frame_mindist, times, all_hbonds, cutoff):
    """결과 그래프 생성"""

    # Occupancy bar chart
    fig, ax = plt.subplots(figsize=(12, 5))
    top20 = sorted_contacts[:20]
    labels = [f"{rn}{ri}" for (ri, rn), _ in top20]
    occs = [c / total_frames * 100 for _, c in top20]
    colors = ['#c0392b' if o > 80 else '#e67e22' if o > 50 else '#3498db' for o in occs]
    ax.barh(range(len(labels)), occs, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Occupancy (%)")
    ax.set_title(f"{prefix} - Ligand Contact Occupancy ({cutoff}A cutoff)")
    ax.invert_yaxis()
    ax.axvline(x=80, color='red', linestyle='--', alpha=0.3)
    ax.axvline(x=50, color='orange', linestyle='--', alpha=0.3)
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{prefix}_contact_occupancy.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> {outpath}")

    # Distance timeseries (상위 5)
    fig, ax = plt.subplots(figsize=(12, 5))
    top5 = sorted_contacts[:5]
    colors_line = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for idx, ((ri, rn), _) in enumerate(top5):
        dists = per_frame_mindist[(ri, rn)]
        ax.plot(times, dists, label=f"{rn}{ri}",
                alpha=0.8, linewidth=0.6, color=colors_line[idx])
    ax.axhline(y=cutoff, color='red', linestyle='--', alpha=0.5, label=f'{cutoff}A cutoff')
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Min Distance to Ligand (A)")
    ax.set_title(f"{prefix} - Residue-Ligand Minimum Distance")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{prefix}_mindist_timeseries.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> {outpath}")

    # H-bond count timeseries
    if all_hbonds:
        hb_per_frame = {}
        for hb in all_hbonds:
            hb_per_frame[hb['frame']] = hb_per_frame.get(hb['frame'], 0) + 1
        counts = [hb_per_frame.get(f, 0) for f in range(len(times))]

        fig, ax = plt.subplots(figsize=(12, 3))
        ax.fill_between(times, counts, alpha=0.3, color='darkgreen')
        ax.plot(times, counts, color='darkgreen', linewidth=0.5)
        # running average
        if len(counts) > 10:
            dt = times[1] - times[0] if len(times) > 1 else 1
            window = max(1, int(1.0 / dt))
            if len(counts) > window:
                avg = np.convolve(counts, np.ones(window)/window, mode='valid')
                ax.plot(times[:len(avg)], avg, color='black', linewidth=1.5,
                        alpha=0.8, label='1 ns avg')
                ax.legend(fontsize=8)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("H-bond Count")
        ax.set_title(f"{prefix} - Protein-Ligand Hydrogen Bonds")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        outpath = os.path.join(outdir, f"{prefix}_hbond_timeseries.png")
        plt.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  -> {outpath}")


# ============================================================
# main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Protein-Ligand 접촉 분석 (MDAnalysis)")
    parser.add_argument("--prefix", default=None,
                        help="파일 prefix (미지정 시 자동 감지)")
    parser.add_argument("-s", "--tpr", default=None, help="TPR 파일 (미지정 시 prefix.tpr)")
    parser.add_argument("-f", "--xtc", default=None,
                        help="trajectory 파일 (미지정 시 noPBC.xtc > .xtc 순서로 자동 선택)")
    parser.add_argument("-l", "--ligand", default=None,
                        help="리간드 선택 문자열 (미지정 시 자동 감지, 예: 'resname UNL')")
    parser.add_argument("-c", "--cutoff", type=float, default=5.0,
                        help="접촉 거리 cutoff in A (default: 5.0)")
    parser.add_argument("-o", "--outdir", default=None,
                        help="출력 디렉토리 (default: {prefix}_contacts)")
    parser.add_argument("--start", type=float, default=None,
                        help="분석 시작 시간 (ns)")
    parser.add_argument("--stop", type=float, default=None,
                        help="분석 종료 시간 (ns)")

    args = parser.parse_args()

    # --- prefix 자동 감지 ---
    if args.prefix is None:
        args.prefix = find_prefix()
        if args.prefix is None:
            print("[ERROR] .tpr + .xtc 파일을 찾을 수 없습니다. --prefix를 지정하세요.")
            exit(1)
        print(f"  자동 감지 prefix: {args.prefix}")

    PREFIX = args.prefix

    # --- 파일 경로 결정 ---
    tpr = args.tpr or f"{PREFIX}.tpr"
    xtc = args.xtc or find_trajectory(PREFIX)
    outdir = args.outdir or f"{PREFIX}_contacts"

    if not os.path.exists(tpr):
        print(f"[ERROR] TPR 파일 없음: {tpr}")
        exit(1)
    if xtc is None or not os.path.exists(xtc):
        print(f"[ERROR] trajectory 파일 없음")
        exit(1)

    print("\n" + "=" * 60)
    print("  Protein-Ligand Contact Analysis")
    print("=" * 60)
    print(f"  Prefix  : {PREFIX}")
    print(f"  TPR     : {tpr}")
    print(f"  XTC     : {xtc}")
    print(f"  Cutoff  : {args.cutoff} A")
    print(f"  Output  : {outdir}")

    # --- Universe 로드 ---
    print("\n[ trajectory 로딩 ]")
    u = mda.Universe(tpr, xtc)
    total_time_ns = len(u.trajectory) * u.trajectory.dt / 1000.0
    print(f"  Atoms   : {len(u.atoms)}")
    print(f"  Frames  : {len(u.trajectory)}")
    print(f"  dt      : {u.trajectory.dt:.1f} ps")
    print(f"  Total   : {total_time_ns:.1f} ns")

    # --- 리간드 자동 감지 ---
    if args.ligand is None:
        detected = detect_ligand_resname(u)
        if detected:
            args.ligand = f"resname {detected}"
            print(f"  Ligand  : {args.ligand} (자동 감지)")
        else:
            print("[ERROR] 리간드를 자동 감지할 수 없습니다. -l 옵션으로 지정하세요.")
            print(f"  사용 가능한 resname: {sorted(set(u.atoms.resnames))}")
            exit(1)
    else:
        print(f"  Ligand  : {args.ligand}")

    if args.start is not None or args.stop is not None:
        s = args.start or 0
        e = args.stop or total_time_ns
        print(f"  구간    : {s:.1f} - {e:.1f} ns")

    print("=" * 60)

    # --- 접촉 분석 ---
    print("\n[ 접촉 분석 ]")
    contact_counts, total_frames, per_frame_mindist, times = analyze_contacts(
        u, args.ligand, args.cutoff, args.start, args.stop
    )

    # --- 수소결합 분석 ---
    print("\n[ 수소결합 분석 ]")
    hbonds_pl, hbonds_lp = analyze_hbonds(u, args.ligand, args.start, args.stop)

    # --- 결과 출력 ---
    print("\n[ 결과 출력 ]")
    write_results(outdir, PREFIX, contact_counts, total_frames,
                  per_frame_mindist, times, hbonds_pl, hbonds_lp,
                  u, args.ligand, args.cutoff)

    print(f"\n{'='*60}")
    print(f"  완료! 결과: {outdir}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
