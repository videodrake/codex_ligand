# PyMOL 시각화 가이드 — EGFR-MYO1D 최종 결과

## 필요한 파일 (HPC에서 로컬로 복사)

```bash
# 구조 파일
scp hpc:codex_ligand2/input/receptors/3GT8_raw.pdb .
scp hpc:codex_ligand2/input/receptors/EGFR_170-200.pdb .

# PPI 도킹 대표 모델 (top-ranked)
scp hpc:codex_ligand2/output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/docking_3GT8_raw_ext_beta_meander/final_result/Rank01*.pdb .
scp hpc:codex_ligand2/output/workflow_a/phase2_ppi_docking/EGFR_170-200/prod_seed0/docking_EGFR_170-200_ext_beta_meander/final_result/Rank01*.pdb .

# Focused Vina 도킹 포즈 (PKT07, PKT34)
scp -r hpc:codex_ligand2/output/workflow_b/phase3_focused_docking/runs/3GT8_raw/3GT8_raw_PKT07/ .
scp -r hpc:codex_ligand2/output/workflow_b/phase3_focused_docking/runs/EGFR_170-200/EGFR_170-200_PKT34/ .

# fpocket 포켓 구조
scp -r hpc:codex_ligand2/output/workflow_b/phase2_pocket_analysis/3GT8_raw/3GT8_raw_out/ .
scp -r hpc:codex_ligand2/output/workflow_b/phase2_pocket_analysis/EGFR_170-200/EGFR_170-200_out/ .
```

---

## Figure 1: PPI 결합 부위 전체 조감도

**목적**: EGFR C-lobe 위에 MYO1D가 어떻게 붙는지 보여주기

```python
# PyMOL script
load Rank01*3GT8*.pdb, ppi_complex

# Chain 색상
color lightblue, chain A      # EGFR
color salmon, chain B          # MYO1D

# 표면 표시
show surface, chain A
set transparency, 0.7, chain A
show cartoon, chain B

# PPI hotspot 잔기 강조 (EGFR 측, 3/3 상태 공통)
select egfr_hotspots, chain A and resi 941+977+993+986+940+980+992+937+982
color red, egfr_hotspots
show sticks, egfr_hotspots

# MYO1D active face 강조 (sheet 8/9)
select myo1d_active, chain B and resi 961+962+963+964+968+969+970+971+972
color orange, myo1d_active
show sticks, myo1d_active

# 라벨
label egfr_hotspots and name CA, "%s%s" % (resn, resi)

# 시점
orient
zoom ppi_complex

# 저장
png figure1_ppi_overview.png, width=2400, height=1800, dpi=300, ray=1
```

---

## Figure 2: PKT34 — PPI 근접 allosteric pocket (핵심 발견)

**목적**: WF-B 핵심 발견인 PKT34가 PPI에서 9.4Å 거리에 있음을 보여주기

```python
# PyMOL script
load EGFR_170-200.pdb, receptor

# 수용체 표면
show surface, receptor
color white, receptor
set transparency, 0.5

# PPI hotspot 잔기 (EGFR 측)
select ppi_patch, resi 941+977+993+986+940+980+992+937+982
color red, ppi_patch

# PKT34 포켓 잔기 (62개)
select pkt34, resi 718+720+721+722+723+726+728+743+745+748+791+792+793+796+797+799+800+803+835+836+837+838+839+841+842+844+854+855+858+859+864+865+866+867+868+869+870+871+872+873+874+875+876+877+878+879+880+881+882+885+886+889+891+895+896+899+906+913+914+920+921+924
color cyan, pkt34
show sticks, pkt34

# PKT34 centroid 표시 (53.7, 50.4, 37.4)
pseudoatom pkt34_center, pos=[53.683, 50.357, 37.367]
show spheres, pkt34_center
color cyan, pkt34_center
set sphere_scale, 1.5, pkt34_center

# PPI patch centroid (52.8, 49.2, 28.1) — EGFR_170-200 값
pseudoatom ppi_center, pos=[52.8, 49.2, 28.1]
show spheres, ppi_center
color red, ppi_center
set sphere_scale, 1.5, ppi_center

# 거리 표시 (9.4Å)
distance pkt34_ppi_dist, pkt34_center, ppi_center
set dash_color, yellow, pkt34_ppi_dist
set dash_width, 3
set label_size, 20

# Vina 도킹 포즈 (있으면)
# load EGFR_170-200_PKT34/173940/EGFR_170-200_PKT34_173940_seed1.pdbqt, pose_173940
# show sticks, pose_173940
# color green, pose_173940

orient
zoom pkt34, 15
png figure2_pkt34_allosteric.png, width=2400, height=1800, dpi=300, ray=1
```

---

## Figure 3: PKT07 — PPI rim pocket (최고 스코어)

**목적**: PKT07이 PPI 가장자리에 위치함을 보여주기

```python
# PyMOL script
load 3GT8_raw.pdb, receptor

show surface, receptor
color white, receptor
set transparency, 0.5

# PPI hotspot 잔기
select ppi_patch, resi 941+977+993+986+940+980+992+937+982
color red, ppi_patch

# PKT07 포켓 잔기 (42개)
select pkt07, resi 716+718+719+720+721+722+723+724+726+728+737+738+739+740+741+743+744+745+766+775+777+788+790+791+792+793+794+796+797+800+837+841+842+844+846+854+855+856+858+877+1001+1002
color blue, pkt07
show sticks, pkt07

# PKT07 centroid (-16.3, -7.6, -66.3)
pseudoatom pkt07_center, pos=[-16.316, -7.557, -66.272]
show spheres, pkt07_center
color blue, pkt07_center
set sphere_scale, 1.5, pkt07_center

# PPI patch centroid (3GT8_raw: -2.2, 1.1, -57.6)
pseudoatom ppi_center, pos=[-2.2, 1.1, -57.6]
show spheres, ppi_center
color red, ppi_center
set sphere_scale, 1.5, ppi_center

# 거리 표시 (18.7Å)
distance pkt07_ppi_dist, pkt07_center, ppi_center
set dash_color, yellow, pkt07_ppi_dist

orient
zoom pkt07, 15
png figure3_pkt07_rim.png, width=2400, height=1800, dpi=300, ray=1
```

---

## Figure 4: WF-A vs WF-B 비교 — ATP 포켓 vs 신규 포켓

**목적**: blind Vina(WF-A)로 찾은 ATP 포켓 vs PPI-first(WF-B)로 찾은 PKT34 비교

```python
# PyMOL script (EGFR_170-200 기준)
load EGFR_170-200.pdb, receptor

show surface, receptor
color white, receptor
set transparency, 0.6

# PPI 패치 (빨간색)
select ppi_patch, resi 941+977+993+986+940+980+992+937+982
color red, ppi_patch

# WF-A ATP 포켓 P004 — PPI 14.5Å but ATP (배제됨)
# ATP 포켓 위치 (대략적 — 실제 centroid는 valid_sites.csv에서 확인)
# 주요 ATP 잔기: CYS797, ASP855, LYS745 등
select atp_region, resi 745+766+790+791+792+793+794+795+796+797+854+855+856+858
color orange, atp_region
show sticks, atp_region

# WF-B PKT34 — PPI 9.4Å, allosteric (신규 발견)
select pkt34, resi 718+720+721+722+723+726+728+743+745+791+792+793+796+797+799+800+835+836+837+841+842+844+854+855+858+859+864+865+866+867+868+869+870+871+872+873+874+875+876+877+878+879+880+881+882+885+886+889+891+895+896+899+906+913+914+920+921+924
color cyan, pkt34

orient
png figure4_atp_vs_pkt34.png, width=2400, height=1800, dpi=300, ray=1
```

---

## Figure 5: Orientation Filter 시각화

**목적**: Orientation filter가 올바른 방향(active face → receptor)의 모델만 선택했음을 보여주기

```python
# Pass 모델 vs Fail 모델 비교
# Rank PDB 2개를 로드 (orientation_filter_log.csv에서 pass/fail 모델 ID 확인)

# 예시:
load pass_model.pdb, pass_orientation
load fail_model.pdb, fail_orientation

# MYO1D active face 표시
select active_pass, pass_orientation and chain B and resi 961-964+968-972
select active_fail, fail_orientation and chain B and resi 961-964+968-972

color green, active_pass    # 올바른 방향: active face가 receptor 쪽
color red, active_fail      # 잘못된 방향: active face가 반대쪽

show spheres, active_pass
show spheres, active_fail
```

---

## 확인해야 할 사항 (시각적 검증)

### 필수 확인

1. **PKT34가 PPI 인터페이스와 공간적으로 가까운지**
   - PKT34 centroid (53.7, 50.4, 37.4)와 PPI patch centroid (52.8, 49.2, 28.1)
   - 9.4Å 거리가 시각적으로 확인되는지
   - PKT34가 ATP 포켓과 겹치지 않는지 (allosteric이어야 함)

2. **PKT07이 PPI rim에 있는지**
   - PPI 패치 가장자리에 포켓이 위치하는지
   - Hotspot 잔기(ILE941 등)와의 공간적 관계

3. **PPI hotspot의 공간적 패치 연속성**
   - ILE941, ARG977, THR993, ARG986가 하나의 연속된 표면 패치인지
   - 또는 분산된 잔기인지 (연속이면 druggable target으로서 가치 높음)

4. **3개 상태 구조 중첩**
   - 3GT8_raw, EGFR_160-185, EGFR_170-200를 align
   - PPI 패치 위치가 상태 간 보존되는지
   - PKT07/PKT34에 해당하는 포켓이 각 상태에서 유사한 위치인지

### 선택 확인

5. **Vina 도킹 포즈 품질**
   - PKT34 내부에 리간드가 잘 들어가는지
   - 수소결합/소수성 상호작용 패턴

6. **Orientation filter pass vs fail 비교**
   - Active face가 receptor를 향하는지 (pass) vs 반대쪽인지 (fail)

---

## 3상태 중첩 명령어

```python
# 3상태 중첩
load 3GT8_raw.pdb, state_3gt8
load EGFR_160-185.pdb, state_160
load EGFR_170-200.pdb, state_170

# Kinase domain으로 정렬 (잔기 699-1007)
align state_160 and resi 699-1007, state_3gt8 and resi 699-1007
align state_170 and resi 699-1007, state_3gt8 and resi 699-1007

# 색상
color lightblue, state_3gt8
color lightorange, state_160
color lightgreen, state_170

# PPI hotspot 공통 잔기 (3/3 상태)
select hotspots, resi 941+977+993+986+940+980+992+937+982
show sticks, hotspots
color red, hotspots

orient
png figure_3state_overlay.png, width=2400, height=1800, dpi=300, ray=1
```
