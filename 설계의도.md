# Workflow A: 설계 의도 중심 상세 보고

**EGFR-MYO1D Standard Production Pipeline**
Lab Meeting | March 2026

---

## 요약: 왜 이렇게 만들었는가

이 파이프라인은 "EGFR 표면에서 MYO1D 부착을 약물로 교란할 수 있는 포켓은 어디인가"라는 질문에 답하기 위해 설계되었다. 이 질문에는 두 가지 독립적인 하위 질문이 포함되어 있고, Workflow A는 이 둘을 병렬로 답한 뒤 사후적으로 수렴을 평가하는 구조를 취한다.

**Vina Branch (Step 1-6):** "EGFR 표면 어디에 소분자가 잘 붙는가?"를 묻는다. Blind docking으로 편향 없이 전체 표면을 탐색하고(Step 1), raw 결과를 구조화된 테이블로 변환하고(Step 2), 좌표 공간의 데이터를 잔기 공간으로 번역하고(Step 3), 산발적 포즈를 재현 가능한 포켓으로 그룹화하고(Step 4), 포켓의 특성을 다양한 관점에서 요약하고(Step 5), 3개 receptor state에 걸친 포켓의 구조적 안정성을 평가한다(Step 6).

**PPI Branch (Step 7-8):** "EGFR 표면 어디에 MYO1D가 붙는가?"를 묻는다. PyRosetta global blind PPI docking으로 300K개 도킹 가설을 생성하고(Step 7), 이를 모델 수준이 아닌 잔기 수준의 receptor-side patch 증거로 응축한 뒤, 3개 state에 걸친 patch robustness를 평가하고, LightDock으로 독립 검증한다(Step 8).

**통합 (Step 9-10):** 두 branch의 증거가 공간적으로 겹치는지를 3축(Vina 증거 50점 + PPI 공간 관계 20점 + Cross-Receptor 일관성 30점) 체계로 정량화하여 STRONG/MODERATE/WEAK 판정을 내리고(Step 9), 결과를 사람이 읽을 수 있는 보고서로 정리하고 출력 무결성을 검증한다(Step 10).

전체를 관통하는 설계 원칙은 다음과 같다.

**"편향 없이 시작하라."** Vina도 PPI도 blind docking으로 시작한다. 사전 지식으로 탐색 범위를 제한하지 않는다.

**"교차 검증 구조를 내장하라."** 다중 리간드(cross-chemical consensus), 다중 receptor state(cross-state robustness), 다중 방법론(PyRosetta + LightDock)이 각각 독립적인 검증 축을 형성한다. 이 교차 검증은 사후 분석이 아니라 Step 1의 입력 설계 단계에서부터 의도된 것이다.

**"모든 중간 결과를 machine-readable하게 보존하라."** 파이프라인의 모든 단계가 CSV를 생성하고, 모든 결정(클러스터링 병합, 필터 탈락)의 근거가 로그로 남는다. 이를 통해 파라미터 sensitivity 분석이나 post-hoc 재분석을 도킹을 다시 돌리지 않고 할 수 있다.

**"결론이 아니라 증거를 남겨라."** Verdict의 STRONG 판정은 "이 포켓이 실제로 MYO1D를 교란한다"는 뜻이 아니라, "이 포켓에 대한 계산적 증거가 여러 축에서 수렴한다"는 뜻이다. 모든 판정에는 축별 점수 분해가 제공되어, 왜 그 판정이 나왔는지를 투명하게 추적할 수 있다.

---

## Step 1: Vina Blind Docking — 설계 의도

**해결하려는 문제:**

우리 프로젝트의 질문 중 하나가 "EGFR 표면 어디에 소분자가 잘 붙는가"이다. 이걸 답하는 전통적인 방법은 도킹인데, 도킹에는 두 가지 접근이 있다. **Focused docking**은 "여기가 포켓일 것 같다"는 사전 지식을 넣고 그 부위만 집중 탐색하는 것이고, **blind docking**은 사전 지식 없이 단백질 전체 표면을 열어두는 것이다.

우리가 blind를 선택한 이유가 명확하다. 이 프로젝트는 **알려진 포켓을 확인하는 게 아니라, 포켓 후보를 발견하는 단계**이기 때문이다. Focused로 시작하면 "어디를 볼 것인가"라는 편향이 처음부터 들어간다. Blind로 시작해야 나중에 "왜 이 포켓이 나왔는가"를 데이터로 설명할 수 있다.

**exhaustiveness를 384로 올린 이유:**

Vina의 탐색 알고리즘은 Iterated Local Search이다. 랜덤 perturbation → BFGS local optimization → Metropolis acceptance를 반복하는데, exhaustiveness는 이 독립 탐색을 몇 회 병렬로 돌리냐는 파라미터이다. Vina 소프트웨어 기본값은 8인데, 이건 "결합 부위를 이미 알고, 20-30Å 정도의 작은 box에서 탐색할 때" 기준이다.

우리 search box는 EGFR kinase domain 전체를 감싸는 70Å 이상이다. 탐색 공간이 기본 가정보다 수십 배 넓으니까, 같은 exhaustiveness로는 표면의 넓은 영역이 제대로 샘플링되지 않는다. 특히 얕거나 작은 포켓은 에너지 landscape에서 좁은 basin을 형성하는데, 탐색 횟수가 부족하면 이런 포켓을 아예 못 찾고 지나칠 수 있다. 384는 기본값의 48배로, blind docking 문헌에서 권장하는 높은 수준의 탐색 강도이다.

**리간드를 3종 쓴 설계 의도:**

이건 단순히 "더 많이 하면 좋겠다"가 아니라, **cross-chemical consensus**라는 검증 논리를 파이프라인에 내장하기 위해서이다.

한 리간드가 어떤 포켓에 들어갔다고 해서 그 포켓이 druggable하다고 바로 말할 수 없다. 그 리간드의 화학적 성질(전하, 크기, 소수성)이 우연히 그 포켓에 맞았을 수 있기 때문이다. 하지만 화학적으로 다른 3종 리간드가 **독립적으로** 같은 포켓에 수렴했다면, 그 포켓의 druggability는 특정 화학 구조에 의존하지 않는 일반적인 성질일 가능성이 높다. 이 정보는 나중에 Step 5의 Drug-Pocket Map에서 `n_ligand` 컬럼으로 정량화된다.

**Receptor state를 3개 쓴 설계 의도:**

이것도 파이프라인 전체의 핵심 설계와 직결된다. 이 프로젝트는 "어디에 포켓이 있는가"만 묻는 게 아니라, **"그 포켓이 receptor state가 변해도 유지되는가"**를 묻고 있다.

X-ray 결정구조(3GT8_raw) 하나만 쓰면, 결정화 조건에서 고정된 단일 conformation만 보게 된다. 단백질은 용액에서 움직이고, C-lobe 표면의 loop나 helix 배향이 달라지면 포켓이 열리거나 닫힐 수 있다. MD 시뮬레이션에서 뽑은 2개 클러스터 대표 구조(EGFR_160-185, EGFR_170-200)를 추가하면, 동적 구조 변이(conformational heterogeneity) 속에서도 일관되게 존재하는 포켓을 골라낼 수 있다.

이것이 나중에 Step 6에서 "3개 state 모두에서 같은 위치에 포켓이 나타났는가"를 비교하는 cross-receptor 분석으로 직접 이어진다. 결국 Step 1에서 3 receptor × 3 ligand = 9 조합으로 설계한 것 자체가, Step 6의 교차 비교를 가능하게 하기 위한 사전 투자인 셈이다.

**n_poses를 500으로 설정한 이유:**

Blind docking에서 포즈는 표면 전체에 흩어진다. Focused docking이면 20-30개 pose로 한 포켓 안의 미세한 배향 차이를 충분히 볼 수 있지만, blind에서는 포즈들이 여러 포켓에 분산되기 때문에, 한 포켓에 통계적으로 의미 있는 수의 포즈가 모이려면 총 포즈 수가 충분해야 한다.

500이면 한 조합에서 major pocket 3-5개에 각각 수십 개씩 포즈가 배정될 수 있고, 이 정도면 Step 4 클러스터링에서 "이 포켓은 pose가 많이 모인 재현 가능한 포켓이다 vs 이건 한두 개 포즈만 들어간 노이즈다"를 통계적으로 구분할 수 있다.

**Deterministic seed를 쓴 이유:**

이건 재현성(reproducibility) 원칙이다. Vina의 Monte Carlo 탐색은 본질적으로 확률적인데, 시드를 고정하지 않으면 같은 입력으로 돌려도 매번 다른 결과가 나온다. 연구에서 "이 포켓은 안정적으로 발견된다"고 주장하려면, 최소한 같은 조건에서 같은 결과가 나온다는 걸 보장해야 한다.

`derive_docking_seed(base, receptor_id, ligand_id)` 함수가 receptor와 ligand ID의 조합으로 시드를 결정론적으로 생성하기 때문에, 어떤 조합이든 다시 돌리면 동일한 결과를 재현할 수 있다. 동시에 조합마다 시드가 다르기 때문에 조합 간에는 독립적인 탐색이 보장된다.

**정리:**

Step 1의 설계 의도는 "편향 없이, 충분한 밀도로, 재현 가능하게, 교차 비교가 가능한 형태로 raw docking 데이터를 생성하라"이다. blind docking, 높은 exhaustiveness, 다중 리간드, 다중 receptor state, deterministic seed — 이 모든 선택이 각각 이유가 있고, 그 이유가 파이프라인의 다운스트림 단계(특히 Step 5의 consensus, Step 6의 cross-receptor 비교)에서 활용된다.

---

## Step 2: Pose Parsing — 설계 의도

**해결하려는 문제:**

이 파이프라인의 핵심 설계 원칙 중 하나가 **"모든 중간 결과를 machine-readable하게 유지한다"**이다. Step 2는 이 원칙이 처음 구현되는 지점이다.

Vina 출력 PDBQT는 본질적으로 **구조 파일**이지 **분석 파일**이 아니다. 에너지 값이 REMARK 줄에 텍스트로 묻혀있고, 좌표는 원자 단위로 흩어져 있다. 이 상태로는 "3GT8에서 가장 에너지가 좋은 포즈 10개의 centroid를 비교하라" 같은 질문에 답하려면 매번 스크립트를 새로 짜야 한다.

일반적인 도킹 연구에서는 이 단계를 대충 넘기는 경우가 많다. 결과 PDBQT를 PyMOL에서 열어보고, 에너지 좋은 것 몇 개 골라서 수동으로 분석하는 식이다. 단일 receptor, 단일 ligand면 그래도 되는데, 우리는 **3 receptor × 3 ligand = 9 조합, 4,500 포즈**를 다루고 있고, 이걸 나중에 cross-receptor 비교까지 해야 한다. 수동 접근은 불가능하다.

**설계 결정:**

그래서 PDBQT를 **단일 CSV 테이블로 평탄화(flatten)**하는 걸 파이프라인의 첫 후처리 단계로 넣었다. 각 포즈가 테이블의 한 행이 되고, 모든 메타데이터(어떤 receptor, 어떤 ligand, 어떤 세팅)가 컬럼으로 들어간다.

이렇게 하면 이후 모든 분석이 이 테이블 위에서 **컬럼 추가** 방식으로 진행된다. Step 3에서 contact_residues 컬럼이 붙고, Step 4에서 pocket_id 컬럼이 붙고, 테이블이 점점 풍부해지는 구조이다. 데이터가 한 번도 분절되지 않고 하나의 테이블에서 계속 축적되니까, 나중에 "포켓 X에 속한 포즈들의 평균 에너지"나 "receptor A vs B에서 같은 포켓에 들어간 리간드 비교" 같은 질문에 필터 한 줄로 답할 수 있다.

**centroid를 왜 여기서 미리 계산하냐면:**

이것도 의도적인 결정인데, 원자 좌표 수백 개를 매 분석마다 다시 읽는 것보다, 한 번 centroid로 요약해서 테이블에 넣어두는 게 효율적이다. centroid는 "이 포즈가 3D 공간 어디에 있는가"를 점 하나로 대표하는데, Step 4 클러스터링에서 바로 이 좌표를 거리 계산에 쓴다. 원본 원자 좌표가 필요한 경우를 위해 raw_pose_file 경로도 함께 기록해두었다.

**coverage 파일을 왜 따로 만드냐면:**

`vina_postprocess_coverage.csv`는 "요청한 500개 pose 중 실제로 몇 개가 파싱됐는가"를 기록한다. 이것도 설계 원칙의 일부인데, 이 파이프라인은 **silent failure를 허용하지 않는다**. Vina가 중간에 죽었거나 파일이 잘렸으면, coverage가 100% 미만으로 찍혀서 다운스트림에서 잡힌다. 이게 없으면 "왜 이 receptor-ligand 조합만 포즈가 적지?"를 한참 뒤에야 알게 된다.

**정리:**

Step 2의 설계 의도는 "분석 가능한 구조로 만들어라, 그리고 그 과정에서 아무것도 잃지 마라"이다. 화려한 과학이 아니라 **데이터 엔지니어링**이지만, 이 테이블이 파이프라인 전체의 backbone이기 때문에 여기서 구조를 잘못 잡으면 뒤가 다 무너진다.

---

## Step 3: 접촉 잔기 추출 — 설계 의도

**해결하려는 문제:**

Step 2까지 끝나면 4,500개 포즈가 테이블에 정리되어 있는데, 각 포즈에 대해 알 수 있는 건 **에너지(affinity)와 위치(centroid)** 뿐이다. 문제는 이 정보만으로는 생물학적 해석이 불가능하다는 것이다.

"이 포즈의 affinity가 -7.2 kcal/mol이고 centroid가 (-38, 2, -74)에 있다" — 구조생물학에서 의미 있는 단위는 좌표가 아니라 **잔기(residue)**이다. "A-loop의 PRO848 근처에 붙는다"가 해석 가능한 정보이다.

그래서 이 단계의 설계 의도는 **포즈를 좌표 공간에서 잔기 공간으로 변환**하는 것이다. 솔직히 이건 도킹 후처리에서 당연히 하는 일이다. 포즈가 어떤 잔기와 접촉하는지 안 보는 도킹 분석은 없다. 다만 이 파이프라인에서는 4,500개 포즈 전부에 대해 일괄로, 일관된 기준으로, machine-readable한 형태로 이걸 해야 한다는 점이 수동 분석과 다르다.

**왜 거리 기반 접촉 판정을 쓰냐면:**

이 단계에서는 의도적으로 **가장 단순한 방법, 즉 거리 cutoff**를 선택했다. 이 단계의 목적은 "이 잔기가 리간드와 특정 종류의 상호작용을 하는가"가 아니라, **"이 잔기가 리간드 근처에 물리적으로 존재하는가"**이다. 4,500개 포즈에 대해 일괄 적용해야 하고, 결과가 다운스트림의 클러스터링 병합(Step 4)과 포켓 요약(Step 5)으로 직접 흘러가야 하니까, 복잡한 상호작용 분류보다는 robust하고 일관된 기하학적 판정이 적합하다.

**cutoff를 4.0Å로 잡은 근거:**

4.0Å은 비공유 상호작용의 일반적인 접촉 거리이다. van der Waals 반경의 합(~3.4-3.8Å)에 약간의 여유를 둔 수준이고, 이보다 짧으면 수소결합이나 직접 접촉, 이보다 길면 2차 접촉권(second shell)으로 들어간다. PPI 분석에서는 보통 8Å Cα distance를 쓰지만, 소분자-수용체 접촉은 원자 수준이 더 적절하기 때문에 4.0Å이 표준적인 선택이다. 이 값은 config에서 조정 가능하게 설계했다.

**최소 거리를 기록하는 이유:**

각 접촉 잔기에 대해 단순히 "접촉이다/아니다"만 기록하는 게 아니라, **min_distance_A**(그 잔기의 가장 가까운 원자까지의 최소 거리)를 함께 저장한다. 지금 당장은 4.0Å cutoff로 binary 판정만 해도 되지만, 나중에 distance-weighted 분석이나 cutoff 재조정을 **원본 데이터를 다시 계산하지 않고** 할 수 있다. 이것도 파이프라인의 설계 원칙 중 하나인 "raw evidence를 보존하라"의 구현이다.

**이 단계가 파이프라인 전체에서 왜 중요하냐면:**

이 시점을 기준으로 데이터의 성격이 바뀐다. Step 2까지는 **물리화학적 데이터**(에너지, 좌표)였는데, Step 3 이후부터는 **구조생물학적 데이터**(잔기)가 된다. 이 전환이 없으면 뒤의 분석 — 포켓 병합(잔기 Jaccard로 병합 판정), 포켓 요약(top_residues), cross-receptor 비교(잔기 중첩) — 이 불가능하다.

그리고 결과적으로, 이 잔기 데이터가 Step 9 Verdict에서 PPI branch의 interface 잔기와 비교되는 데도 쓰인다. 처음부터 Verdict를 위해 설계한 건 아니지만, 잔기라는 공통 단위로 데이터를 변환해두었기 때문에 나중에 두 branch를 같은 기준으로 비교할 수 있게 된 것이다.

---

## Step 4: Pocket Clustering — 설계 의도

**해결하려는 문제:**

Step 3까지 끝나면 4,500개 포즈 각각이 "어디에 있고, 어떤 잔기와 접촉하는지" 알려주는 테이블이 완성된다. 하지만 4,500개 개별 포즈를 하나하나 보는 건 의미가 없다. 과학적으로 의미 있는 단위는 **포즈가 아니라 포켓**이다.

blind docking을 하면 포즈가 표면 전체에 흩어지는데, 어떤 영역에는 포즈가 밀집하고 어떤 영역에는 띄엄띄엄 있다. 밀집한 영역이 바로 "리간드가 반복적으로 선호하는 포켓"이다. 이 단계는 그 밀집 패턴을 자동으로 식별해서, **포즈 레벨 데이터를 포켓 레벨로 올리는** 단계이다.

**왜 기존 클러스터링 알고리즘을 그대로 안 썼냐면:**

k-means는 클러스터 수(k)를 미리 정해야 하는데, receptor마다 포켓 수가 다르고 사전에 알 수 없다. DBSCAN은 밀도 기반이라 k가 필요 없지만, 도킹 포즈의 밀도 분포가 균일하지 않아서 파라미터 튜닝이 까다롭다.

그래서 **affinity 순으로 시드를 생성하는 탐욕적(greedy) 접근**을 선택했다. 에너지가 좋은 포즈부터 시드로 삼고, 그 시드에서 cutoff 이내의 포즈를 흡수하는 방식이다. 이렇게 하면 "에너지적으로 가장 유리한 영역"부터 자연스럽게 포켓으로 잡히고, 클러스터 수를 사전에 지정할 필요가 없다.

**3단계 구조로 만든 이유:**

**① Seed 생성** — 포즈를 affinity 순으로 스캔하면서, 기존 시드에서 7Å 이내면 그 시드에 할당하고, 아니면 새 시드를 만든다. 이게 초기 포켓 후보를 빠르게 잡는 단계이다.

**② 반복 할당** — 시드 생성 후, 모든 포즈를 가장 가까운 시드에 재할당하고 시드 centroid를 업데이트한다. 이걸 수렴까지(또는 max 15회) 반복한다. 초기 시드 위치가 첫 번째 포즈의 centroid에 의존하는데, 해당 포켓에 속한 다른 포즈들까지 고려하면 포켓 중심이 더 정확한 위치로 이동하기 때문이다.

**③ 후처리 병합** — ①②만으로는 **같은 포켓이 여러 개로 쪼개지는 문제**가 생길 수 있다. 큰 포켓의 양쪽 끝에 시드가 각각 생기면, 하나의 포켓이 두 개로 분절된다. 이걸 해소하기 위해 잔기 기반 병합을 추가했다. 두 포켓의 접촉 잔기 Jaccard가 ≥ 0.25이거나 Overlap coefficient가 ≥ 0.4이면, 또는 centroid 거리가 8Å 이내이면 병합한다. Union-Find 자료구조가 전이적 폐쇄(A-B 병합 + B-C 병합 → A-C 자동 병합)를 처리한다.

**centroid 거리와 잔기 기반 메트릭을 둘 다 쓰는 이유:**

각각의 한계를 보완하기 위해서이다. centroid 거리만 쓰면 잔기는 같은데 centroid가 미묘하게 다른 경우를 놓치고, 잔기 Jaccard만 쓰면 접촉 잔기 추출이 불완전한 경우를 처리하지 못한다. centroid fallback 8Å은 "잔기 데이터가 부족해도, 물리적으로 이 정도 가까우면 같은 포켓으로 봐야 한다"는 안전망이다.

**클러스터링 파라미터를 전부 기록하는 이유:**

`vina_clustering_parameters.json`에 모든 파라미터를 스냅샷으로, `vina_clustering_merge_log.csv`에 병합 사유를 기록한다. **"결정의 근거를 보존하라"** 원칙의 구현이다. "cutoff를 6Å로 바꾸면 포켓이 몇 개 더 나올까?" 같은 sensitivity 분석이 가능해진다.

---

## Step 5: Pocket Summary — 설계 의도

**해결하려는 문제:**

Step 4에서 4,500개 포즈가 포켓 단위로 그룹화되었다. 하지만 "포켓 3에 포즈 47개가 속해있다"는 것만으로는 아직 해석이 안 된다. 우리가 알고 싶은 건 **"이 포켓은 얼마나 강하고, 얼마나 수렴되어 있고, 화학적으로 얼마나 일반적인가"**이다.

이 단계의 설계 의도는 포즈 집합을 **포켓의 특성을 기술하는 통계**로 요약하는 것이다.

**어떤 통계를 왜 선택했냐면:**

**centroid_spread**는 포즈들이 포켓 중심으로부터 얼마나 퍼져있는지를 나타낸다. spread < 3Å이면 잘 정의된 포켓, > 5Å이면 경계가 모호하거나 여러 sub-site를 포함할 수 있다. 에너지가 좋은 포켓이라도 spread가 크면 "정말 하나의 잘 정의된 결합 부위인가"에 의문을 제기할 수 있다.

**best_affinity와 mean_affinity를 둘 다 기록**하는 이유가 있다. best_affinity만 보면 outlier 한 개에 속을 수 있다. mean_affinity가 함께 있으면 "이 포켓 전체가 에너지적으로 유리한가"를 판단할 수 있다.

**n_ligand**는 이 파이프라인에서 특히 중요한 지표이다. Step 1에서 설계한 cross-chemical consensus 논리가 여기서 숫자로 나타난다. 3이면 3종 리간드 모두 독립적으로 같은 포켓을 찾아간 것이고, Step 9 Verdict의 vina_consensus_pts 점수로 직접 반영된다.

**Drug-Pocket Map을 별도로 만드는 이유:**

`vina_drug_pocket_map.csv`는 포켓 중심이 아니라 **리간드 중심**의 관점이다. 각 (receptor, ligand) 쌍에서 dominant_pocket_fraction을 계산하는데, > 0.8이면 monodal(단일 결합 모드), 0.5-0.8이면 bimodal(두 포켓 경쟁), < 0.5이면 promiscuous(다수 부위 분산)이다.

이 정보가 필요한 이유는, Step 6의 cross-receptor 비교에서 "리간드 A가 3GT8에서는 포켓 1을 선호하는데 EGFR_160-185에서는 포켓 2로 바뀌었다"는 식의 state-dependent binding mode 변화를 포착할 수 있기 때문이다.

**Residue occupancy를 별도로 만드는 이유:**

pocket_table의 union_contact_residues는 "한 번이라도 접촉한 잔기" 전부를 포함하니까 노이즈가 섞인다. `vina_pocket_residue_occupancy.csv`는 "이 잔기는 포즈의 90%에서 접촉한다 vs 저 잔기는 5%에서만"을 구분할 수 있게 해주고, 90% 잔기가 그 포켓의 진짜 핵심이다.

---

## Step 6: Cross-Receptor 포켓 비교 — 설계 의도

**해결하려는 문제:**

Step 5까지 끝나면 receptor state마다 독립적인 포켓 프로파일이 완성되어 있다. 문제는 **이 포켓들이 서로 같은 포켓인지 다른 포켓인지 아직 모른다는 것**이다. pocket_id는 각 receptor 안에서 독립적으로 부여되기 때문에, 번호만 봐서는 대응 관계를 알 수 없다.

이 단계의 설계 의도는 **"receptor state가 달라도 같은 결합 부위인가"를 정량적으로 판정**하는 것이다. 이 프로젝트의 근본 질문 중 하나가 "EGFR의 동적 구조 변이 속에서도 일관되게 존재하는 포켓이 있는가"이기 때문이다.

**메트릭을 4개 쓰는 이유:**

**centroid_distance**는 가장 직관적이지만, 포켓 형태가 약간 달라져서 centroid가 이동한 경우 같은 포켓을 놓칠 수 있다. **residue_jaccard**는 비슷한 크기 포켓 비교에 적합하지만, 큰 포켓과 작은 포켓 비교에서 분모(합집합)가 커서 낮게 나올 수 있다. **residue_overlap_coefficient**는 부분집합 관계를 포착하고, **shared_ligands**는 기하학적 유사성 이상의 화학적 증거를 제공한다. 어떤 하나만으로는 false positive와 false negative를 동시에 막을 수 없기 때문에 네 가지를 조합한다.

**Same-Patch Candidate 판정 기준:**

`centroid_dist < 8Å AND (Jaccard ≥ 0.3 OR Overlap ≥ 0.5)` — centroid 거리는 필수 조건(잔기 ID 우연 겹침 방지), 잔기 기준은 Jaccard와 Overlap 중 하나만 충족하면 됨(둘의 맹점이 다르므로). 통계적으로, EGFR C-lobe 표면적 ~8000Å²에서 무작위로 8Å 이내에 겹칠 확률은 ~5-10%이고 잔기 중첩 조건까지 동시 충족해야 하므로, same-patch candidate는 p < 0.05 수준의 유의미한 일치이다.

**Bootstrap 안정성을 선택적으로 넣은 이유:**

200회 리샘플링은 계산 비용이 있는데, 대부분의 경우 cross-receptor 비교 자체가 이미 충분한 robustness 증거를 준다. **3개 독립 receptor state에서 같은 포켓이 나타난다는 것 자체가 일종의 자연적 bootstrap**이기 때문이다. Bootstrap은 단일 state 내 클러스터링 안정성을 추가로 확인하고 싶을 때의 보조 도구이다.

**이 단계의 위치:**

Step 6은 Vina branch의 최종 산출물을 만드는 단계이다. 동시에 Vina branch가 답할 수 있는 것의 한계가 명확해지는 지점이기도 하다. "어디에 소분자가 잘 붙는가, 그리고 그게 state에 걸쳐 안정한가"까지는 답할 수 있지만, **"그 포켓이 MYO1D 부착과 관련이 있는가"는 답할 수 없다**. 이 질문에 답하려면 PPI branch(Step 7-8)가 필요하다.

---

## Step 7: PyRosetta Global Blind PPI Docking — 설계 의도

**해결하려는 문제:**

Step 1-6은 "EGFR 표면 어디에 소분자가 잘 붙는가"를 답했다. 하지만 이 프로젝트의 최종 질문은 **"MYO1D 부착을 교란할 수 있는 포켓은 어디인가"**이다. 이 질문에 답하려면 "MYO1D가 EGFR 어디에 붙는가"를 먼저 알아야 한다.

Step 7-8은 Vina branch와 완전히 독립적인 PPI branch를 구성한다. 둘 다 독립적으로 실행한 뒤 Step 9 Verdict에서 사후적으로 수렴을 평가하는 것이 Workflow A의 핵심 설계이다. 이렇게 하면 두 branch가 서로를 오염시키지 않고, 수렴 자체가 증거가 된다.

**왜 PyRosetta를 쓰냐면:**

소분자 도킹과 단백질-단백질 도킹은 근본적으로 다른 문제이다. Vina는 소분자 torsion 자유도에 최적화된 경험적 scoring이고, 단백질 도킹은 6자유도 전체 탐색과 물리 기반(physics-based) scoring이 필요하다. PyRosetta의 ref2015 energy function(fa_atr, fa_rep, fa_sol, fa_elec, hbond, rama, fa_dun, ref의 가중합)이 단백질-단백질 interface의 상보성을 더 정확하게 평가할 수 있다.

**Global blind을 선택한 이유:**

Vina와 같은 논리이다. MYO1D가 EGFR 어디에 붙는지 확정되지 않은 상태에서, 특정 부위를 가정하고 시작하면 편향이 들어간다. RigidBodyPerturbMover(360°, 100Å)가 모든 가능한 배향을 탐색한다.

**왜 MYO1D 전체가 아니라 beta-meander(955-1006)를 쓰냐면:**

이건 시행착오에서 나온 결정이다. 처음에는 TH1 domain 전체를 도킹 입력으로 시도했는데, 큰 도메인이 EGFR 표면 전체에 비특이적으로 붙어서 노이즈가 너무 많았다. 실험적으로 MYO1D의 beta-meander 영역(β-sheet 8-12)이 EGFR과의 상호작용에 관여한다는 근거가 있어서, 이 영역을 도킹 입력으로 좁혔다.

잔기 범위도 조정이 있었다. 처음 pilot에서는 962-1006으로 잘랐는데, VAL962가 chain의 첫 번째 잔기가 되면서 인공적인 N-terminal charge와 과도한 backbone freedom 때문에 모든 도킹 모델에서 100% 접촉 잔기로 나타나는 artifact가 생겼다. 이걸 해결하기 위해 시작점을 ~955로 연장해서 VAL962가 더 이상 chain 끝이 아니도록 만들었다. 마찬가지로 receptor도 C-lobe fragment(45 잔기)에서 full kinase domain(~280 잔기)으로 확장했는데, N-lobe가 없으면 C-lobe 표면의 정전기 landscape와 steric environment가 달라져서 도킹 결과가 왜곡되기 때문이다.

**20K × 5 seeds 설계:**

단일 100K 대신 20K × 5 seeds를 선택한 이유가 세 가지이다. 첫째, **수렴 진단** — 5개 독립 시드에서 같은 interface patch가 반복되면 수렴된 것이고, 특정 시드에서만 나타나면 artifact이다. 둘째, **패치 재현성 평가** — "이 잔기가 5개 시드 중 몇 개에서 interface에 참여하는가"라는 cross-seed robustness를 정량화할 수 있다. 셋째, **HPC 작업 관리** — 실패한 시드만 재실행하면 된다.

**6-Step 내부 파이프라인:**

**Relax** — 입력 구조의 steric clash 해소. relaxed_cache에 캐싱하여 반복 사용 시 재계산 방지.

**Global Docking** — 실제 6-DOF 탐색. SlideIntoContact가 필수인 이유는, 이 없이 DockMCMProtocol을 호출하면 모든 dG가 0.0으로 나오는 역사적 버그가 있었기 때문이다.

**Scoring & Filtering** — v2.0 2-pass 설계. Pass 1에서 total_score percentile로 대부분 걸러내고, Pass 2에서 생존 모델에만 비용 큰 정밀 메트릭을 계산한다. 300K 전부에 정밀 scoring을 하면 비용이 감당 안 되기 때문이다. min_survivors=50은 과도한 필터링 방지 안전장치이다.

**Clustering** — L_RMSD greedy. Vina의 centroid 기반과 다른 이유는, 200+ 잔기 파트너의 전체 배향을 비교해야 하기 때문이다. cluster_top_n = auto(15-35)는 모델 수에 따라 적절한 클러스터 수가 다르기 때문이다.

**Refine & Select** — 대표 구조당 10회 미세 섭동(0.1Å/1.0°)으로 local minimum을 정밀하게 찾는다. Round-robin 선택은 최종 20 모델이 여러 클러스터에서 고르게 나오도록 보장한다.

**Viz & Report** — PyMOL 스크립트와 PPI Validation Report로 시각적 검증 경로를 내장한다. 비물리적 배향(membrane 쪽 도킹 등)을 바로 잡을 수 있다.

---

## Step 8: PPI Interface 추출 & Cross-State Patch — 설계 의도

**해결하려는 문제:**

Step 7에서 300K개 도킹 모델이 생성되었다. 하지만 우리가 궁극적으로 원하는 건 "어떤 모델이 최고인가"가 아니라 **"EGFR의 어떤 잔기가 반복적으로 MYO1D interface에 참여하는가"**이다.

이 관점 전환이 핵심이다. 모델 수준에서 잔기 수준으로, "best model"에서 **"robust patch"**로. 단백질-단백질 도킹의 개별 모델은 신뢰도가 낮지만, 300K 모델에서 같은 잔기가 반복적으로 interface에 나타난다면, 그건 에너지 landscape의 실제 신호일 가능성이 높다.

**데이터 변환 흐름을 이 순서로 설계한 이유:**

**300K raw → scored_all_models.csv** — 탈락 모델도 기록하여 post-hoc 재분석을 가능하게 한다.

**scored_all_models → pyrosetta_decoy_scores.csv** — standardize_scores.py가 5개 시드 × 3개 state의 결과를 하나의 표준 스키마로 통합한다. 이 시점부터 cross-seed 비교가 가능해진다.

**decoy_scores → ppi_pyrosetta_residues.csv** — 핵심적 관점 전환. Binding_Residues_A/B를 파싱하여 잔기 단위 테이블로 풀어내고, N-lobe / C-lobe 분류를 한다. MYO1D 상호작용이 C-lobe에서 예상되므로, N-lobe 접촉은 steric context이고 C-lobe 접촉이 실제 후보이다.

**residues → ppi_hotspot_residues.csv** — 빈도를 세서 반복적으로 interface에 참여하는 hotspot 잔기를 식별한다. 여러 시드, 여러 클러스터에서 독립적으로 나타나는 잔기가 진짜 hotspot이다. 20K × 5 seeds 설계가 이 cross-seed robustness를 가능하게 한다.

**hotspot → ppi_patch_state_robustness.csv** — 3개 state에 걸친 patch robustness 평가. Vina Step 6과 평행한 논리이되, "소분자 포켓"이 아니라 "PPI patch"의 state 안정성을 본다.

**LightDock 독립 검증의 설계 의도:**

단일 방법론의 결과만 믿는 건 위험하다. LightDock는 PyRosetta와 완전히 다른 탐색 알고리즘(Glowworm Swarm Optimization)과 scoring function을 쓴다. 두 방법이 같은 잔기를 지목하면 방법론 특이적 artifact가 아닌 실제 신호이다. cross_method_convergence.csv의 method_agreement 컬럼이 `both`, `pyrosetta_only`, `lightdock_only`로 구분한다.

LightDock을 secondary로 위치시킨 이유는, PyRosetta가 더 많은 모델(300K), 더 정밀한 scoring(ref2015), 더 체계적인 filtering(v2.0 2-pass)을 거치기 때문이다. LightDock은 독립 검증 역할이지 대체가 아니다.

**Orientation filter가 왜 필요하냐면:**

이건 MYO1D beta-meander의 구조적 특성에서 오는 고유한 문제이다. beta-meander는 5개 연속 β-strand로 이루어진 평평한 β-sheet 리본인데, 이 구조는 두 면이 있다. 실험적으로 검증된 결합면(active face, sheet 8/9의 side chain이 receptor를 향하는 배향)과 그 반대면(back face)이다.

문제는 이 sheet가 매우 얇아서(두 면 사이 거리 ~5-7Å) 8Å contact cutoff를 쓰면 양면의 잔기가 동시에 "접촉"으로 잡힌다는 것이다. 즉, beta-meander가 뒤집혀서 back face가 receptor를 향하고 있어도 contact counting만으로는 올바른 배향과 구분이 안 된다.

그래서 sheet 8/9의 active-face normal vector와 receptor 방향 vector의 dot product를 계산해서, 양수면 pass(active face가 receptor를 향함), 음수면 fail(뒤집힘), 0 근처면 ambiguous로 분류한다. fail 모델이 잔기 빈도 계산에 들어가면 노이즈가 되니까, pass 모델만 consensus에 사용한다.

**이 단계의 위치:**

Step 8은 PPI branch의 최종 산출물을 만든다. 이 시점에서 두 branch가 준비된다. Vina: "이 포켓은 여기에 있고, 이 잔기와 접촉하고, state에 걸쳐 보존된다." PPI: "이 receptor-side patch는 이 잔기들로 구성되고, state에 걸쳐 robust하다." 다음 Step 9에서 이 두 산출물이 만난다.

---

## Step 9: 3축 Verdict 통합 판정 — 설계 의도

**해결하려는 문제:**

Step 6의 Vina 증거와 Step 8의 PPI 증거가 독립적으로 준비되었다. 남은 질문은 **"이 둘이 겹치는가, 그리고 그 증거가 얼마나 강한가"**이다. 이 질문에 "겹친다/안 겹친다"로 답하면 안 된다. 증거의 강도를 정량화해야 한다.

**왜 3축 체계인가:**

하나의 종합 점수가 아니라 3개 독립 축으로 분리한 이유는, 각 축이 서로 다른 종류의 증거를 대표하기 때문이다.

**축 1: Vina 증거 (50점)** — "소분자가 이 포켓에 잘 붙는가?" 5개 sub-score로 구성: vina_affinity_pts(에너지), vina_convergence_pts(포즈 수렴), vina_consensus_pts(다중 리간드 합의), vina_stability_pts(리샘플링 안정성), vina_diversity_pts(리간드 다양성). 5개를 따로 둔 이유는, **어떤 증거가 부족한지를 투명하게 보여주기 위해서**이다. affinity는 좋은데 convergence가 낮으면 "샘플링 부족 가능성"을, convergence는 높은데 consensus가 낮으면 "화학 특이적 포켓 가능성"을 시사한다.

**축 2: PPI 공간 관계 (20점)** — "이 포켓이 MYO1D interface 근처인가?" 이 축이 프로젝트의 핵심이다. ppi_spatial_pts(centroid 거리), ppi_overlap_pts(잔기 Jaccard), ppi_reproducibility_pts(PPI 잔기 재현성)로 구성. PPI 증거 자체가 약하면, spatial이 좋아도 신뢰하기 어렵기 때문에 reproducibility를 별도로 둔다.

**축 3: Cross-Receptor (30점)** — "이 증거가 구조 변이에도 유지되는가?" 이 축에 30점이라는 큰 비중을 준 이유는, 이 프로젝트가 state-comparison pipeline이기 때문이다. 단일 구조에서의 좋은 결과보다, 여러 구조에서의 일관된 결과가 더 가치 있다.

**적응적 점수 배분:**

PPI 데이터가 있을 때 50+20+30=100, 없을 때 60+0+40=100으로 자동 조정된다. PPI 축을 0으로 처리하면 만점이 달라져서 비교가 불공정해지므로, 가중치를 재배분하여 만점을 항상 100으로 유지한다.

**판정 기준:**

STRONG ≥ 55점. 이 기준의 의미를 거꾸로 생각하면 이해가 된다. Vina 축 만점(50)만으로는 55에 못 미치고, Cross-Receptor 만점(30)만으로도 못 미친다. **STRONG이 되려면 최소 2개 축에서 의미 있는 점수를 받아야 한다.** 단일 증거 유형만으로는 STRONG을 줄 수 없고, 다중 증거가 수렴해야 한다는 원칙이다.

**"증거 분류이지 타당성 판정이 아니다":**

STRONG은 "이 포켓이 실제로 MYO1D를 교란한다"는 뜻이 아니라, "계산적 증거가 여러 축에서 수렴한다"는 뜻이다. STRONG도 반드시 PyMOL에서 시각적으로 검증해야 한다. 반대로 WEAK도 cryptic pocket일 수 있다. Verdict는 "어디를 먼저 봐야 하는가"의 우선순위이지, "맞다/틀리다"의 판정이 아니다.

---

## Step 10: Report 생성 & 출력 검증 — 설계 의도

**해결하려는 문제:**

Step 9까지 끝나면 과학적 분석은 완료이다. 하지만 두 가지 문제가 남는다. CSV 파일을 열어야만 결과를 볼 수 있다는 것, 그리고 출력 파일이 정말 온전한지 아무도 확인하지 않았다는 것이다.

**Report (report.py):**

project_report.txt의 5개 섹션이 **결과를 읽는 사람의 질문 순서**를 따르도록 설계했다. (1) 입력 요약 → (2) Vina 포켓 결과 → (3) PPI interface 잔기 → (4) 교차 방법 비교 → (5) 최종 판정 및 권장사항. "무엇이 입력이었고 → 소분자 결과 → PPI 결과 → 둘이 겹치는지 → 최종 판단"이 자연스럽게 따라간다.

combined_residue_evidence.csv는 Vina 접촉 잔기와 PPI interface 잔기를 **잔기 단위로 통합**한다. valid_sites.csv가 포켓 단위인 반면, 때로는 "잔기 987이 왜 중요한가?"라는 잔기 수준 질문에 직접 답해야 할 때가 있다.

vina_consensus_sites.csv는 verdict.py에서 생성된다(report.py가 아님). 이건 의도적 역할 분리이다. verdict.py는 판정을, report.py는 그 판정의 포맷팅을 담당한다.

**Validate (validate.py):**

4개 검증 그룹(8개 함수)이 **문제를 빠르게 좁혀갈 수 있는** 구조로 설계되었다.

**(8.1) 파일 존재 + ID 일관성 + 추적성 + coverage** — 파이프라인이 중간에 죽었는지의 가장 기본적인 체크.

**(8.2) CSV 스키마 회귀 검사** — 필수 컬럼 누락, 데이터 타입 불일치. 코드 수정 시 출력 스키마를 의도치 않게 바꿀 수 있으므로 특히 중요하다.

**(8.3) 잔기 번호 일관성 + 알려진 변이 확인** — 3개 receptor state의 번호 체계가 맞는지. 3GT8은 PDB 번호, MD 클러스터는 다른 오프셋을 가질 수 있고, 이게 어긋나면 cross-receptor 비교와 PPI overlap 판정이 전부 틀어진다.

**(8.4) 핸드오프 준비 확인** — Workflow B나 외부 전달을 위한 파일/메타데이터 완비 여부.

C1-C10 체크는 pipeline_manager.py의 PPI validation report(Step 7 내부 품질 체크)에 해당하며, validate.py와는 별개이다.

종료 코드 0(통과)/1(경고)/2(실패)는 PBS 스크립트에서 후속 작업 중단을 자동화하기 위한 설계이다.

**Step View의 설계 의도:**

step_view.py가 step1_vina_raw/ ~ step7_validate/ 파생 폴더를 만드는 것은, **canonical output과 interpretation view를 분리**하기 위해서이다. 원본은 `output/egfr_myo1d_vina/`에 그대로 있고(source of truth), step view는 탐색 편의성만 추가하는 파생 뷰이다.

**정리:**

Step 10은 "계산 결과를 사람의 언어로 번역하라"(report)와 "조용히 깨지는 것을 방지하라"(validate)를 담당한다. 이 프로젝트가 반복적으로 강조하는 원칙 — **"결론이 아니라 증거를 남겨라"** — 의 마지막 구현이다.
