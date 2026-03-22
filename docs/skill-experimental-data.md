# Skill: 실험 데이터 통합

## 적용 범위
- Group: 1
- Features: F-1

## 도메인 지식

### ATP Binding Site 정의 (확정)
- 핵심 5잔기: Lys745, Thr790, Gln791, Met793, Asp855 — ATP와 직접 접촉
- 확장 37잔기: `region_definitions.py`의 `ATP_SITE_RESIDUES` — P-loop(719-723), hinge(788-796), catalytic loop(831-837), A-loop(854-858) 포함
- 포켓의 ATP site 판정: 접촉 잔기 중 ATP_SITE_RESIDUES 비율 > 50% → `is_atp_site = True`
- 실험적 근거: 약물 처리 시 ATP 결합은 유지되면서 kinase 활성만 소실 → ATP site 포켓은 이 프로젝트에서 false positive

### Ko et al. Sheet 잔기 (확정)
| Sheet | 잔기 범위 | 실험 결과 | 파이프라인 용도 |
|-------|-----------|-----------|----------------|
| 8 | 961-964 (4개) | Ala sub → 기능 소실 | Active face, orientation filter |
| 9 | 968-972 (5개) | Ala sub → 기능 소실 | Active face, orientation filter |
| 10 | 975-984 (10개) | Ala sub → WT-level | Neutral — hotspot 출현 시 WARNING |
| 11 | 985-991 (7개) | Ala sub → WT-level | Neutral — hotspot 출현 시 WARNING |
| 12 | 998-1004 (7개) | Ala sub → 기능 소실 | Structural support 추정, 모니터링 |

- Active face = sheet 8 | sheet 9 = 9잔기 (Ko et al.에서 direct contact 확인)
- Ko et al. 체크: hotspot에 active face 3개 미만 → FAIL (PPI 예측이 실험과 괴리)

### 리간드 다양성 기준
- 쌍별 Tanimoto (Morgan fp, radius=2, nBits=2048)
- < 0.4: "화학적으로 다양" → cross-chemical consensus 유효
- 0.4-0.7: "중간" → 유효하나 해석 시 주의
- \> 0.7: "화학적으로 유사" → consensus 가치 재평가 필요
- 리간드 3종: 173940, 97806, VAX-C12_0 (SDF 포맷, input/ 디렉토리)

### Sheet 접촉 정보 해석
- `contacts_sheet_8_9` 높음: 약물이 active face를 직접 교란 → orthosteric 가능성
- `contacts_sheet_12` 높음: structural support 영역 교란 → 간접 불안정화
- `contacts_sheet_10_11`만 높음: 실험적으로 기능 영향 없는 영역 → 교란 효과 낮음

## 주의 사항
- `is_atp_site` 경계(40-60% 겹침) 포켓은 `is_atp_site_borderline` 추가 태깅 필요
- Ko et al. 체크의 FAIL은 Workflow B 전체를 중단시키므로, 테스트에서 정상 CSV도 반드시 PASS 확인
- sheet 잔기 번호는 MYO1D(chain B) 기준이며, EGFR(chain A) 잔기 번호와 혼동하지 않을 것
- `exclusion_reason` 컬럼 추가 시 기존 valid_sites.csv 파서가 깨지지 않는지 반드시 확인
