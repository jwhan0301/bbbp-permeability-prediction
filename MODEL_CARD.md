# BBBP Random Forest 모델 카드

## 모델 개요

이 프로젝트는 분자의 SMILES를 입력받아 혈액-뇌 장벽(Blood-Brain Barrier, BBB)을 투과하는 경향이 있는지 `BBB+` 또는 `BBB-`로 분류하는 교육용 예시입니다.

- 입력: 유효한 SMILES 문자열 한 개
- 출력: `BBB+` 또는 `BBB-` 분류, 보정되지 않은 BBB+ 모델 점수, 가장 유사한 Train 분자와 Tanimoto similarity
- 사용 목적: 분자 표현, 분류 모델, 평가 지표와 간단한 예측 앱을 학습하기 위한 baseline
- 현재 Streamlit 앱 모델: Day 3 Morgan fingerprint Random Forest
- 저장 형식 버전: `bundle_version=1`
- 저장 모델: `models/bbbp_random_forest.joblib`

Day 6에서 Train-only 개선 후보를 비교했지만, 기존 Test를 다시 사용해 선택하지 않기 위해 앱 모델은 교체하지 않았습니다.

## 데이터

데이터는 MoleculeNet에 포함된 BBBP 공개 데이터이며, 프로젝트에서는 DeepChem의 공개 URL에서 실행할 때 불러옵니다.

- 데이터 URL: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- 원본 행 수: 2,050
- RDKit SMILES 파싱 성공: 2,039
- RDKit SMILES 파싱 실패: 11
- 동일 라벨 중복으로 모델 대상에서 제외: 54행
- 라벨 충돌로 모델 대상에서 제외: 20행
- 최종 모델용 분자: 1,965
- 최종 BBB+: 1,500개(76.34%)
- 최종 BBB-: 465개(23.66%)

원본 행은 저장소에 복사하거나 이유 없이 삭제하지 않았습니다. 파싱 실패 분자는 구조 특징을 계산할 수 없어 모델 학습 대상에서 제외했습니다. Canonical SMILES가 같고 라벨도 같으면 원본 순서의 첫 행만 유지했습니다. 같은 canonical SMILES에 서로 다른 라벨이 있으면 임의로 라벨을 고치지 않고 해당 충돌 그룹 전체를 모델 대상에서 제외했습니다.

최종 1,965개를 `stratified 80:20 train-test split`, `random_state=42`로 나눴습니다.

- Train: 1,572개(BBB- 372, BBB+ 1,200)
- Test: 393개(BBB- 93, BBB+ 300)

## 모델과 분자 표현

현재 앱 모델의 입력은 descriptor가 아니라 Morgan fingerprint입니다.

- Morgan radius: 2
- Fingerprint 크기: 2,048 bit
- 생성 방식: RDKit `FingerprintGenerator.GetMorganGenerator`
- Fingerprint 값: 0과 1
- 모델: `RandomForestClassifier`
- `n_estimators=300`
- `criterion="gini"`
- `max_depth=None`
- `max_features="sqrt"`
- `min_samples_split=2`
- `min_samples_leaf=1`
- `bootstrap=True`
- `class_weight=None`
- `random_state=42`
- `n_jobs=-1`
- 분류 threshold: 0.5
- Descriptor 사용 여부: 현재 앱 모델에서는 사용하지 않음

앱의 전처리와 예측 흐름은 다음과 같습니다.

```text
SMILES 입력
→ 빈 값·길이·RDKit 파싱 검사
→ 입체화학 정보를 가능한 범위에서 유지한 canonical SMILES 생성
→ Morgan fingerprint (1, 2048) 생성
→ 저장된 Random Forest에서 BBB+ 점수 계산
→ 0.5 기준 BBB+/BBB- 분류
→ Train fingerprint와 Tanimoto similarity 계산
→ 가장 유사한 Train 분자와 경험적 유사도 경고 표시
```

BBB+ 점수가 정확히 0.5이면 scikit-learn의 클래스 선택과 맞추기 위해 BBB-로 분류합니다.

## 평가 결과

### 현재 앱에서 사용하는 기존 Random Forest

아래 값은 Day 3에 봉인한 고정 Test 393개에서 한 번 평가한 결과입니다.

| 지표 | 값 |
|---|---:|
| Accuracy | 0.882952 |
| F1-score | 0.927445 |
| ROC-AUC | 0.895986 |
| Balanced Accuracy | 0.774946 |
| MCC | 0.654306 |
| Sensitivity(BBB+ recall) | 0.980000 |
| Specificity(BBB- recall) | 0.569892 |

Confusion matrix의 클래스 순서는 실제 `[BBB-, BBB+]` × 예측 `[BBB-, BBB+]`입니다.

|  | 예측 BBB- | 예측 BBB+ |
|---|---:|---:|
| 실제 BBB- | TN 53 | FP 40 |
| 실제 BBB+ | FN 6 | TP 294 |

이 모델은 BBB+를 찾는 Sensitivity는 높았지만 BBB-를 찾는 Specificity는 상대적으로 낮았습니다.

### Day 6 Balanced RF 개선 후보

Day 6의 후보 B는 기존 Random Forest와 같은 Morgan fingerprint와 Random Forest 설정을 사용하고 `class_weight="balanced"`만 추가했습니다. Day 3 Train 1,572개만 사용한 동일 5-fold 교차검증에서 사전 기준에 따라 개선 후보로 선택됐습니다.

| 지표 | 5-fold 평균 ± 표준편차 |
|---|---:|
| Balanced Accuracy | 0.837782 ± 0.050489 |
| MCC | 0.707550 ± 0.084956 |
| Sensitivity | 0.952500 ± 0.013693 |
| Specificity | 0.723063 ± 0.095080 |

기존 후보 A의 Train-only CV와 비교하면 Specificity는 0.615784에서 0.723063으로 높아졌고 Sensitivity는 0.974167에서 0.952500으로 낮아졌습니다. 즉 BBB- 식별 개선과 BBB+ 식별 감소가 함께 나타난 trade-off입니다.

기존 Test 393개는 후보 선택이나 후보 재평가에 사용하지 않았습니다. 외부 데이터 또는 새로운 미사용 평가셋에서 검증하기 전이므로 Balanced RF를 현재 앱 모델로 교체하지 않았습니다.

## 구조 유사도 기능

Tanimoto similarity는 두 Morgan fingerprint에서 켜진 bit의 겹침 정도를 0~1로 나타낸 구조 유사도입니다. 앱은 입력 분자가 학습 데이터의 어떤 분자와 비슷한지 사용자가 확인할 수 있도록 가장 유사한 Train 분자를 보여줍니다.

경고선을 임의의 0.5나 0.7로 정하지 않았습니다. Train 1,572개 각각에 대해 자기 자신을 제외한 가장 가까운 다른 Train 분자의 similarity를 계산했습니다. 이 분포의 10번째 백분위수인 **0.3000**을 경험적 경고선으로 사용합니다.

입력 분자의 최대 Train similarity가 0.3000보다 낮으면 학습 데이터의 일반적인 구조 유사도 범위보다 낮다는 주의 문구를 표시합니다. 이 값은 현재 Train 데이터에 의존하며 예측 신뢰확률, 정확도 또는 통계적 신뢰수준을 보장하지 않습니다.

## 의도된 사용

- 학습 및 연구용 BBB permeability baseline
- SMILES, Morgan fingerprint와 분류 모델 학습 예시
- Accuracy 외에 Sensitivity, Specificity, Balanced Accuracy와 MCC를 비교하는 예시
- 후보 물질의 예측과 학습 데이터 구조 유사도를 대략적으로 살펴보는 데모

## 사용하면 안 되는 경우

- 임상적 판단
- 환자 치료 결정
- 실제 신약 후보의 안전성 또는 효능 보증
- 사람 또는 동물 실험을 대신하는 용도
- 실제 BBB 투과 확률, 뇌 농도 또는 logBB를 정량적으로 추정하는 용도

## 알려진 한계

- 최종 모델용 분자가 1,965개인 비교적 작은 데이터셋입니다.
- BBB+가 76.34%로 클래스가 불균형합니다.
- 현재 앱 모델의 고정 Test Specificity는 0.569892로 BBB- 식별력이 상대적으로 낮습니다.
- 한 번의 stratified random split만 최종 Test 평가에 사용했습니다.
- Scaffold split, 반복 split과 외부 데이터 검증이 부족합니다.
- 실험 조건, 농도, 대사, efflux/influx transporter와 같은 생물학적 정보를 입력으로 사용하지 않습니다.
- 학습 데이터와 구조적으로 다른 분자에서는 모델의 동작이 더 불확실할 수 있습니다.
- Random Forest 점수는 확률 보정을 하지 않았습니다.

## 재현 방법

저장소 루트에서 새 가상환경을 만들고 다음 순서로 실행합니다.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts\smoke_test.py
python -m streamlit run app.py
```

`scripts/smoke_test.py`는 저장된 모델과 유사도 artifact를 불러와 정상 입력, 잘못된 입력, 기존 예측 점수, fingerprint 설정과 구조 유사도 계산을 검사합니다.
