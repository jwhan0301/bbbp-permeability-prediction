# SMILES 기반 혈액–뇌 장벽 투과성 예측

SMILES로 표현된 분자가 혈액–뇌 장벽(Blood–Brain Barrier, BBB)을 투과하는 경향이 있는지 분류하는 교육·연구용 기준모델 프로젝트입니다.

## 30초 요약

- MoleculeNet BBBP 데이터 2,050행을 점검하고, RDKit이 해석할 수 있으며 중복·라벨 충돌 처리를 통과한 1,965개 분자를 모델링에 사용했습니다.
- 현재 앱은 RDKit Morgan fingerprint(`radius=2`, `fpSize=2048`)와 랜덤 포레스트(Random Forest)를 사용합니다.
- 현재 앱 모델의 고정 시험 데이터 ROC-AUC는 **0.895986**입니다.
- BBB+ 민감도는 **0.980000**이지만 BBB- 특이도는 **0.569892**로, BBB- 식별력이 상대적으로 낮습니다.
- Streamlit 데모는 BBB 예측과 함께 가장 유사한 학습 분자와 타니모토 유사도(Tanimoto similarity)를 보여줍니다.

## 앱 화면

![Streamlit BBBP 예측 데모](docs/images/bbbp_streamlit_demo.png)

SMILES를 입력하면 BBB+/BBB- 예측과 BBB+ 모델 점수를 표시합니다. 가장 유사한 학습 분자와 타니모토 유사도도 함께 보여주지만, 이 유사도는 예측 신뢰확률이나 정확도를 보장하지 않습니다.

## 프로젝트 흐름

```text
SMILES 입력
→ RDKit으로 분자 구조 해석
→ Morgan fingerprint 생성
→ Random Forest 예측
→ BBB 점수와 가장 유사한 학습 분자 표시
```

- **SMILES:** 원자와 결합을 문자로 나타낸 분자 구조 표현입니다.
- **Morgan fingerprint:** 분자 주변 구조의 특징을 0과 1로 바꾼 2,048비트 분자 지문입니다.
- **타니모토 유사도:** 두 분자 지문에서 켜진 비트가 얼마나 겹치는지 0~1로 나타낸 구조 유사도입니다.

## 핵심 결과

현재 Streamlit 앱은 3일 차에 만든 기존 랜덤 포레스트를 사용합니다. 최종 모델용 분자 1,965개를 클래스 비율이 유지되도록 80:20으로 나누고(`random_state=42`), 학습 데이터 1,572개로 학습한 뒤 따로 보관한 시험 데이터 393개에서 한 번 평가했습니다. 분류 기준은 0.5입니다.

| 지표 | 고정 시험 결과 |
|---|---:|
| 정확도(Accuracy) | 0.882952 |
| F1 점수 | 0.927445 |
| ROC-AUC | 0.895986 |
| 균형 정확도(Balanced Accuracy) | 0.774946 |
| MCC | 0.654306 |
| 민감도(BBB+ 재현율) | 0.980000 |
| 특이도(BBB- 재현율) | 0.569892 |
| TN / FP / FN / TP | 53 / 40 / 6 / 294 |

BBB+는 잘 찾았지만 BBB-를 BBB+로 잘못 분류한 경우가 상대적으로 많았습니다. BBB+가 전체 모델용 데이터의 76.34%이므로 정확도만 보면 성능을 실제보다 좋게 판단할 수 있습니다. 또한 이 수치는 한 번의 무작위 분할에서 얻은 결과이며, 다른 데이터 분할이나 외부 데이터에서 달라질 수 있습니다.

## 6일 차 개선 후보

클래스 가중치를 적용한 Balanced RF는 **기존 학습 데이터 1,572개만 사용한 5겹 교차검증 후보**입니다. 아래 값은 같은 교차검증 조건에서 비교한 평균이며, 위 고정 시험 결과와 직접 비교하는 수치가 아닙니다.

| 학습 데이터 전용 교차검증 지표 | 기존 RF | Balanced RF |
|---|---:|---:|
| 균형 정확도 | 0.794975 | 0.837782 |
| MCC | 0.674409 | 0.707550 |
| 민감도 | 0.974167 | 0.952500 |
| 특이도 | 0.615784 | 0.723063 |

Balanced RF는 균형 정확도·MCC·특이도가 높아졌지만 민감도는 낮아졌습니다. 기존 시험 데이터는 후보 선택이나 재평가에 사용하지 않았습니다. 아직 외부 데이터 검증을 하지 않았으므로 앱 모델도 Balanced RF로 교체하지 않았습니다.

## 빠른 실행 방법

Windows PowerShell에서 저장소 루트로 이동한 뒤 다음 명령을 실행합니다. 패키지 설치에는 인터넷 연결이 필요합니다.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts\smoke_test.py
python -m streamlit run app.py
```

Streamlit이 출력한 주소(보통 `http://localhost:8501`)를 브라우저에서 엽니다. macOS나 Linux에서는 가상환경 활성화 명령만 `source .venv/bin/activate`로 바꾸면 됩니다.

저장된 앱 모델을 실행하는 데 BBBP 원본 CSV는 필요하지 않습니다. 모델과 분석을 처음부터 다시 만드는 방법은 [전체 실험 기록](docs/EXPERIMENT_LOG.md)과 [재현성 검증 기록](docs/REPRODUCIBILITY.md)을 참고하세요.

## 주요 폴더 구조

| 경로 | 역할 |
|---|---|
| `app.py` | SMILES 한 개를 입력받는 Streamlit 앱 |
| `src/` | 분자 특징 생성, 예측, 유사도 계산 코드 |
| `scripts/` | 모델·유사도 자료 재생성과 작동 검사 스크립트 |
| `notebooks/` | 데이터 확인부터 모델 검증까지의 단계별 노트북 |
| `models/` | 현재 앱 모델과 학습 분자 유사도 참조자료 |
| `results/` | 실제 실행으로 생성한 성능표와 그림 |
| `docs/` | 전체 실험 일지, 재현성 기록과 앱 이미지 |
| `MODEL_CARD.md` | 모델 데이터·설정·평가·사용 한계 |
| `LICENSE` | 이 저장소에서 직접 작성한 코드의 MIT License |
| `requirements.txt` | 실행을 확인한 Python 패키지 버전 |

## 한계

- 최종 모델용 분자는 1,965개로 작고 BBB+가 76.34%인 불균형 데이터입니다.
- 현재 앱 모델의 BBB- 특이도는 0.569892로 BBB- 식별력이 낮습니다.
- 최종 평가는 `random_state=42`인 한 번의 무작위 분할 결과입니다.
- 화학적으로 유사한 골격(scaffold)이 학습과 시험 데이터에 나뉘어 결과가 낙관적일 수 있습니다.
- 골격 기준 분할과 외부 데이터 검증이 부족합니다.
- 수송체, 대사, 농도와 실험 조건 같은 생물학적 정보를 입력에 사용하지 않습니다.
- 이 모델과 앱은 임상 판단이나 실제 신약 개발 의사결정에 사용할 수 없습니다.

## 상세 문서

- [모델 카드](MODEL_CARD.md): 데이터, 모델 설정, 평가 결과와 사용 한계
- [1~6일 차 전체 실험 기록](docs/EXPERIMENT_LOG.md): 전처리, 특징 생성, 학습, 검증, 오류 분석과 생성 파일
- [독립 환경 재현성 검증](docs/REPRODUCIBILITY.md): 새 환경 설치와 실행 확인 결과
- [실행 가능한 단계별 노트북](notebooks/): 분석 코드를 처음부터 실행할 수 있는 노트북
- [결과 파일](results/): 실제 실행으로 생성한 CSV와 PNG

## 데이터 출처와 라이선스

- 데이터: MoleculeNet BBBP, DeepChem의 [BBBP 불러오기 코드](https://github.com/deepchem/deepchem/blob/master/deepchem/molnet/load_function/bbbp_datasets.py)에 명시된 공개 CSV를 실행 시 불러옵니다.
- 원 논문: Martins et al. (2012), [A Bayesian Approach to in Silico Blood-Brain Barrier Penetration Modeling](https://doi.org/10.1021/ci300124c)
- BBBP 원본 CSV는 이 저장소에 포함하지 않습니다. 데이터의 권리와 이용 조건은 원 데이터 제공처의 조건을 따릅니다.
- 이 저장소에서 직접 작성한 코드에는 [MIT License](LICENSE)가 적용됩니다.
- 제3자 데이터와 라이브러리는 각각의 원래 라이선스 및 이용 조건을 따르며, 이 저장소의 MIT License는 BBBP 원본 데이터에 적용되지 않습니다.
