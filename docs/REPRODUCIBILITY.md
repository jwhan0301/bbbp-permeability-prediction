# 독립 재현성 검증 기록

검증일: 2026-08-25

이 문서는 GitHub 공개 전에 첫 커밋을 실제로 `git clone`하고, 기존 프로젝트의 `.venv`를 사용하지 않은 새 Python 환경에서 설치와 실행을 확인한 기록입니다.

## 검증 환경

- 운영체제: Windows 11 (`10.0.26200`, AMD64)
- Python: 3.12.13
- 패키지 목록: 저장소의 `requirements.txt`
- 검증 커밋: `b63d086d9d0ab3d799e827ad46d96e1541a2ac01`
- 검증 방법: 로컬 저장소의 `main` 브랜치를 `git clone --no-local`로 새 폴더에 복제한 뒤, 복제본 안에서 새 가상환경을 생성했습니다.
- 복제본 상태: 추적 파일 52개, 미커밋 변경 0개
- 제외한 항목: 기존 `.venv`, `work/`, Python 캐시, notebook checkpoint, 원본 BBBP CSV, 환경변수 파일과 사용자별 설정

## 실행한 검증

| 검증 항목 | 결과 | 확인 내용 |
|---|---:|---|
| 새 가상환경 생성 | PASS | 기존 프로젝트 `.venv`와 분리된 Python 3.12.13 환경 생성 |
| `requirements.txt` 설치 | PASS | 고정된 패키지 버전을 새 환경에 설치 |
| 패키지 의존성 검사 | PASS | `pip check` 결과 `No broken requirements found.` |
| 자동 smoke test | PASS | 입력 오류 처리, 기존 예측 점수, 유사도 범위, train 자기 구조 검색을 확인 |
| 저장 모델 불러오기 | PASS | `models/bbbp_random_forest.joblib`을 새 환경에서 불러옴 |
| 정상 SMILES 예측 | PASS | 카페인 예시를 BBB+로 예측했고 점수는 `0.9700` |
| 빈 입력 처리 | PASS | `SMILES가 비어 있습니다.` 안내를 확인 |
| 잘못된 SMILES 처리 | PASS | `유효한 SMILES로 해석할 수 없습니다.` 안내를 확인 |
| 유사도 기능 | PASS | 카페인 입력의 최근접 train 분자가 카페인이며 Tanimoto similarity가 `1.0000`임을 확인 |
| Streamlit 서버 시작 | PASS | 별도 포트에서 health endpoint와 첫 화면 모두 HTTP 200 응답 |
| Streamlit 실제 화면 | PASS | 브라우저에서 카페인 결과, 구조 이미지, 유사도, 빈 입력 오류와 잘못된 SMILES 오류를 확인 |
| Streamlit 서버 종료 | PASS | 검증 후 서버가 더 이상 실행 중이지 않음을 확인 |
| 임시 드라이브 연결 해제 | PASS | 설치 경로 단축에 사용한 임시 연결이 남지 않음을 확인 |
| 보호 파일 무결성 | PASS | notebook, `results/`, `models/`의 36개 파일을 SHA-256으로 대조했으며 불일치 0개 |
| 공개 제외 파일 검사 | PASS | 원본 `BBBP.csv`, 비밀정보 파일, `.venv`, 캐시가 Git 추적 대상에 없음을 확인 |
| 임시 검증 폴더 삭제 | 미완료 | 안전 정책이 재귀 삭제 명령을 실행 전에 차단함. 모든 임시 폴더는 Git에서 제외되는 `work/` 아래에만 있음 |

자동 smoke test의 실제 출력은 다음과 같았습니다.

```text
Smoke test 통과: 입력 오류 처리, 기존 점수 불변, 유사도 범위, Train 입력 1.0, 자기 자신 제외
CCO BBB+ 예측 점수: 1.0000
Morgan fingerprint 설정: radius=2, fpSize=2048
```

## 설치 중 발생한 문제와 해결

앞선 clean-copy 설치 시도는 프로젝트의 전체 경로가 길어서 Jupyter 관련 패키지 파일을 만드는 단계에서 Windows 경로 길이 문제로 실패했습니다. 이것은 프로젝트 코드나 패키지 버전 충돌 때문에 발생한 오류가 아닙니다.

실제 fresh clone 검증에서는 복제본을 Windows 임시 드라이브 문자 `R:`로 짧게 연결한 뒤 설치했습니다. 설치, `pip check`, smoke test, 예측 함수와 Streamlit 실제 화면 검증이 모두 통과했으며 검증 후 임시 드라이브 연결도 해제했습니다. 임시 연결이 다음 명령 창에 유지되지 않아, 설치와 검증 명령마다 연결하고 바로 해제하는 방식으로 처리했습니다.

일반 사용자는 저장소를 `C:\projects\bbbp`처럼 짧은 경로에 clone하면 같은 문제를 피할 수 있습니다.

## 새 환경에서 실행하는 최소 명령

저장소 루트에서 다음 순서로 실행합니다.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
python scripts\smoke_test.py
python -m streamlit run app.py
```

저장된 모델을 사용하는 smoke test와 앱 실행에는 BBBP 원본 CSV를 저장소에 복사할 필요가 없습니다. 모델을 처음부터 다시 만드는 `scripts/build_model.py`를 실행할 때만 공식 데이터 URL에 접근할 인터넷 연결이 필요합니다.

## 검증 범위와 남은 확인

- 이번 검증은 기존 모델을 다시 학습하거나 교체하지 않았습니다.
- notebook, `results/`, `models/`의 기존 파일은 수정하지 않았습니다.
- 첫 커밋을 대상으로 한 실제 fresh clone과 독립 가상환경 검증은 완료했습니다.
- 아직 GitHub 원격 저장소에는 올리지 않았으므로, 공개 후 GitHub URL에서 다시 clone하는 검증은 별도 단계입니다.
- 안전 정책 때문에 남은 검증 복사본과 새 가상환경은 모두 Git에서 제외되는 로컬 `work/` 폴더 아래에 있으며 공개 파일 수에 포함되지 않습니다.
