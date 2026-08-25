"""저장된 Day 3 Random Forest bundle을 이용한 단일 SMILES 추론."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np

from src.features import MORGAN_FP_SIZE, MORGAN_RADIUS, featurize_smiles
from src.similarity import find_most_similar_train


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL_PATH = PROJECT_ROOT / "models" / "bbbp_random_forest.joblib"


def load_model_bundle(model_path: str | Path = DEFAULT_MODEL_PATH) -> dict:
    """joblib bundle을 읽고 학습·추론 fingerprint 설정이 같은지 검사합니다."""
    path = Path(model_path)
    if not path.exists():
        raise FileNotFoundError(
            f"저장된 모델이 없습니다: {path}. 먼저 python scripts/build_model.py를 실행하세요."
        )

    bundle = joblib.load(path)
    required_keys = {"model", "fingerprint", "labels", "threshold", "training", "reproduction"}
    missing_keys = required_keys.difference(bundle)
    if missing_keys:
        raise ValueError(f"모델 bundle에 필요한 정보가 없습니다: {sorted(missing_keys)}")
    if bundle["fingerprint"]["radius"] != MORGAN_RADIUS:
        raise ValueError("저장 모델과 현재 코드의 Morgan radius가 다릅니다.")
    if bundle["fingerprint"]["fp_size"] != MORGAN_FP_SIZE:
        raise ValueError("저장 모델과 현재 코드의 fingerprint 크기가 다릅니다.")
    if float(bundle["threshold"]) != 0.5:
        raise ValueError("저장 모델의 threshold가 Day 3 기준 0.5와 다릅니다.")
    return bundle


def predict_smiles(
    smiles: Any,
    bundle: dict | None = None,
    model_path: str | Path = DEFAULT_MODEL_PATH,
    similarity_reference: dict | None = None,
) -> dict:
    """SMILES 하나의 BBB+ 점수와 0.5 기준 예측을 안전한 사전 형태로 반환합니다."""
    result = {
        "original_smiles": smiles,
        "canonical_smiles": None,
        "predicted_label": None,
        "predicted_class": None,
        "bbb_positive_score": None,
        "threshold": 0.5,
        "is_valid": False,
        "error_message": "",
        "similarity": None,
    }

    features = featurize_smiles(smiles)
    if not features["is_valid"]:
        result["error_message"] = features["error_message"]
        return result

    try:
        active_bundle = bundle if bundle is not None else load_model_bundle(model_path)
        threshold = float(active_bundle["threshold"])
        model = active_bundle["model"]
        positive_column = int(np.where(model.classes_ == 1)[0][0])
        score = float(model.predict_proba(features["fingerprint"])[0, positive_column])
        # scikit-learn predict와 같게 0.5 동점은 클래스 0(BBB-)으로 둡니다.
        predicted_label = int(score > threshold)
    except Exception as exc:
        result["error_message"] = f"예측을 수행할 수 없습니다: {exc}"
        return result

    similarity = None
    if similarity_reference is not None:
        try:
            similarity = find_most_similar_train(
                features["fingerprint"], similarity_reference
            )
        except Exception as exc:
            result["error_message"] = f"구조 유사도를 계산할 수 없습니다: {exc}"
            return result

    result.update(
        {
            "canonical_smiles": features["canonical_smiles"],
            "predicted_label": predicted_label,
            "predicted_class": active_bundle["labels"][predicted_label],
            "bbb_positive_score": score,
            "threshold": threshold,
            "is_valid": True,
            "error_message": "",
            "similarity": similarity,
        }
    )
    return result
