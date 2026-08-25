"""Day 3 Random Forest를 공식 BBBP 데이터부터 재현하고 검증 후 저장합니다."""

from __future__ import annotations

from pathlib import Path
import hashlib
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.features import (  # noqa: E402
    MORGAN_FP_SIZE,
    MORGAN_RADIUS,
    canonical_smiles,
    molecules_to_morgan_matrix,
    parse_smiles,
)


DATA_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv"
MODEL_PATH = PROJECT_ROOT / "models" / "bbbp_random_forest.joblib"
RANDOM_STATE = 42
TEST_SIZE = 0.20
THRESHOLD = 0.5
RF_PARAMETERS = {
    "n_estimators": 300,
    "random_state": RANDOM_STATE,
    "n_jobs": -1,
}
EXPECTED_METRICS = {"accuracy": 0.8830, "f1_score": 0.9274, "roc_auc": 0.8960}
EXPECTED_CONFUSION = np.array([[53, 40], [6, 294]])


def stable_hash(values) -> str:
    """행 순서를 확인할 수 있는 짧은 재현성 식별자를 만듭니다."""
    text = "\n".join(map(str, values))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_and_clean_bbbp():
    """Day 2와 같은 규칙으로 최종 모델용 1,965개를 만듭니다."""
    original_df = pd.read_csv(DATA_URL)
    original_df["source_index"] = original_df.index

    parsed = original_df["smiles"].map(parse_smiles)
    audit_df = original_df.copy()
    audit_df["mol"] = parsed.map(lambda item: item[0])
    audit_df["parse_error"] = parsed.map(lambda item: item[1])
    audit_df["parse_ok"] = audit_df["mol"].notna()

    candidate_df = audit_df.loc[audit_df["parse_ok"]].copy()
    candidate_df["canonical_smiles"] = candidate_df["mol"].map(canonical_smiles)
    label_counts = candidate_df.groupby("canonical_smiles")["p_np"].nunique()
    conflicting_keys = set(label_counts.loc[label_counts > 1].index)
    final_model_df = (
        candidate_df.loc[~candidate_df["canonical_smiles"].isin(conflicting_keys)]
        .sort_values("source_index")
        .drop_duplicates(subset="canonical_smiles", keep="first")
        .reset_index(drop=True)
    )

    parse_failures = int((~audit_df["parse_ok"]).sum())
    conflict_exclusions = int(candidate_df["canonical_smiles"].isin(conflicting_keys).sum())
    duplicate_exclusions = len(candidate_df) - conflict_exclusions - len(final_model_df)
    assert len(original_df) == 2050
    assert parse_failures == 11
    assert conflict_exclusions == 20
    assert duplicate_exclusions == 54
    assert len(final_model_df) == 1965
    assert len(original_df) == parse_failures + conflict_exclusions + duplicate_exclusions + len(final_model_df)
    return final_model_df


def reproduce_and_save_model(model_path: Path = MODEL_PATH) -> dict:
    """Day 3 test 결과가 일치할 때만 학습된 train-only 모델 bundle을 저장합니다."""
    final_model_df = load_and_clean_bbbp()
    train_df, test_df = train_test_split(
        final_model_df,
        test_size=TEST_SIZE,
        stratify=final_model_df["p_np"],
        random_state=RANDOM_STATE,
    )
    train_df = train_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)
    y_train = train_df["p_np"].to_numpy()
    y_test = test_df["p_np"].to_numpy()

    assert len(train_df) == 1572 and len(test_df) == 393
    assert train_df["p_np"].value_counts().sort_index().to_dict() == {0: 372, 1: 1200}
    assert test_df["p_np"].value_counts().sort_index().to_dict() == {0: 93, 1: 300}
    assert set(train_df["canonical_smiles"]).isdisjoint(set(test_df["canonical_smiles"]))

    X_train = molecules_to_morgan_matrix(train_df["mol"].tolist())
    X_test = molecules_to_morgan_matrix(test_df["mol"].tolist())
    assert X_train.shape == (1572, MORGAN_FP_SIZE)
    assert X_test.shape == (393, MORGAN_FP_SIZE)
    assert set(np.unique(X_train)).issubset({0, 1})
    assert set(np.unique(X_test)).issubset({0, 1})

    model = RandomForestClassifier(**RF_PARAMETERS)
    model.fit(X_train, y_train)
    prediction = model.predict(X_test)
    positive_column = int(np.where(model.classes_ == 1)[0][0])
    score = model.predict_proba(X_test)[:, positive_column]
    metrics = {
        "accuracy": accuracy_score(y_test, prediction),
        "f1_score": f1_score(y_test, prediction, pos_label=1, zero_division=0),
        "roc_auc": roc_auc_score(y_test, score),
    }
    matrix = confusion_matrix(y_test, prediction, labels=[0, 1])

    for metric_name, expected_value in EXPECTED_METRICS.items():
        assert round(metrics[metric_name], 4) == expected_value, (
            f"{metric_name} 재현 실패: {metrics[metric_name]:.6f} != {expected_value:.4f}"
        )
    assert np.array_equal(matrix, EXPECTED_CONFUSION), f"confusion matrix 재현 실패: {matrix.tolist()}"

    bundle = {
        "bundle_version": 1,
        "model": model,
        "feature_type": "Morgan bit vector",
        "fingerprint": {
            "radius": MORGAN_RADIUS,
            "fp_size": MORGAN_FP_SIZE,
            "generator": "RDKit FingerprintGenerator",
            "is_bit_vector": True,
        },
        "labels": {0: "BBB-", 1: "BBB+"},
        "threshold": THRESHOLD,
        "random_forest_parameters": model.get_params(deep=True),
        "training": {
            "train_molecules": len(train_df),
            "test_molecules": len(test_df),
            "test_in_training": False,
            "split": "stratified 80:20",
            "random_state": RANDOM_STATE,
            "train_source_order_sha256": stable_hash(train_df["source_index"]),
            "test_source_order_sha256": stable_hash(test_df["source_index"]),
        },
        "data_source": DATA_URL,
        "reproduction": {
            "matches_day3": True,
            "metrics": metrics,
            "confusion_matrix": matrix.tolist(),
            "saved_model_reload_prediction_match": True,
        },
        "score_note": "The BBB+ score is an uncalibrated Random Forest score, not a clinical probability.",
    }

    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, model_path, compress=3)

    loaded_bundle = joblib.load(model_path)
    loaded_model = loaded_bundle["model"]
    loaded_prediction = loaded_model.predict(X_test)
    loaded_positive_column = int(np.where(loaded_model.classes_ == 1)[0][0])
    loaded_score = loaded_model.predict_proba(X_test)[:, loaded_positive_column]
    assert np.array_equal(loaded_prediction, prediction)
    max_score_difference = float(np.max(np.abs(loaded_score - score)))
    assert np.allclose(loaded_score, score, rtol=0.0, atol=1e-12)
    assert round(roc_auc_score(y_test, loaded_score), 4) == EXPECTED_METRICS["roc_auc"]
    assert loaded_bundle["training"]["train_molecules"] == 1572
    assert loaded_bundle["training"]["test_molecules"] == 393

    return {
        "model_path": str(model_path),
        "model_bytes": model_path.stat().st_size,
        "train_rows": len(train_df),
        "test_rows": len(test_df),
        "metrics": metrics,
        "confusion_matrix": matrix.tolist(),
        "reload_prediction_match": True,
        "reload_max_score_difference": max_score_difference,
    }


def main():
    summary = reproduce_and_save_model()
    print("Day 3 Random Forest 재현 및 저장 검증 통과")
    print(f"Train/Test: {summary['train_rows']:,}/{summary['test_rows']:,}")
    print(
        "Test metrics: "
        f"Accuracy={summary['metrics']['accuracy']:.4f}, "
        f"F1={summary['metrics']['f1_score']:.4f}, "
        f"ROC-AUC={summary['metrics']['roc_auc']:.4f}"
    )
    print(f"Confusion matrix: {summary['confusion_matrix']}")
    print(f"저장 후 재로드 예측 일치: {summary['reload_prediction_match']}")
    print(f"재로드 점수 최대 절대차: {summary['reload_max_score_difference']:.3e}")
    print(f"모델 파일: {summary['model_path']} ({summary['model_bytes'] / 1024 / 1024:.2f} MiB)")


if __name__ == "__main__":
    main()
