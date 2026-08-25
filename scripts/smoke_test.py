"""Day 5 입력 처리와 저장 모델 추론 경로의 빠른 회귀 검사."""

from __future__ import annotations

from pathlib import Path
import hashlib
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.features import MORGAN_FP_SIZE, MORGAN_RADIUS, featurize_smiles  # noqa: E402
from src.predict import DEFAULT_MODEL_PATH, load_model_bundle, predict_smiles  # noqa: E402
from src.similarity import (  # noqa: E402
    DEFAULT_SIMILARITY_REFERENCE_PATH,
    compute_leave_one_out_nearest_neighbors,
    load_similarity_reference,
)


EXPECTED_MODEL_SHA256 = "e0e6bd289cba4d98930106ca40b9ee38de505c9c9dd0775b5eee31673c94ec6c"
FIXED_EXAMPLE_SCORES = {
    "CCO": 1.0,
    "CC(=O)O": 0.9733333333333334,
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O": 0.97,
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    bundle = load_model_bundle(DEFAULT_MODEL_PATH)
    similarity_reference = load_similarity_reference(DEFAULT_SIMILARITY_REFERENCE_PATH)

    assert file_sha256(Path(DEFAULT_MODEL_PATH)) == EXPECTED_MODEL_SHA256
    assert bundle["reproduction"]["confusion_matrix"] == [[53, 40], [6, 294]]
    assert round(bundle["reproduction"]["metrics"]["roc_auc"], 4) == 0.8960

    for valid_smiles in ["CCO", "CC(=O)O"]:
        features = featurize_smiles(valid_smiles)
        assert features["is_valid"], valid_smiles
        assert features["fingerprint"].shape == (1, MORGAN_FP_SIZE)
        assert set(np.unique(features["fingerprint"])).issubset({0, 1})

        prediction = predict_smiles(valid_smiles, bundle=bundle)
        assert prediction["is_valid"], prediction
        assert 0.0 <= prediction["bbb_positive_score"] <= 1.0
        assert prediction["predicted_label"] in {0, 1}
        assert prediction["threshold"] == 0.5

    for fixed_smiles, expected_score in FIXED_EXAMPLE_SCORES.items():
        without_similarity = predict_smiles(fixed_smiles, bundle=bundle)
        with_similarity = predict_smiles(
            fixed_smiles,
            bundle=bundle,
            similarity_reference=similarity_reference,
        )
        assert without_similarity["bbb_positive_score"] == expected_score
        assert with_similarity["bbb_positive_score"] == expected_score
        assert with_similarity["similarity"] is not None
        assert 0.0 <= with_similarity["similarity"]["maximum_similarity"] <= 1.0

    included_train_smiles = similarity_reference["train_canonical_smiles"][0]
    included_result = predict_smiles(
        included_train_smiles,
        bundle=bundle,
        similarity_reference=similarity_reference,
    )
    assert included_result["is_valid"]
    assert included_result["similarity"]["maximum_similarity"] == 1.0

    nearest_scores, nearest_indices = compute_leave_one_out_nearest_neighbors(
        similarity_reference["train_fingerprint_matrix"]
    )
    assert np.all(nearest_indices != np.arange(len(nearest_indices)))
    assert np.all((nearest_scores >= 0.0) & (nearest_scores <= 1.0))
    assert np.array_equal(
        nearest_indices, similarity_reference["nearest_neighbor_indices"]
    )

    invalid_inputs = ["", "   ", "this_is_not_smiles", None, "C" * 5000]
    for invalid_input in invalid_inputs:
        result = predict_smiles(invalid_input, bundle=bundle)
        assert not result["is_valid"]
        assert result["bbb_positive_score"] is None
        assert result["error_message"]
        assert result["similarity"] is None

    repeated_first = predict_smiles("CCO", bundle=bundle)
    reloaded_bundle = load_model_bundle(DEFAULT_MODEL_PATH)
    repeated_second = predict_smiles("CCO", bundle=reloaded_bundle)
    assert repeated_first["predicted_label"] == repeated_second["predicted_label"]
    assert repeated_first["bbb_positive_score"] == repeated_second["bbb_positive_score"]
    assert bundle["reproduction"]["saved_model_reload_prediction_match"] is True

    assert bundle["fingerprint"]["radius"] == MORGAN_RADIUS == 2
    assert bundle["fingerprint"]["fp_size"] == MORGAN_FP_SIZE == 2048
    assert similarity_reference["fingerprint"]["radius"] == MORGAN_RADIUS
    assert similarity_reference["fingerprint"]["fp_size"] == MORGAN_FP_SIZE
    app_source = (PROJECT_ROOT / "app.py").read_text(encoding="utf-8")
    build_source = (PROJECT_ROOT / "scripts" / "build_model.py").read_text(encoding="utf-8")
    assert "from src.predict import" in app_source
    assert "from src.features import" in build_source

    print(
        "Smoke test 통과: 입력 오류 처리, 기존 점수 불변, 유사도 범위, "
        "Train 입력 1.0, 자기 자신 제외"
    )
    print(f"CCO BBB+ 예측 점수: {repeated_first['bbb_positive_score']:.4f}")
    print(f"공용 Morgan 설정: radius={MORGAN_RADIUS}, fpSize={MORGAN_FP_SIZE}")


if __name__ == "__main__":
    main()
