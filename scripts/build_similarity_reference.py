"""Day 3 Train 분자를 재현해 앱용 구조 유사도 참조자료를 만듭니다."""

from __future__ import annotations

from pathlib import Path
import hashlib
import sys

import numpy as np
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_model import (  # noqa: E402
    RANDOM_STATE,
    TEST_SIZE,
    load_and_clean_bbbp,
    stable_hash,
)
from src.features import MORGAN_FP_SIZE, MORGAN_RADIUS, molecules_to_morgan_matrix  # noqa: E402
from src.similarity import (  # noqa: E402
    DEFAULT_SIMILARITY_REFERENCE_PATH,
    create_similarity_reference,
    load_similarity_reference,
    save_similarity_reference,
)


MODEL_PATH = PROJECT_ROOT / "models" / "bbbp_random_forest.joblib"
EXPECTED_MODEL_SHA256 = "e0e6bd289cba4d98930106ca40b9ee38de505c9c9dd0775b5eee31673c94ec6c"
EXPECTED_TRAIN_SOURCE_HASH = "7a21d868677d85e44b64087819ec59c7a4979545193463ae36bcd97a11c5a38a"


def file_sha256(path: Path) -> str:
    """파일 내용을 바꾸지 않고 SHA-256을 계산합니다."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_similarity_reference(output_path: Path = DEFAULT_SIMILARITY_REFERENCE_PATH) -> dict:
    """공식 데이터부터 동일 Train을 재현해 참조자료를 저장하고 다시 검사합니다."""
    if file_sha256(MODEL_PATH) != EXPECTED_MODEL_SHA256:
        raise ValueError("기존 Day 3 앱 모델 파일의 SHA-256이 변경되었습니다.")

    final_model_df = load_and_clean_bbbp()
    train_df, test_df = train_test_split(
        final_model_df,
        test_size=TEST_SIZE,
        stratify=final_model_df["p_np"],
        random_state=RANDOM_STATE,
    )
    train_df = train_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)
    assert len(train_df) == 1572 and len(test_df) == 393
    assert set(train_df["source_index"]).isdisjoint(set(test_df["source_index"]))
    assert stable_hash(train_df["source_index"]) == EXPECTED_TRAIN_SOURCE_HASH

    fingerprint_matrix = molecules_to_morgan_matrix(train_df["mol"].tolist())
    assert fingerprint_matrix.shape == (1572, MORGAN_FP_SIZE)
    assert MORGAN_RADIUS == 2 and MORGAN_FP_SIZE == 2048
    assert set(np.unique(fingerprint_matrix)).issubset({0, 1})

    reference = create_similarity_reference(train_df, fingerprint_matrix)
    reference["training"]["train_source_order_sha256"] = stable_hash(
        train_df["source_index"]
    )
    saved_path = save_similarity_reference(reference, output_path)
    reloaded = load_similarity_reference(saved_path)
    assert np.array_equal(
        reloaded["train_fingerprint_matrix"], reference["train_fingerprint_matrix"]
    )
    assert np.all(
        reloaded["nearest_neighbor_indices"] != np.arange(len(train_df))
    )
    return {
        "path": str(saved_path),
        "bytes": saved_path.stat().st_size,
        "train_rows": len(train_df),
        "distribution": reloaded["nearest_neighbor_distribution"],
    }


def main():
    summary = build_similarity_reference()
    values = summary["distribution"]
    print("Train 구조 유사도 참조자료 생성 및 재로드 검증 통과")
    print(f"Train 행 수: {summary['train_rows']:,}")
    print(f"Morgan 설정: radius={MORGAN_RADIUS}, fpSize={MORGAN_FP_SIZE}")
    print(
        "자기 제외 최근접 Tanimoto 분포: "
        f"min={values['min']:.4f}, p10={values['p10']:.4f}, "
        f"median={values['median']:.4f}, p90={values['p90']:.4f}, "
        f"max={values['max']:.4f}"
    )
    print(f"저장 파일: {summary['path']} ({summary['bytes'] / 1024 / 1024:.2f} MiB)")


if __name__ == "__main__":
    main()
