"""Train Morgan fingerprint와 입력 분자의 Tanimoto 구조 유사도를 계산합니다."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np

from src.features import MORGAN_FP_SIZE, MORGAN_RADIUS


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SIMILARITY_REFERENCE_PATH = (
    PROJECT_ROOT / "models" / "bbbp_train_similarity_reference.joblib"
)


def compute_leave_one_out_nearest_neighbors(fingerprint_matrix: np.ndarray):
    """각 Train 행에서 자기 자신을 제외한 가장 가까운 행과 유사도를 구합니다."""
    matrix = np.asarray(fingerprint_matrix, dtype=np.uint8)
    if matrix.ndim != 2 or matrix.shape[1] != MORGAN_FP_SIZE:
        raise ValueError(f"fingerprint 행렬은 (n, {MORGAN_FP_SIZE}) shape이어야 합니다.")
    if len(matrix) < 2:
        raise ValueError("자기 자신을 제외하려면 Train fingerprint가 2개 이상 필요합니다.")
    if not set(np.unique(matrix)).issubset({0, 1}):
        raise ValueError("fingerprint 행렬은 0과 1로만 구성되어야 합니다.")

    integer_matrix = matrix.astype(np.uint16, copy=False)
    intersections = integer_matrix @ integer_matrix.T
    bit_counts = integer_matrix.sum(axis=1, dtype=np.int32)
    unions = bit_counts[:, None] + bit_counts[None, :] - intersections
    similarities = np.divide(
        intersections,
        unions,
        out=np.zeros(intersections.shape, dtype=float),
        where=unions > 0,
    )
    np.fill_diagonal(similarities, -1.0)
    nearest_indices = similarities.argmax(axis=1)
    nearest_scores = similarities[np.arange(len(matrix)), nearest_indices]

    assert np.all(nearest_indices != np.arange(len(matrix)))
    assert np.all((nearest_scores >= 0.0) & (nearest_scores <= 1.0))
    return nearest_scores, nearest_indices


def create_similarity_reference(train_df, fingerprint_matrix: np.ndarray) -> dict:
    """앱에 필요한 Train 식별정보, fingerprint와 경험적 경고 기준을 묶습니다."""
    matrix = np.asarray(fingerprint_matrix, dtype=np.uint8)
    required_columns = {"source_index", "canonical_smiles"}
    missing_columns = required_columns.difference(train_df.columns)
    if missing_columns:
        raise ValueError(f"Train 표에 필요한 열이 없습니다: {sorted(missing_columns)}")
    if matrix.shape != (len(train_df), MORGAN_FP_SIZE):
        raise ValueError("Train 행 수와 fingerprint 행렬 shape이 일치하지 않습니다.")

    names = []
    for _, row in train_df.iterrows():
        raw_name = row.get("name", "")
        raw_number = row.get("num", "")
        if isinstance(raw_name, str) and raw_name.strip():
            names.append(raw_name.strip())
        elif str(raw_number).strip() and str(raw_number).lower() != "nan":
            names.append(f"BBBP 데이터 번호 {raw_number}")
        else:
            names.append(f"원본 행 {int(row['source_index'])}")

    nearest_scores, nearest_indices = compute_leave_one_out_nearest_neighbors(matrix)
    distribution = {
        "min": float(np.min(nearest_scores)),
        "p10": float(np.percentile(nearest_scores, 10)),
        "median": float(np.median(nearest_scores)),
        "p90": float(np.percentile(nearest_scores, 90)),
        "max": float(np.max(nearest_scores)),
    }
    return {
        "artifact_version": 1,
        "fingerprint": {"radius": MORGAN_RADIUS, "fp_size": MORGAN_FP_SIZE},
        "training": {"train_molecules": int(len(train_df)), "random_state": 42},
        "train_names": names,
        "train_source_indices": train_df["source_index"].to_numpy(dtype=np.int64),
        "train_canonical_smiles": train_df["canonical_smiles"].astype(str).tolist(),
        "train_fingerprint_matrix": matrix,
        "nearest_neighbor_indices": nearest_indices.astype(np.int64),
        "nearest_neighbor_similarities": nearest_scores.astype(float),
        "nearest_neighbor_distribution": distribution,
        "warning_threshold": distribution["p10"],
        "warning_rule": (
            "입력의 최대 Train Tanimoto similarity가 Train 내부 자기 제외 최근접 "
            "유사도 분포의 10번째 백분위수보다 낮을 때 경고"
        ),
    }


def save_similarity_reference(
    reference: dict,
    path: str | Path = DEFAULT_SIMILARITY_REFERENCE_PATH,
) -> Path:
    """유사도 참조자료를 압축된 작은 joblib 파일로 저장합니다."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(reference, output_path, compress=3)
    return output_path


def load_similarity_reference(
    path: str | Path = DEFAULT_SIMILARITY_REFERENCE_PATH,
) -> dict:
    """유사도 참조자료를 읽고 학습·예측 Morgan 설정과 일치하는지 확인합니다."""
    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(
            f"유사도 참조자료가 없습니다: {input_path}. "
            "먼저 python scripts/build_similarity_reference.py를 실행하세요."
        )
    reference = joblib.load(input_path)
    required_keys = {
        "fingerprint",
        "train_names",
        "train_source_indices",
        "train_canonical_smiles",
        "train_fingerprint_matrix",
        "nearest_neighbor_indices",
        "nearest_neighbor_similarities",
        "nearest_neighbor_distribution",
        "warning_threshold",
    }
    missing_keys = required_keys.difference(reference)
    if missing_keys:
        raise ValueError(f"유사도 참조자료에 필요한 정보가 없습니다: {sorted(missing_keys)}")
    if reference["fingerprint"]["radius"] != MORGAN_RADIUS:
        raise ValueError("유사도 참조자료와 현재 코드의 Morgan radius가 다릅니다.")
    if reference["fingerprint"]["fp_size"] != MORGAN_FP_SIZE:
        raise ValueError("유사도 참조자료와 현재 코드의 fingerprint 크기가 다릅니다.")
    matrix = np.asarray(reference["train_fingerprint_matrix"])
    if matrix.shape != (len(reference["train_canonical_smiles"]), MORGAN_FP_SIZE):
        raise ValueError("유사도 참조자료의 행 수 또는 fingerprint shape이 잘못되었습니다.")
    if len(reference["nearest_neighbor_indices"]) != len(matrix):
        raise ValueError("자기 제외 최근접 이웃 결과의 행 수가 잘못되었습니다.")
    return reference


def find_most_similar_train(query_fingerprint: np.ndarray, reference: dict) -> dict:
    """입력 fingerprint와 가장 유사한 Train 분자 한 개를 반환합니다."""
    query = np.asarray(query_fingerprint, dtype=np.uint8).reshape(-1)
    matrix = np.asarray(reference["train_fingerprint_matrix"], dtype=np.uint8)
    if query.shape != (MORGAN_FP_SIZE,):
        raise ValueError(f"입력 fingerprint는 {MORGAN_FP_SIZE}개 bit여야 합니다.")
    if matrix.ndim != 2 or matrix.shape[1] != MORGAN_FP_SIZE:
        raise ValueError("Train fingerprint 행렬 shape이 잘못되었습니다.")
    if reference["fingerprint"]["radius"] != MORGAN_RADIUS:
        raise ValueError("학습과 유사도 계산의 Morgan radius가 다릅니다.")
    if reference["fingerprint"]["fp_size"] != MORGAN_FP_SIZE:
        raise ValueError("학습과 유사도 계산의 fingerprint 크기가 다릅니다.")

    intersections = np.count_nonzero(matrix & query, axis=1)
    unions = np.count_nonzero(matrix | query, axis=1)
    similarities = np.divide(
        intersections,
        unions,
        out=np.zeros(len(matrix), dtype=float),
        where=unions > 0,
    )
    best_index = int(np.argmax(similarities))
    maximum_similarity = float(similarities[best_index])
    warning_threshold = float(reference["warning_threshold"])
    assert 0.0 <= maximum_similarity <= 1.0
    return {
        "train_position": best_index,
        "train_name": reference["train_names"][best_index],
        "train_source_index": int(reference["train_source_indices"][best_index]),
        "train_canonical_smiles": reference["train_canonical_smiles"][best_index],
        "maximum_similarity": maximum_similarity,
        "warning_threshold": warning_threshold,
        "below_reference_range": maximum_similarity < warning_threshold,
    }
