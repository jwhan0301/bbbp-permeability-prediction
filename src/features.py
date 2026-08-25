"""학습과 추론에서 함께 사용하는 SMILES 및 Morgan fingerprint 함수."""

from __future__ import annotations

from typing import Any

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors


# Day 3에서 사용한 설정입니다. 학습 스크립트와 앱은 이 값을 직접 가져옵니다.
MORGAN_RADIUS = 2
MORGAN_FP_SIZE = 2048
MAX_SMILES_LENGTH = 4096
DESCRIPTOR_NAMES = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "FractionCSP3",
]


def parse_smiles(smiles: Any, max_length: int = MAX_SMILES_LENGTH):
    """SMILES를 검사해 ``(RDKit Mol, 오류 메시지)``를 반환합니다."""
    if smiles is None:
        return None, "SMILES가 입력되지 않았습니다."
    if not isinstance(smiles, str):
        return None, "SMILES는 문자열로 입력해야 합니다."

    cleaned = smiles.strip()
    if not cleaned:
        return None, "SMILES가 비어 있습니다."
    if len(cleaned) > max_length:
        return None, f"SMILES가 너무 깁니다. {max_length:,}자 이하로 입력해 주세요."

    # 잘못된 입력의 자세한 RDKit 로그는 UI에 노출하지 않고 이 호출 동안만 막습니다.
    log_blocker = rdBase.BlockLogs()
    try:
        mol = Chem.MolFromSmiles(cleaned)
    except Exception:
        mol = None
    finally:
        del log_blocker

    if mol is None:
        return None, "유효한 SMILES로 해석할 수 없습니다."
    return mol, ""


def canonical_smiles(mol: Chem.Mol) -> str:
    """RDKit 분자에서 입체화학 정보를 가능한 범위에서 유지한 canonical SMILES를 만듭니다."""
    if mol is None:
        raise ValueError("canonical SMILES를 만들 RDKit 분자가 없습니다.")
    value = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    if not value:
        raise ValueError("canonical SMILES 생성에 실패했습니다.")
    return value


def make_morgan_generator():
    """Day 3와 같은 Morgan bit-vector generator를 만듭니다."""
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS,
        fpSize=MORGAN_FP_SIZE,
    )


def generate_morgan_fingerprint(mol: Chem.Mol):
    """RDKit 분자 한 개를 Day 3 설정의 Morgan bit vector로 바꿉니다."""
    if mol is None:
        raise ValueError("fingerprint를 만들 RDKit 분자가 없습니다.")
    return make_morgan_generator().GetFingerprint(mol)


def fingerprint_to_numpy(fingerprint) -> np.ndarray:
    """Morgan fingerprint를 모델 입력용 ``(1, 2048)`` NumPy 배열로 바꿉니다."""
    row = np.zeros(MORGAN_FP_SIZE, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fingerprint, row)
    return row.reshape(1, -1)


def molecules_to_morgan_matrix(molecules) -> np.ndarray:
    """여러 RDKit 분자를 ``(분자 수, 2048)`` 행렬로 묶습니다."""
    rows = [fingerprint_to_numpy(generate_morgan_fingerprint(mol))[0] for mol in molecules]
    if not rows:
        return np.empty((0, MORGAN_FP_SIZE), dtype=np.uint8)
    return np.asarray(rows, dtype=np.uint8)


def calculate_descriptors(mol: Chem.Mol) -> dict:
    """Day 2와 같은 순서로 한 분자의 8개 descriptor를 계산합니다."""
    if mol is None:
        raise ValueError("descriptor를 계산할 RDKit 분자가 없습니다.")
    return {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "NumHDonors": Lipinski.NumHDonors(mol),
        "NumHAcceptors": Lipinski.NumHAcceptors(mol),
        "NumRotatableBonds": Lipinski.NumRotatableBonds(mol),
        "RingCount": Lipinski.RingCount(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
    }


def molecules_to_descriptor_matrix(molecules) -> np.ndarray:
    """여러 분자를 Day 2 순서의 ``(분자 수, 8)`` 숫자 행렬로 바꿉니다."""
    rows = [
        [descriptor_values[name] for name in DESCRIPTOR_NAMES]
        for descriptor_values in (calculate_descriptors(mol) for mol in molecules)
    ]
    if not rows:
        return np.empty((0, len(DESCRIPTOR_NAMES)), dtype=float)
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape[1] != len(DESCRIPTOR_NAMES):
        raise ValueError("descriptor 열 수가 8개가 아닙니다.")
    return matrix


def featurize_smiles(smiles: Any) -> dict:
    """검증부터 canonical SMILES와 Morgan 배열 생성까지 한 번에 수행합니다."""
    mol, error_message = parse_smiles(smiles)
    if mol is None:
        return {
            "is_valid": False,
            "mol": None,
            "canonical_smiles": None,
            "fingerprint": None,
            "error_message": error_message,
        }

    try:
        canonical = canonical_smiles(mol)
        fingerprint = fingerprint_to_numpy(generate_morgan_fingerprint(mol))
    except Exception as exc:
        return {
            "is_valid": False,
            "mol": None,
            "canonical_smiles": None,
            "fingerprint": None,
            "error_message": f"분자 특징을 만들 수 없습니다: {exc}",
        }

    return {
        "is_valid": True,
        "mol": mol,
        "canonical_smiles": canonical,
        "fingerprint": fingerprint,
        "error_message": "",
    }
