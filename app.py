"""저장된 Day 3 Random Forest를 사용하는 단일 페이지 Streamlit 데모."""

from pathlib import Path

import streamlit as st
from rdkit.Chem import Draw

from src.features import parse_smiles
from src.predict import DEFAULT_MODEL_PATH, load_model_bundle, predict_smiles
from src.similarity import (
    DEFAULT_SIMILARITY_REFERENCE_PATH,
    load_similarity_reference,
)


st.set_page_config(page_title="BBBP SMILES 예측 데모", page_icon="🧠", layout="centered")


@st.cache_resource
def get_model_bundle():
    """앱 실행 중 모델을 한 번만 디스크에서 읽습니다."""
    return load_model_bundle(DEFAULT_MODEL_PATH)


@st.cache_resource
def get_similarity_reference():
    """앱 실행 중 Train 유사도 참조자료를 한 번만 디스크에서 읽습니다."""
    return load_similarity_reference(DEFAULT_SIMILARITY_REFERENCE_PATH)


st.title("BBBP SMILES 예측 데모")
st.write("SMILES 하나를 입력하면 Day 3 Random Forest가 BBB+ 예측 점수를 계산합니다.")

examples = {
    "에탄올 예시": "CCO",
    "아세트산 예시": "CC(=O)O",
    "카페인 예시": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "직접 입력": "",
}
selected_example = st.selectbox("입력 예시 선택", list(examples))
smiles_input = st.text_input(
    "SMILES",
    value=examples[selected_example],
    placeholder="예: CCO",
    help="버튼을 누르기 전에는 예측하지 않습니다.",
)

if st.button("Predict", type="primary"):
    try:
        bundle = get_model_bundle()
        similarity_reference = get_similarity_reference()
        prediction = predict_smiles(
            smiles_input,
            bundle=bundle,
            similarity_reference=similarity_reference,
        )
    except FileNotFoundError:
        st.error("저장된 모델이 없습니다. 먼저 `python scripts/build_model.py`를 실행해 주세요.")
    except Exception as exc:
        st.error(f"모델을 불러올 수 없습니다: {exc}")
    else:
        if not prediction["is_valid"]:
            st.error(prediction["error_message"] or "유효한 SMILES로 해석할 수 없습니다.")
        else:
            mol, _ = parse_smiles(smiles_input)
            structure_image = Draw.MolToImage(mol, size=(500, 350))
            st.image(structure_image, caption="입력 분자의 RDKit 2차원 구조")

            st.subheader("예측 결과")
            left, right = st.columns(2)
            left.metric("예측 분류", prediction["predicted_class"])
            right.metric("BBB+ 예측 점수", f"{prediction['bbb_positive_score']:.3f}")
            st.write(f"**Canonical SMILES:** `{prediction['canonical_smiles']}`")
            st.write(f"**사용한 결정 기준:** {prediction['threshold']:.1f}")
            st.caption(
                "BBB+ 예측 점수는 보정(calibration)하지 않은 Random Forest 출력입니다. "
                "실제 BBB 투과 확률, 뇌 농도 또는 logBB가 아닙니다."
            )

            similarity = prediction["similarity"]
            st.subheader("학습 데이터와의 구조 유사도")
            st.metric("가장 높은 Tanimoto similarity", f"{similarity['maximum_similarity']:.3f}")
            st.write(f"**가장 유사한 Train 분자:** {similarity['train_name']}")
            st.write(
                f"**Train 원본 행 번호:** {similarity['train_source_index']} "
                "(0부터 시작)"
            )
            st.write(
                "**가장 유사한 Train canonical SMILES:** "
                f"`{similarity['train_canonical_smiles']}`"
            )

            similar_mol, _ = parse_smiles(similarity["train_canonical_smiles"])
            if similar_mol is not None:
                similar_image = Draw.MolToImage(similar_mol, size=(500, 350))
                st.image(similar_image, caption="가장 유사한 Train 분자의 RDKit 2차원 구조")

            if similarity["below_reference_range"]:
                st.warning(
                    "입력 분자는 학습 데이터의 일반적인 구조 유사도 범위보다 낮습니다. "
                    "예측을 주의해서 해석하세요."
                )
            st.caption(
                "Tanimoto similarity는 Morgan fingerprint로 비교한 구조 유사도이며, "
                "예측 신뢰확률이나 정확도 보장이 아닙니다. 낮은 유사도 경고 기준 "
                f"{similarity['warning_threshold']:.3f}은 Train 각 분자에서 자기 자신을 "
                "제외한 최근접 유사도 분포의 10번째 백분위수입니다. 이 기준은 경험적이고 "
                "현재 데이터에 의존합니다."
            )

st.warning(
    "연구 및 교육용 데모입니다. 실제 약물개발·임상·의학적 판단에 사용할 수 없습니다. "
    "고정 test에서 BBB- Specificity는 0.5699로 BBB+ Sensitivity 0.9800보다 낮았습니다. "
    "학습 데이터와 구조가 크게 다른 분자의 결과는 특히 신뢰하기 어렵습니다."
)
st.caption(
    "모델: Train 1,572개로 학습한 Day 3 Morgan fingerprint Random Forest · "
    "threshold=0.5 · 앱에서는 재학습하지 않음"
)
