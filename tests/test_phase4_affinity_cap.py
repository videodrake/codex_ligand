from egfr_pipeline.phase4.perturbation_scoring import (
    AFFINITY_CAP_FRACTION,
    apply_affinity_cap,
)


def test_affinity_cap_enforces_fraction_and_guards():
    row = {
        "mechanistic_class": "ligandable_but_ppi_irrelevant_candidate",
        "A1_ppi_interface": "0.1",
        "A2_druggability": "0.9",
        "A3_perturbation_relevance": "0.8",
        "A4_state_robustness": "0.1",
    }
    raw_score = 0.1 * 0.30 + 0.9 * 0.25 + 0.8 * 0.30 + 0.1 * 0.15

    capped = apply_affinity_cap(row, raw_score)
    bio_component = 0.1 * 0.30 + 0.1 * 0.15
    capped_affinity = capped - bio_component

    assert capped < raw_score
    assert capped_affinity / capped <= AFFINITY_CAP_FRACTION + 1e-9
    assert apply_affinity_cap(row, raw_score, cap_fraction=0.0) == round(bio_component, 4)
    assert apply_affinity_cap(row, raw_score, cap_fraction=1.0) == round(raw_score, 4)
