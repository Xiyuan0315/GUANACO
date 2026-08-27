from guanaco.utils.search import ranked_substring_matches


def test_exact_feature_match_is_not_truncated_by_earlier_substring_matches():
    genes = [f"C2CD4D{i}" for i in range(12)] + ["CD48", "CD4", "CD40LG"]

    matches = ranked_substring_matches(genes, "cd4", limit=10)

    assert matches[0] == "CD4"
    assert len(matches) == 10


def test_display_labels_can_be_ranked_by_underlying_feature_names():
    labels = ["RNA · CD48", "ADT · CD4", "RNA · CD40LG"]
    features = ["CD48", "CD4", "CD40LG"]

    matches = ranked_substring_matches(
        labels,
        "CD4",
        match_values=features,
    )

    assert matches == ["ADT · CD4", "RNA · CD48", "RNA · CD40LG"]
