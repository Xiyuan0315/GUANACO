"""Unit tests for the pure helper functions behind the notebook (`gc.pl`) and
marimo control-bar APIs. These cover input coercion and basis/option resolution
in isolation -- no plotting, no AnnData heavy lifting."""

from types import SimpleNamespace

from guanaco import marimo as gm
from guanaco import widget as gw


def _adata_with_obsm(*keys):
    """Minimal stand-in exposing only the ``.obsm`` mapping these helpers read."""
    return SimpleNamespace(obsm={k: None for k in keys})


# --------------------------------------------------------------------------- #
# widget._as_list
# --------------------------------------------------------------------------- #
def test_as_list_wraps_single_string_in_one_element_list():
    assert gw._as_list("CD8A") == ["CD8A"]


def test_as_list_maps_none_to_empty_list():
    assert gw._as_list(None) == []


def test_as_list_copies_sequence_into_a_new_list():
    source = ("CD8A", "CD4")
    result = gw._as_list(source)
    assert result == ["CD8A", "CD4"]
    assert isinstance(result, list)


# --------------------------------------------------------------------------- #
# widget._resolve_basis
# --------------------------------------------------------------------------- #
def test_resolve_basis_returns_exact_obsm_key_when_present():
    adata = _adata_with_obsm("X_umap")
    assert gw._resolve_basis(adata, "X_umap") == "X_umap"


def test_resolve_basis_prefixes_scanpy_style_short_name_with_x():
    adata = _adata_with_obsm("X_umap")
    assert gw._resolve_basis(adata, "umap") == "X_umap"


def test_resolve_basis_passes_unknown_basis_through_unchanged():
    adata = _adata_with_obsm("X_umap")
    assert gw._resolve_basis(adata, "tsne") == "tsne"


# --------------------------------------------------------------------------- #
# marimo._coerce
# --------------------------------------------------------------------------- #
def test_coerce_keeps_value_when_it_is_a_valid_option():
    assert gm._coerce("b", ["a", "b", "c"]) == "b"


def test_coerce_uses_fallback_when_value_invalid_and_fallback_valid():
    assert gm._coerce("z", ["a", "b"], fallback="a") == "a"


def test_coerce_uses_first_option_when_value_and_fallback_both_invalid():
    assert gm._coerce("z", ["a", "b"], fallback="zz") == "a"


def test_coerce_returns_none_for_empty_options():
    assert gm._coerce("z", []) is None


# --------------------------------------------------------------------------- #
# marimo._default_basis
# --------------------------------------------------------------------------- #
def test_default_basis_prefers_umap_over_other_embeddings():
    adata = _adata_with_obsm("X_pca", "X_umap", "X_tsne")
    assert gm._default_basis(adata) == "X_umap"


def test_default_basis_falls_back_to_first_embedding_when_no_preference_matches():
    adata = _adata_with_obsm("X_custom", "X_other")
    assert gm._default_basis(adata) == "X_custom"


def test_default_basis_returns_none_when_no_embeddings_exist():
    assert gm._default_basis(_adata_with_obsm()) is None
