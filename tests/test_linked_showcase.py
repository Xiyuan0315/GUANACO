"""Exercise the gallery through Dash's real HTTP callback endpoints."""

import sys
import base64
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

DEMO_DIR = Path(__file__).resolve().parents[1] / "examples" / "linked_views"


@pytest.fixture(scope="module")
def gallery():
    sys.path.insert(0, str(DEMO_DIR))
    try:
        from showcase import create_app

        yield create_app()
    finally:
        sys.path.remove(str(DEMO_DIR))


def invoke(app, client, output, values, *, state=None, changed=None):
    callback = app.callback_map[output]
    outputs = callback["output"]
    if isinstance(outputs, list):
        outputs = [{"id": value.component_id, "property": value.component_property} for value in outputs]
    else:
        outputs = {"id": outputs.component_id, "property": outputs.component_property}
    response = client.post("/_dash-update-component", json={
        "output": output, "outputs": outputs,
        "inputs": [{**item, "value": values.get(f"{item['id']}.{item['property']}")}
                   for item in callback["inputs"]],
        "state": [{**item, "value": state} for item in callback["state"]],
        "changedPropIds": changed or list(values),
    })
    assert response.status_code == 200, response.get_data(as_text=True)
    return response.get_json()["response"]


def navigate(app, client, slug, resets=0):
    output = next(key for key in app.callback_map if "example-content.children" in key)
    return invoke(app, client, output, {
        "gallery-location.search": f"?demo={slug}", "reset-example.n_clicks": resets,
    })


def components(tree):
    if isinstance(tree, dict):
        if "props" in tree:
            yield tree
        for value in tree.values():
            yield from components(value)
    elif isinstance(tree, list):
        for value in tree:
            yield from components(value)


def component(tree, id):
    return next(node["props"] for node in components(tree) if node["props"].get("id") == id)


INTERACTIONS = [
    ("cells", "gallery-cells", "cells", "selectedData",
     {"points": [{"id": f"cell-{i:04d}"} for i in range(30)]}, ("dotplot", "heatmap")),
    ("rna-protein", "gallery-multiomics", "rna", "selectedData",
     {"points": [{"id": f"cell-{i:04d}"} for i in range(30)]}, ("protein",)),
    ("spatial", "gallery-spatial", "neighborhoods", "clickData",
     {"points": [{"x": "Pyramidal_layer", "y": "Hippocampus"}]}, ("locations", "cooccurrence")),
    ("pathways", "gallery-pathways", "pathways", "clickData",
     {"points": [{"customdata": ["T-cell activation"]}]}, ("genes",)),
    ("atac", "gallery-atac", "gene_overview", "clickData",
     {"points": [{"x": "MS4A1", "y": "B"}]}, ("peaks",)),
]


def test_gallery_serves_layout_dependencies_and_styles(gallery):
    client = gallery.server.test_client()
    for route in ("/", "/_dash-layout", "/_dash-dependencies", "/assets/gallery.css"):
        assert client.get(route).status_code == 200
    assert len(gallery.callback_map) == 13


def test_page_spinner_only_tracks_navigation_not_plot_updates(gallery):
    client = gallery.server.test_client()
    layout = client.get("/_dash-layout").get_json()
    loading = component(layout, "example-loading")
    assert loading["target_components"] == {"example-content": "children"}
    assert loading["delay_show"] == 250
    assert loading["show_initially"] is False
    assert gallery.config.update_title is None
    navigation = next(value for key, value in gallery.callback_map.items()
                      if "example-content.children" in key)
    assert navigation["inputs"] == [
        {"id": "gallery-location", "property": "search"},
        {"id": "reset-example", "property": "n_clicks"},
    ]


@pytest.mark.parametrize("slug,prefix,source,event,payload,targets", INTERACTIONS)
def test_each_link_updates_details_and_resets_without_leaking(
    gallery, slug, prefix, source, event, payload, targets,
):
    visitor_a = gallery.server.test_client()
    initial = navigate(gallery, visitor_a, slug)
    assert initial[f"nav-{slug}"]["className"] == "demo-link active"
    ids = [node["props"]["id"] for node in components(initial) if "id" in node["props"]]
    assert len(ids) == len(set(ids))
    assert [id for id in ids if id.endswith("-state")] == [f"{prefix}-state"]
    empty = component(initial, f"{prefix}-state")["data"]
    assert empty == {"sources": {}}

    selection = invoke(gallery, visitor_a, f"{prefix}-state.data", {
        f"{prefix}-view-{source}.{event}": payload,
    }, state=empty)[f"{prefix}-state"]["data"]
    assert selection["sources"][source]["members"]
    for target in targets:
        id = f"{prefix}-view-{target}"
        if slug == "rna-protein":
            # This target is now a browser restyle, not a server figure response.
            assert f"{id}.figure" not in gallery.callback_map
            entry = gallery.callback_map[f"{id}-highlight.data"]
            assert "callback" not in entry
            assert {item["property"] for item in entry["inputs"]} == {
                "clickData", "selectedData", "figure",
            }
            assert component(initial, id)["figure"]["data"][0]["ids"]
            continue
        updated = invoke(gallery, visitor_a, f"{id}.figure", {
            f"{prefix}-state.data": selection,
        })[id]["figure"]
        assert updated["data"]
        assert updated != component(initial, id)["figure"]
        if slug == "spatial" and target == "locations":
            original = component(initial, id)["figure"]
            assert updated["layout"]["images"] == original["layout"]["images"]
            for axis in ("xaxis", "yaxis"):
                assert updated["layout"][axis]["range"] == original["layout"][axis]["range"]
            assert {trace["name"] for trace in updated["data"]} == {"Hippocampus", "Pyramidal_layer"}
        # Clearing selection restores the full, original view.
        cleared = invoke(gallery, visitor_a, f"{id}.figure", {
            f"{prefix}-state.data": empty,
        })[id]["figure"]
        assert cleared == component(initial, id)["figure"]

    # A new visitor and a reset both receive pristine stores AND original figures,
    # even though the notebook runtime mirrors the most recently selected IDs.
    visitor_b = gallery.server.test_client()
    assert navigate(gallery, visitor_b, slug) == initial
    reset = navigate(gallery, visitor_a, slug, resets=1)
    assert component(reset, f"{prefix}-state")["data"] == empty
    assert reset["example-content"]["children"]["props"]["key"] == f"{slug}-1"
    for target in targets:
        id = f"{prefix}-view-{target}"
        assert component(reset, id)["figure"] == component(initial, id)["figure"]


def test_unknown_example_falls_back_without_reflecting_input(gallery):
    page = navigate(gallery, gallery.server.test_client(), "unknown")
    assert page["nav-cells"]["className"] == "demo-link active"


def test_spatial_gallery_preserves_original_matrix_and_tissue_overlay(gallery):
    page = navigate(gallery, gallery.server.test_client(), "spatial")
    spatial = ad.read_h5ad(DEMO_DIR / "data" / "visium_hne_spatial.h5ad")
    assert spatial.n_obs == 2688
    assert spatial.X.nnz == 0  # No expression matrix is shipped with this example.
    groups = list(spatial.obs["cluster"].cat.categories)
    assert len(groups) == 15
    heatmap = component(page, "gallery-spatial-view-neighborhoods")["figure"]["data"][0]
    assert heatmap["x"] == heatmap["y"] == groups
    z = heatmap["z"]
    if isinstance(z, dict):
        z = np.frombuffer(base64.b64decode(z["bdata"]), dtype=z["dtype"]).reshape(15, 15)
    np.testing.assert_allclose(z, spatial.uns["cluster_nhood_enrichment"]["zscore"])
    image = component(page, "gallery-spatial-view-locations")["figure"]["layout"]["images"][0]
    assert image["source"].startswith("data:image/png;base64,")
    assert image["layer"] == "below"
    original_image = spatial.uns["spatial"]["V1_Adult_Mouse_Brain"]["images"]["hires"]
    assert (image["sizey"], image["sizex"]) == original_image.shape[:2]


@pytest.mark.skipif(
    not Path("/Users/xiyuanzhang/Documents/GUANACO_v2/data/visium_hne_spatial.h5ad").is_file(),
    reason="Original notebook data unavailable",
)
def test_spatial_export_matches_original_notebook_data():
    original = ad.read_h5ad(
        "/Users/xiyuanzhang/Documents/GUANACO_v2/data/visium_hne_spatial.h5ad", backed="r",
    )
    try:
        excerpt = ad.read_h5ad(DEMO_DIR / "data" / "visium_hne_spatial.h5ad")
        assert excerpt.obs.equals(original.obs[["cluster"]])
        np.testing.assert_array_equal(excerpt.obsm["spatial"], original.obsm["spatial"])
        for key in ("cluster_nhood_enrichment", "cluster_co_occurrence"):
            for name, values in original.uns[key].items():
                np.testing.assert_array_equal(excerpt.uns[key][name], values)
        for library, metadata in original.uns["spatial"].items():
            assert excerpt.uns["spatial"][library]["scalefactors"] == metadata["scalefactors"]
            for name, image in metadata["images"].items():
                np.testing.assert_array_equal(excerpt.uns["spatial"][library]["images"][name], image)
    finally:
        original.file.close()


def test_synthetic_spatial_data_never_searches_local_files(monkeypatch):
    # This loader is also used by notebooks, so synthetic mode must be explicit.
    sys.path.insert(0, str(DEMO_DIR))
    try:
        import demo_data

        def forbid(*args, **kwargs):
            raise AssertionError("The public synthetic demo must not discover local files")

        monkeypatch.setattr(demo_data, "_first_existing", forbid)
        pairs, spatial, curves, source = demo_data.load_spatial_relationship_demo(synthetic=True)
        assert source is None
        assert spatial.n_obs == 400
        assert not pairs.empty and not curves.empty
    finally:
        sys.path.remove(str(DEMO_DIR))
