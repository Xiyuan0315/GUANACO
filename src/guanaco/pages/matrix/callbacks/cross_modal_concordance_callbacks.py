"""Callbacks for the paired cross-modal concordance exploratory plot."""

from collections import OrderedDict
from threading import Lock

import numpy as np
from dash import Input, Output, State, no_update

from guanaco.pages.matrix.plots.cross_modal_concordance import (
    GROUP_CORRELATION,
    RELATIVE_SKEW,
    analyze_concordance,
    build_concordance_embedding,
    build_feature_scatter,
    build_group_summary_heatmap,
    calculate_group_correlations,
    empty_concordance_figure,
)
from guanaco.pages.visualizations.registry import EXPLORATION_WORKSPACE


_ACTIVE_TAB = "cross-modal-concordance-tab"
_VIEW_MODES = {RELATIVE_SKEW, GROUP_CORRELATION}


_LINKED_SELECTION_JS = """
function(leftSelection, rightSelection, comparison, mainView, filteredData) {
    const callbackContext = window.dash_clientside.callback_context;
    const triggered = (callbackContext.triggered || [])[0];
    if (!triggered) return null;
    const prop = triggered.prop_id || '';

    if (prop.endsWith('concordance-committed-features.data') ||
        prop.endsWith('concordance-main-view-trigger.data') ||
        prop.endsWith('global-filtered-data.data')) {
        return null;
    }

    let selection = null;
    if (prop.endsWith('concordance-embedding-plot.selectedData')) {
        selection = leftSelection;
    } else if (prop.endsWith('concordance-feature-plot.selectedData')) {
        selection = rightSelection;
    } else {
        return window.dash_clientside.no_update;
    }

    if (!selection || !selection.points || selection.points.length === 0) {
        return null;
    }
    const seen = new Set();
    const cellIds = [];
    for (const point of selection.points) {
        const custom = point.customdata;
        const cellId = Array.isArray(custom) ? custom[0] : custom;
        if (cellId == null) continue;
        const key = String(cellId);
        if (!seen.has(key)) {
            seen.add(key);
            cellIds.push(key);
        }
    }
    return cellIds.length ? {cell_ids: cellIds} : null;
}
"""


_LINKED_HIGHLIGHT_JS = """
function(selection, _leftFigure, _rightFigure) {
    const noUpdate = window.dash_clientside.no_update;
    const ids = selection && selection.cell_ids ? new Set(selection.cell_ids) : null;

    function highlight(graphId) {
        const wrapper = document.getElementById(graphId);
        if (!wrapper || !window.Plotly) return;
        const graph = wrapper.classList.contains('js-plotly-plot')
            ? wrapper
            : wrapper.querySelector('.js-plotly-plot');
        if (!graph || !graph.data) return;

        const traceIndices = [];
        const selections = [];
        for (let traceIndex = 0; traceIndex < graph.data.length; traceIndex++) {
            const trace = graph.data[traceIndex];
            if (!trace.customdata || String(trace.mode || '').indexOf('markers') < 0) {
                continue;
            }
            traceIndices.push(traceIndex);
            if (ids === null) {
                selections.push(null);
                continue;
            }
            const selectedPoints = [];
            for (let pointIndex = 0; pointIndex < trace.customdata.length; pointIndex++) {
                const custom = trace.customdata[pointIndex];
                const cellId = Array.isArray(custom) ? custom[0] : custom;
                if (cellId != null && ids.has(String(cellId))) {
                    selectedPoints.push(pointIndex);
                }
            }
            selections.push(selectedPoints);
        }
        if (traceIndices.length) {
            window.Plotly.restyle(graph, {selectedpoints: selections}, traceIndices);
        }
    }

    highlight('__LEFT_ID__');
    highlight('__RIGHT_ID__');
    return noUpdate;
}
"""


def _feature_options(source, query, current):
    matches = source.search_all_features(query, limit=20)
    current_values = _feature_list(current)
    values = list(dict.fromkeys(current_values + matches))
    return [{"label": value, "value": value} for value in values]


def _feature_list(values):
    if isinstance(values, str):
        values = [values]
    return list(dict.fromkeys(str(value) for value in (values or []) if value))


def _committed_feature_sets(comparison):
    feature_a = _feature_list((comparison or {}).get("feature_a"))
    feature_b = _feature_list((comparison or {}).get("feature_b"))
    if not feature_a or not feature_b:
        raise ValueError("Select at least one feature in both sets.")
    return feature_a, feature_b


def _feature_set_label(features, role):
    if len(features) == 1:
        return features[0]
    return f"Feature set {role} score (n={len(features)})"


def _selected_positions(source, filtered_data, selected_cells):
    raw_indices = (filtered_data or {}).get("cell_indices")
    if raw_indices is None:
        positions = np.arange(len(source.cell_ids), dtype=np.int64)
    else:
        positions = np.asarray(raw_indices, dtype=np.int64)
        positions = positions[(positions >= 0) & (positions < len(source.cell_ids))]

    if selected_cells:
        selected_positions = source.cell_ids.get_indexer(
            np.asarray(selected_cells, dtype=str)
        )
        selected_positions = selected_positions[selected_positions >= 0]
        selected_mask = np.zeros(len(source.cell_ids), dtype=bool)
        selected_mask[selected_positions] = True
        positions = positions[selected_mask[positions]]
    return positions


def _pair_data(source, feature_a, feature_b, filtered_data, selected_cells):
    feature_a = _feature_list(feature_a)
    feature_b = _feature_list(feature_b)
    if not feature_a or not feature_b:
        raise ValueError("Select at least one feature in both sets.")
    positions = _selected_positions(source, filtered_data, selected_cells)
    if len(positions) < 3:
        raise ValueError("At least three cells must remain after filtering.")

    score_a = source.score_features(feature_a)
    score_b = source.score_features(feature_b)
    analysis = analyze_concordance(
        np.asarray(score_a)[positions],
        np.asarray(score_b)[positions],
    )
    return source.base_adata, positions, analysis


def register_cross_modal_concordance_callbacks(
    app,
    source,
    prefix,
    *,
    make_cache_key=None,
    hash_list_signature=None,
    cached_figure_get=None,
    cached_figure_set=None,
):
    """Register the paired-feature concordance callbacks."""

    pair_cache = OrderedDict()
    group_score_cache = OrderedDict()
    memo_lock = Lock()

    def figure_cache_key(
        kind, feature_a, feature_b, filtered_data, cells_hash, **kwargs
    ):
        if make_cache_key is None:
            return None

        filtered_cells_hash = (filtered_data or {}).get("cell_indices_hash")
        if filtered_cells_hash is None and hash_list_signature is not None:
            filtered_cells_hash = hash_list_signature(
                (filtered_data or {}).get("cell_indices")
            )

        return make_cache_key(
            kind,
            source.base_adata,
            prefix=prefix,
            feature_a=feature_a,
            feature_b=feature_b,
            filtered_cells=filtered_cells_hash,
            selected_cells=cells_hash,
            **kwargs,
        )

    def get_cached_figure(cache_key):
        if cache_key is None or cached_figure_get is None:
            return None
        return cached_figure_get(cache_key)

    def cache_figure(cache_key, figure):
        if cache_key is not None and cached_figure_set is not None:
            cached_figure_set(cache_key, figure)
        return figure

    def memoized_pair_data(
        feature_a,
        feature_b,
        filtered_data,
        cells_hash,
        selected_cells,
    ):
        cache_key = figure_cache_key(
            "cross-modal-pair-data",
            feature_a,
            feature_b,
            filtered_data,
            cells_hash,
        )
        if cache_key is not None:
            with memo_lock:
                cached = pair_cache.get(cache_key)
                if cached is not None:
                    pair_cache.move_to_end(cache_key)
                    return cached

        result = _pair_data(
            source, feature_a, feature_b, filtered_data, selected_cells
        )
        if cache_key is not None:
            with memo_lock:
                pair_cache[cache_key] = result
                pair_cache.move_to_end(cache_key)
                while len(pair_cache) > 8:
                    pair_cache.popitem(last=False)
        return result

    def memoized_group_correlations(
        feature_a,
        feature_b,
        group_by,
        filtered_data,
        cells_hash,
        frame,
        positions,
        analysis,
    ):
        cache_key = figure_cache_key(
            "cross-modal-group-correlation",
            feature_a,
            feature_b,
            filtered_data,
            cells_hash,
            group_by=group_by,
        )
        if cache_key is None:
            return calculate_group_correlations(
                analysis,
                frame.obs[group_by].iloc[positions].to_numpy(),
            )

        with memo_lock:
            cached = group_score_cache.get(cache_key)
            if cached is not None:
                group_score_cache.move_to_end(cache_key)
                return cached
            scores = calculate_group_correlations(
                analysis,
                frame.obs[group_by].iloc[positions].to_numpy(),
            )
            group_score_cache[cache_key] = scores
            group_score_cache.move_to_end(cache_key)
            while len(group_score_cache) > 8:
                group_score_cache.popitem(last=False)
            return scores

    def register_feature_search(role):
        @app.callback(
            Output(f"{prefix}-concordance-feature-{role}", "options"),
            Input(f"{prefix}-concordance-feature-{role}", "search_value"),
            State(f"{prefix}-concordance-feature-{role}", "value"),
        )
        def update_feature_options(query, current):
            return _feature_options(source, query, current)

    register_feature_search("a")
    register_feature_search("b")

    app.clientside_callback(
        """
        function(featureA, featureB) {
            return !(Array.isArray(featureA) && featureA.length > 0 &&
                     Array.isArray(featureB) && featureB.length > 0);
        }
        """,
        Output(f"{prefix}-concordance-update-comparison", "disabled"),
        Input(f"{prefix}-concordance-feature-a", "value"),
        Input(f"{prefix}-concordance-feature-b", "value"),
    )

    @app.callback(
        Output(f"{prefix}-concordance-committed-features", "data"),
        Output(f"{prefix}-concordance-update-status", "children"),
        Input(f"{prefix}-concordance-update-comparison", "n_clicks"),
        State(f"{prefix}-concordance-feature-a", "value"),
        State(f"{prefix}-concordance-feature-b", "value"),
        prevent_initial_call=True,
    )
    def commit_feature_sets(_n_clicks, feature_a, feature_b):
        feature_a = _feature_list(feature_a)
        feature_b = _feature_list(feature_b)
        try:
            if not feature_a or not feature_b:
                raise ValueError("Select at least one feature in both sets.")
            invalid = [
                feature
                for feature in (*feature_a, *feature_b)
                if not source.is_feature(feature)
            ]
            if invalid:
                raise ValueError(f"Feature not found: {invalid[0]}")
            # Validate and warm the score cache before replacing the active plots.
            score_a = source.score_features(feature_a)
            score_b = source.score_features(feature_b)
            analyze_concordance(score_a, score_b)
        except ValueError as exc:
            return no_update, f"Could not update: {exc}"

        return (
            {"feature_a": feature_a, "feature_b": feature_b},
            "",
        )

    app.clientside_callback(
        """
        function(viewMode, groupBy) {
            const callbackContext = window.dash_clientside.callback_context;
            const triggered = (callbackContext.triggered || [])[0];
            const prop = triggered ? (triggered.prop_id || '') : '';
            if (viewMode !== 'group_correlation' &&
                prop.endsWith('concordance-group-by.value')) {
                return window.dash_clientside.no_update;
            }
            return {
                view_mode: viewMode,
                group_by: viewMode === 'group_correlation' ? groupBy : null
            };
        }
        """,
        Output(f"{prefix}-concordance-main-view-trigger", "data"),
        Input(f"{prefix}-concordance-view-mode", "value"),
        Input(f"{prefix}-concordance-group-by", "value"),
    )

    @app.callback(
        Output(f"{prefix}-concordance-embedding-plot", "figure"),
        Input(f"{prefix}-concordance-committed-features", "data"),
        Input(f"{prefix}-concordance-main-view-trigger", "data"),
        Input(f"{prefix}-concordance-embedding", "value"),
        Input(f"{prefix}-global-filtered-data", "data"),
        Input(f"{prefix}-selected-cells-hash", "data"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
        State(f"{prefix}-concordance-view-mode", "value"),
        State(f"{prefix}-concordance-group-by", "value"),
        State(f"{prefix}-selected-cells-store", "data"),
    )
    def update_concordance_embedding(
        comparison,
        _main_view_trigger,
        embedding,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        view_mode,
        group_by,
        selected_cells,
    ):
        if active_workspace != EXPLORATION_WORKSPACE or active_tab != _ACTIVE_TAB:
            return no_update

        try:
            feature_a, feature_b = _committed_feature_sets(comparison)
        except ValueError as exc:
            return empty_concordance_figure(str(exc))

        cache_key = figure_cache_key(
            "cross-modal-concordance-embedding",
            feature_a,
            feature_b,
            filtered_data,
            cells_hash,
            view_mode=view_mode,
            embedding=embedding,
            group_by=group_by if view_mode == GROUP_CORRELATION else None,
        )
        cached_figure = get_cached_figure(cache_key)
        if cached_figure is not None:
            return cached_figure

        try:
            if view_mode not in _VIEW_MODES:
                raise ValueError("Select a valid concordance view.")
            if embedding not in source.embedding_names:
                raise ValueError("Select a valid embedding.")
            if view_mode == GROUP_CORRELATION and group_by not in source.obs_names:
                raise ValueError("Select categorical metadata for Correlation view.")
            frame, positions, analysis = memoized_pair_data(
                feature_a,
                feature_b,
                filtered_data,
                cells_hash,
                selected_cells,
            )
            group_scores = (
                memoized_group_correlations(
                    feature_a,
                    feature_b,
                    group_by,
                    filtered_data,
                    cells_hash,
                    frame,
                    positions,
                    analysis,
                )
                if view_mode == GROUP_CORRELATION
                else None
            )
            figure = build_concordance_embedding(
                np.asarray(source.base_adata.obsm[embedding])[positions],
                analysis,
                frame.obs_names.to_numpy()[positions],
                _feature_set_label(feature_a, "A"),
                _feature_set_label(feature_b, "B"),
                embedding,
                view_mode=view_mode,
                group_correlation=group_scores,
                group_by=group_by,
            )
        except ValueError as exc:
            return empty_concordance_figure(str(exc))

        cache_figure(cache_key, figure)
        return figure

    @app.callback(
        Output(f"{prefix}-concordance-feature-plot", "figure"),
        Input(f"{prefix}-concordance-committed-features", "data"),
        Input(f"{prefix}-concordance-main-view-trigger", "data"),
        Input(f"{prefix}-global-filtered-data", "data"),
        Input(f"{prefix}-selected-cells-hash", "data"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
        State(f"{prefix}-concordance-view-mode", "value"),
        State(f"{prefix}-concordance-group-by", "value"),
        State(f"{prefix}-selected-cells-store", "data"),
    )
    def update_feature_relationship(
        comparison,
        _main_view_trigger,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        view_mode,
        group_by,
        selected_cells,
    ):
        if active_workspace != EXPLORATION_WORKSPACE or active_tab != _ACTIVE_TAB:
            return no_update
        try:
            feature_a, feature_b = _committed_feature_sets(comparison)
        except ValueError as exc:
            return empty_concordance_figure(str(exc))
        cache_key = figure_cache_key(
            "cross-modal-feature-relationship",
            feature_a,
            feature_b,
            filtered_data,
            cells_hash,
            view_mode=view_mode,
            group_by=group_by if view_mode == GROUP_CORRELATION else None,
        )
        cached_figure = get_cached_figure(cache_key)
        if cached_figure is not None:
            return cached_figure

        try:
            if view_mode not in _VIEW_MODES:
                raise ValueError("Select a valid concordance view.")
            if view_mode == GROUP_CORRELATION and group_by not in source.obs_names:
                raise ValueError("Select categorical metadata for Correlation view.")
            frame, positions, analysis = memoized_pair_data(
                feature_a,
                feature_b,
                filtered_data,
                cells_hash,
                selected_cells,
            )
            group_scores = (
                memoized_group_correlations(
                    feature_a,
                    feature_b,
                    group_by,
                    filtered_data,
                    cells_hash,
                    frame,
                    positions,
                    analysis,
                )
                if view_mode == GROUP_CORRELATION
                else None
            )
            figure = build_feature_scatter(
                analysis,
                frame.obs_names.to_numpy()[positions],
                _feature_set_label(feature_a, "A"),
                _feature_set_label(feature_b, "B"),
                view_mode=view_mode,
                group_correlation=group_scores,
                group_by=group_by,
            )
        except ValueError as exc:
            return empty_concordance_figure(str(exc))

        cache_figure(cache_key, figure)
        return figure

    @app.callback(
        Output(f"{prefix}-concordance-group-summary-plot", "figure"),
        Input(f"{prefix}-concordance-committed-features", "data"),
        Input(f"{prefix}-concordance-view-mode", "value"),
        Input(f"{prefix}-concordance-group-by", "value"),
        Input(f"{prefix}-global-filtered-data", "data"),
        Input(f"{prefix}-selected-cells-hash", "data"),
        Input(f"{prefix}-visualization-workspace-tabs", "value"),
        Input(f"{prefix}-exploratory-tabs", "value"),
        State(f"{prefix}-selected-cells-store", "data"),
    )
    def update_group_summary(
        comparison,
        view_mode,
        group_by,
        filtered_data,
        cells_hash,
        active_workspace,
        active_tab,
        selected_cells,
    ):
        if active_workspace != EXPLORATION_WORKSPACE or active_tab != _ACTIVE_TAB:
            return no_update
        try:
            feature_a, feature_b = _committed_feature_sets(comparison)
        except ValueError as exc:
            return empty_concordance_figure(str(exc))
        if view_mode not in _VIEW_MODES:
            return empty_concordance_figure("Select a valid concordance view.")
        if group_by not in source.obs_names:
            return empty_concordance_figure("Select categorical metadata.")

        cache_key = figure_cache_key(
            "cross-modal-group-summary",
            feature_a,
            feature_b,
            filtered_data,
            cells_hash,
            view_mode=view_mode,
            group_by=group_by,
        )
        cached_figure = get_cached_figure(cache_key)
        if cached_figure is not None:
            return cached_figure

        try:
            frame, positions, analysis = memoized_pair_data(
                feature_a,
                feature_b,
                filtered_data,
                cells_hash,
                selected_cells,
            )
            raw_groups = frame.obs[group_by].iloc[positions].to_numpy()
            group_scores = (
                memoized_group_correlations(
                    feature_a,
                    feature_b,
                    group_by,
                    filtered_data,
                    cells_hash,
                    frame,
                    positions,
                    analysis,
                )
                if view_mode == GROUP_CORRELATION
                else None
            )
            figure = build_group_summary_heatmap(
                analysis,
                group_by,
                view_mode=view_mode,
                groups=raw_groups,
                group_correlation=group_scores,
            )
        except ValueError as exc:
            return empty_concordance_figure(str(exc))

        cache_figure(cache_key, figure)
        return figure

    app.clientside_callback(
        _LINKED_SELECTION_JS,
        Output(f"{prefix}-concordance-linked-cells", "data"),
        Input(f"{prefix}-concordance-embedding-plot", "selectedData"),
        Input(f"{prefix}-concordance-feature-plot", "selectedData"),
        Input(f"{prefix}-concordance-committed-features", "data"),
        Input(f"{prefix}-concordance-main-view-trigger", "data"),
        Input(f"{prefix}-global-filtered-data", "data"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        _LINKED_HIGHLIGHT_JS
        .replace("__LEFT_ID__", f"{prefix}-concordance-embedding-plot")
        .replace("__RIGHT_ID__", f"{prefix}-concordance-feature-plot"),
        Output(f"{prefix}-concordance-highlight-link", "data"),
        Input(f"{prefix}-concordance-linked-cells", "data"),
        Input(f"{prefix}-concordance-embedding-plot", "figure"),
        Input(f"{prefix}-concordance-feature-plot", "figure"),
        prevent_initial_call=True,
    )
