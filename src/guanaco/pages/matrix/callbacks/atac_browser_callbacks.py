from dash import Input, Output, State, callback_context, no_update
from dash.exceptions import PreventUpdate

from guanaco.pages.matrix.plots.atac_browser import (
    build_peak_index,
    coerce_region,
    compute_atac_signal,
    format_locus,
    parse_locus,
    plot_atac_browser,
)
from guanaco.pages.matrix.plots.gene_annotation import (
    find_gene_region,
    load_gene_annotation,
    query_gene_models,
)
from guanaco.utils.colors import resolve_discrete_palette
from guanaco.utils.obs_utils import sorted_categories


def _gene_track_label(gene_annotation_path) -> str:
    """Badge text for the gene track: the reference filename without its extension."""
    if not gene_annotation_path:
        return "Genes"
    name = str(gene_annotation_path).split("/")[-1]
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for ext in (".gtf", ".gff3", ".gff"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name or "Genes"


def _range_from_relayout(relayout):
    if not relayout:
        return None
    if "xaxis.range[0]" in relayout and "xaxis.range[1]" in relayout:
        return int(float(relayout["xaxis.range[0]"])), int(float(relayout["xaxis.range[1]"]))
    if "xaxis.range" in relayout and isinstance(relayout["xaxis.range"], (list, tuple)):
        values = relayout["xaxis.range"]
        if len(values) >= 2:
            return int(float(values[0])), int(float(values[1]))
    return None


def register_atac_browser_callbacks(app, adata, prefix, gene_annotation_path=None, color_config=None):
    index = build_peak_index(adata)
    annotation_state = {"loaded": False, "index": None, "error": None}

    def _annotation_index():
        if not gene_annotation_path:
            return None
        if not annotation_state["loaded"]:
            try:
                annotation_state["index"] = load_gene_annotation(gene_annotation_path)
            except Exception as exc:
                annotation_state["error"] = str(exc)
                annotation_state["index"] = None
            annotation_state["loaded"] = True
        return annotation_state["index"]

    @app.callback(
        Output(f"{prefix}-atac-browser-region-store", "data"),
        Output(f"{prefix}-atac-browser-message", "children"),
        Input(f"{prefix}-atac-browser-go", "n_clicks"),
        Input(f"{prefix}-atac-browser-graph", "relayoutData"),
        State(f"{prefix}-atac-browser-locus", "value"),
        State(f"{prefix}-atac-browser-region-store", "data"),
        prevent_initial_call=True,
    )
    def update_atac_region(
        go_clicks,
        relayout,
        locus_value,
        current_region,
    ):
        triggered = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        region = current_region or {"chrom": index.chroms[0], "start": 0, "end": 1_000_000}

        if triggered.endswith("-atac-browser-go.n_clicks"):
            parsed = parse_locus(locus_value)
            if parsed is None:
                annotation = _annotation_index()
                if annotation_state["error"]:
                    return no_update, f"Gene model unavailable: {annotation_state['error']}"
                gene_region = find_gene_region(annotation, locus_value)
                if gene_region is None:
                    return no_update, "Enter a locus like chr1:100,000-200,000 or a gene name present in the GTF/GFF3."
                region = gene_region
            else:
                chrom, start, end = parsed
                if chrom not in index.starts:
                    return no_update, f"{chrom} is not present in the ATAC peak index."
                region = {"chrom": chrom, "start": start, "end": end}
        elif triggered.endswith("-atac-browser-graph.relayoutData"):
            x_range = _range_from_relayout(relayout)
            if x_range is None:
                raise PreventUpdate
            region = {"chrom": region["chrom"], "start": x_range[0], "end": x_range[1]}
        else:
            raise PreventUpdate

        region = coerce_region(index, region)
        # The figure's uirevision is this nav token. On an explicit navigation
        # (Update plot / gene search) we mint a new token so plotly jumps to the
        # target; on a zoom/pan re-render we keep the previous token so plotly
        # preserves the view the user is currently interacting with (no jitter).
        if triggered.endswith("-atac-browser-go.n_clicks"):
            region["nav"] = format_locus(str(region["chrom"]), int(region["start"]), int(region["end"]))
        else:
            region["nav"] = (current_region or {}).get("nav")
        return region, format_locus(str(region["chrom"]), int(region["start"]), int(region["end"]))

    @app.callback(
        Output(f"{prefix}-atac-browser-graph", "figure"),
        Output(f"{prefix}-atac-browser-message", "children", allow_duplicate=True),
        Input(f"{prefix}-atac-browser-region-store", "data"),
        Input(f"{prefix}-single-cell-annotation-dropdown", "value"),
        Input(f"{prefix}-single-cell-label-selection", "value"),
        Input(f"{prefix}-atac-browser-metric", "value"),
        Input(f"{prefix}-atac-browser-yscale", "value"),
        Input(f"{prefix}-discrete-color-map-dropdown", "value"),
        Input(f"{prefix}-selected-cells-store", "data"),
        Input(f"{prefix}-single-cell-tabs", "value"),
        prevent_initial_call="initial_duplicate",
    )
    def render_atac_browser(region, annotation, labels, metric, y_mode, discrete_color_map, selected_cells, active_tab):
        # Only build when the ATAC tab is showing -- switching to it fires this
        # callback (the tab value is an Input), so the figure is still up to date,
        # but a scatter selection on another tab no longer recomputes the signal.
        if active_tab != "peak-browser-tab":
            raise PreventUpdate
        if not region:
            raise PreventUpdate
        metric = metric if metric in {"mean", "detection"} else "mean"
        y_mode = y_mode if y_mode in {"auto", "shared"} else "auto"

        # Track grouping/order/colour all come from the left panel, so the ATAC
        # tracks line up with the heatmap/violin/dotplot for the same dataset.
        group_order = None
        color_map = None
        if annotation and annotation in adata.obs.columns:
            group_order = sorted_categories(adata, annotation)
            palette = resolve_discrete_palette(discrete_color_map, len(group_order), default=color_config)
            if palette:
                color_map = {str(label): palette[i % len(palette)] for i, label in enumerate(group_order)}

        payload = compute_atac_signal(
            adata,
            region,
            selected_cells=selected_cells,
            groupby=annotation,
            labels=labels,
            group_order=group_order,
            metric=metric,
        )
        gene_annotation = _annotation_index()
        gene_models = None
        gene_track_label = _gene_track_label(gene_annotation_path)
        # Only surface a message when the gene model failed to load; on success the
        # status line stays empty (the locus box and track badges say enough).
        message = f"Gene model unavailable: {annotation_state['error']}" if annotation_state["error"] else ""
        if not annotation_state["error"] and gene_annotation is not None:
            normalized_region = payload["region"]
            gene_models = query_gene_models(
                gene_annotation,
                str(normalized_region["chrom"]),
                int(normalized_region["start"]),
                int(normalized_region["end"]),
            )
        return (
            plot_atac_browser(
                payload,
                gene_models=gene_models,
                color_map=color_map,
                gene_track_label=gene_track_label,
                y_mode=y_mode,
                uirevision=region.get("nav") if isinstance(region, dict) else None,
            ),
            message,
        )
