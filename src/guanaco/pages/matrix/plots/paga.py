import math

import numpy as np
import pandas as pd
from dash import html
from plotly.colors import sample_colorscale
from scipy import sparse

try:
    import dash_cytoscape as cyto
except Exception:
    cyto = None

from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.gene_extraction_utils import extract_gene_expression


MAX_RENDERED_CYTOSCAPE_PIE_SLICES = 16
CYTOSCAPE_EMPTY_COLOR = "#B0B7C3"


def _infer_paga_groupby(adata, paga_uns, n_nodes):
    group_key = paga_uns.get("groups")
    if isinstance(group_key, bytes):
        group_key = group_key.decode("utf-8")
    if isinstance(group_key, str) and group_key in adata.obs.columns:
        return group_key

    candidate_columns = []
    for col in adata.obs.columns:
        series = adata.obs[col]
        if hasattr(series, "cat"):
            try:
                n_categories = len(series.cat.categories)
            except Exception:
                continue
            if n_categories == n_nodes:
                candidate_columns.append(col)

    if len(candidate_columns) == 1:
        return candidate_columns[0]

    preferred = ["leiden", "louvain", "clusters", "cluster", "cell_type", "annotation"]
    for col in preferred:
        if col in candidate_columns:
            return col

    raise ValueError(
        "Could not infer the obs column used for PAGA. "
        "Expected adata.uns['paga']['groups'] or a categorical obs column matching the PAGA node count."
    )


def _resolve_paga_positions(connectivities, pos):
    if pos is not None:
        arr = np.asarray(pos, dtype=float)
        if arr.ndim == 2 and arr.shape[1] == 2:
            return arr

    n_nodes = connectivities.shape[0]
    try:
        import networkx as nx

        graph = nx.from_scipy_sparse_array(connectivities)
        layout = nx.spring_layout(graph, seed=0, weight="weight")
        return np.asarray([layout[i] for i in range(n_nodes)], dtype=float)
    except Exception:
        angles = np.linspace(0, 2 * math.pi, n_nodes, endpoint=False)
        return np.column_stack([np.cos(angles), np.sin(angles)])


def _continuous_group_summary(values, group_series, group_categories):
    df = pd.DataFrame({"group": group_series.to_numpy(), "value": np.asarray(values, dtype=float)})
    grouped = df.groupby("group", observed=True)["value"].mean()
    return grouped.reindex(group_categories).to_numpy(dtype=float)


def _resolve_obs_color_mapping(adata, obs_key, categories, fallback_palette, *, prefer_fallback=False):
    if prefer_fallback:
        return {str(cat): fallback_palette[i % len(fallback_palette)] for i, cat in enumerate(categories)}

    colors_key = f"{obs_key}_colors"
    if colors_key in adata.uns and len(adata.uns[colors_key]) >= len(categories):
        return {str(cat): adata.uns[colors_key][i] for i, cat in enumerate(categories)}
    return {str(cat): fallback_palette[i % len(fallback_palette)] for i, cat in enumerate(categories)}


def _edge_width(weight, min_weight, max_weight, min_width=1.0, max_width=10.0):
    if max_weight <= min_weight:
        return (min_width + max_width) / 2.0
    return min_width + (weight - min_weight) / (max_weight - min_weight) * (max_width - min_width)


def _ordered_obs_categories(obs_values):
    return sorted(str(cat) for cat in pd.unique(obs_values.dropna()) if not pd.isna(cat))


def _prepare_paga_context(
    adata,
    *,
    selected_annotation=None,
    selected_labels=None,
    selected_cells=None,
):
    if "paga" not in adata.uns or "connectivities" not in adata.uns["paga"]:
        raise ValueError("PAGA graph not found. Expected adata.uns['paga']['connectivities'].")

    paga_uns = adata.uns["paga"]
    connectivities = paga_uns["connectivities"]
    if not sparse.issparse(connectivities):
        connectivities = sparse.csr_matrix(connectivities)
    connectivities = connectivities.tocoo()

    n_nodes = connectivities.shape[0]
    paga_groupby = _infer_paga_groupby(adata, paga_uns, n_nodes)
    group_series = adata.obs[paga_groupby]
    if not hasattr(group_series, "cat"):
        group_series = group_series.astype("category")
    group_categories = list(group_series.cat.categories)

    if len(group_categories) != n_nodes:
        raise ValueError(
            f"PAGA node count ({n_nodes}) does not match categories in '{paga_groupby}' ({len(group_categories)})."
        )

    positions = _resolve_paga_positions(connectivities.tocsr(), paga_uns.get("pos"))

    selection_mask = np.ones(adata.n_obs, dtype=bool)
    if selected_labels and selected_annotation in adata.obs.columns:
        selection_mask &= adata.obs[selected_annotation].isin(selected_labels).to_numpy()
    if selected_cells:
        selected_cells_set = set(selected_cells)
        selection_mask &= adata.obs_names.isin(selected_cells_set)

    filtered_obs = adata.obs.loc[selection_mask].copy()
    filtered_group_series = filtered_obs[paga_groupby]
    filtered_counts = (
        filtered_group_series.value_counts()
        .reindex(group_categories, fill_value=0)
        .to_numpy(dtype=int)
    )
    node_sizes = (
        14 + 28 * np.sqrt(filtered_counts / filtered_counts.max())
        if filtered_counts.max() > 0
        else np.full(n_nodes, 14.0)
    )
    hover_lines = [
        f"group={group}<br>cells={count}"
        for group, count in zip(group_categories, filtered_counts)
    ]

    visible_weights = [
        float(w)
        for i, j, w in zip(connectivities.row, connectivities.col, connectivities.data)
        if i < j and float(w) >= 0
    ]

    return {
        "paga_uns": paga_uns,
        "connectivities": connectivities,
        "n_nodes": n_nodes,
        "paga_groupby": paga_groupby,
        "group_series": group_series,
        "group_categories": group_categories,
        "positions": positions,
        "selection_mask": selection_mask,
        "filtered_obs": filtered_obs,
        "filtered_group_series": filtered_group_series,
        "filtered_counts": filtered_counts,
        "node_sizes": node_sizes,
        "hover_lines": hover_lines,
        "all_weights": visible_weights,
    }


def _categorical_group_composition(values, group_series, group_categories, category_order=None):
    if category_order is None:
        category_order = _ordered_obs_categories(values)
    series = pd.Series(values, index=group_series.index, dtype="object")
    non_null_values = series.dropna().astype(str)
    total_counts = non_null_values.value_counts()

    ordered_present = [cat for cat in category_order if cat in total_counts.index]
    ordered_present.extend(cat for cat in total_counts.index if cat not in ordered_present)

    effective_categories = list(ordered_present)

    proportions = {}
    detailed_counts = {}
    for group in group_categories:
        group_values = series[group_series == group].dropna()
        if group_values.empty:
            proportions[group] = {}
            detailed_counts[group] = {}
            continue

        counts = group_values.astype(str).value_counts()
        detailed_counts[group] = counts.to_dict()
        total = float(counts.sum())
        group_props = {cat: counts.get(cat, 0) / total * 100.0 for cat in effective_categories if counts.get(cat, 0) > 0}

        proportions[group] = group_props

    return effective_categories, proportions, detailed_counts


def _format_pie_hover_text(group, cell_count, group_props, detail):
    lines = [f"group={group}", f"cells={cell_count}"]
    if not group_props:
        lines.append("no composition data")
        return "<br>".join(lines)

    total = sum(detail.values()) if detail else 0
    for label, percent in group_props.items():
        count = detail.get(label)
        if count is None and total:
            count = int(round(percent / 100.0 * total))
        if count is None:
            lines.append(f"{label}: {percent:.1f}%")
        else:
            lines.append(f"{label}: {count} ({percent:.1f}%)")
    return "<br>".join(lines)


def _format_continuous_hover_text(group, cell_count, value, label):
    if not np.isfinite(value):
        value_text = "NA"
    else:
        value_text = f"{float(value):.4g}"
    return "<br>".join([f"group={group}", f"cells={cell_count}", f"mean {label}={value_text}"])


def _values_to_colors(values, continuous_color_map, na_color=CYTOSCAPE_EMPTY_COLOR):
    values = np.asarray(values, dtype=float)
    finite_mask = np.isfinite(values)
    if not finite_mask.any():
        return [na_color] * len(values)

    finite_values = values[finite_mask]
    vmin = float(np.min(finite_values))
    vmax = float(np.max(finite_values))
    colorscale = resolve_continuous_colorscale(continuous_color_map)
    colors = []
    for value in values:
        if not np.isfinite(value):
            colors.append(na_color)
            continue
        if vmax <= vmin:
            normalized = 0.5
        else:
            normalized = (float(value) - vmin) / (vmax - vmin)
        colors.append(sample_colorscale(colorscale, [float(np.clip(normalized, 0.0, 1.0))])[0])
    return colors


def _continuous_colorbar_legend(values, continuous_color_map, title, na_color=CYTOSCAPE_EMPTY_COLOR):
    values = np.asarray(values, dtype=float)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return html.Div(
            [
                html.Div(title, style={"fontWeight": "bold", "fontSize": "12px", "marginBottom": "6px"}),
                html.Div("No finite values", style={"fontSize": "12px", "color": "#6B7280"}),
            ],
            style={"paddingBottom": "8px", "marginBottom": "8px", "borderBottom": "1px solid #D8DDE6"},
        )

    vmin = float(finite_values.min())
    vmax = float(finite_values.max())
    colorscale = resolve_continuous_colorscale(continuous_color_map)
    sampled = sample_colorscale(colorscale, np.linspace(0, 1, 9))
    gradient = ", ".join(f"{color} {idx * 12.5:.1f}%" for idx, color in enumerate(sampled))
    return html.Div(
        [
            html.Div(title, style={"fontWeight": "bold", "fontSize": "12px", "marginBottom": "6px"}),
            html.Div(
                style={
                    "height": "14px",
                    "width": "100%",
                    "background": f"linear-gradient(to right, {gradient})",
                    "border": "1px solid #AEB6C2",
                    "borderRadius": "2px",
                    "marginBottom": "4px",
                }
            ),
            html.Div(
                [
                    html.Span(f"{vmin:.4g}"),
                    html.Span(f"{vmax:.4g}"),
                ],
                style={"display": "flex", "justifyContent": "space-between", "fontSize": "11px", "color": "#2F3E46"},
            ),
        ],
        style={"paddingBottom": "8px", "marginBottom": "8px", "borderBottom": "1px solid #D8DDE6"},
    )


def _scale_positions_for_cytoscape(positions, width=900.0, height=720.0, padding=70.0):
    positions = np.asarray(positions, dtype=float)
    if positions.size == 0:
        return positions

    x_values = positions[:, 0]
    y_values = positions[:, 1]
    x_span = float(x_values.max() - x_values.min()) if len(x_values) > 1 else 1.0
    y_span = float(y_values.max() - y_values.min()) if len(y_values) > 1 else 1.0
    span = max(x_span, y_span, 1e-9)
    scale = min((width - 2 * padding) / span, (height - 2 * padding) / span)

    x_center = float((x_values.max() + x_values.min()) / 2.0)
    y_center = float((y_values.max() + y_values.min()) / 2.0)
    scaled_x = (x_values - x_center) * scale + width / 2.0
    scaled_y = height / 2.0 - (y_values - y_center) * scale
    return np.column_stack([scaled_x, scaled_y])


def _build_paga_legend_items(color_map):
    if not color_map:
        return None

    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        style={
                            "display": "inline-block",
                            "width": "12px",
                            "height": "12px",
                            "backgroundColor": color,
                            "border": "1px solid #AEB6C2",
                            "borderRadius": "2px",
                            "marginRight": "6px",
                            "verticalAlign": "middle",
                        }
                    ),
                    html.Span(str(label), style={"fontSize": "12px", "color": "#2F3E46"}),
                ],
                style={"display": "flex", "alignItems": "center", "marginBottom": "6px"},
            )
            for label, color in color_map.items()
        ],
        style={"paddingBottom": "8px", "marginBottom": "8px", "borderBottom": "1px solid #D8DDE6"},
    )


def build_paga_cytoscape(
    adata,
    *,
    component_id="paga-cytoscape-view",
    color_mode="obs",
    obs_key=None,
    gene=None,
    continuous_color_map="Viridis",
    discrete_palette=None,
    edge_threshold=0.03,
    selected_annotation=None,
    selected_labels=None,
    selected_cells=None,
):
    if cyto is None:
        return html.Div(
            "dash-cytoscape is not installed. Add `dash-cytoscape` to the environment to render pie-chart PAGA nodes.",
            style={"padding": "16px", "color": "#7A1F1F", "backgroundColor": "#FDECEC", "border": "1px solid #F5C2C7", "borderRadius": "6px"},
        )

    context = _prepare_paga_context(
        adata,
        selected_annotation=selected_annotation,
        selected_labels=selected_labels,
        selected_cells=selected_cells,
    )
    connectivities = context["connectivities"]
    group_categories = context["group_categories"]
    positions = _scale_positions_for_cytoscape(context["positions"])
    filtered_obs = context["filtered_obs"]
    filtered_group_series = context["filtered_group_series"]
    node_sizes = context["node_sizes"]
    filtered_counts = context["filtered_counts"]
    paga_groupby = context["paga_groupby"]

    visible_weights = [
        float(weight)
        for i, j, weight in zip(connectivities.row, connectivities.col, connectivities.data)
        if i < j and float(weight) >= edge_threshold
    ]
    min_weight = min(visible_weights) if visible_weights else float(edge_threshold)
    max_weight = max(visible_weights) if visible_weights else max(float(edge_threshold), 1.0)

    node_ids = {index: f"paga-node-{index}" for index in range(len(group_categories))}
    node_elements = []
    edge_elements = []
    legend_items = None
    continuous_value_label = None

    if color_mode == "gene":
        if not gene:
            raise ValueError("Select a gene to color the PAGA graph.")
        gene_expr = extract_gene_expression(adata, gene, use_cache=True, dtype=np.float32)
        gene_expr = gene_expr[context["selection_mask"]]
        node_color_values = _continuous_group_summary(gene_expr, filtered_group_series, group_categories)
        node_fill_colors = _values_to_colors(node_color_values, continuous_color_map)
        legend_items = _continuous_colorbar_legend(node_color_values, continuous_color_map, f"mean {gene}")
        continuous_value_label = gene
        node_mode = "solid"
    else:
        if not obs_key:
            raise ValueError("Select an obs column to color the PAGA graph.")

        obs_values = filtered_obs[obs_key]
        if pd.api.types.is_numeric_dtype(obs_values):
            node_color_values = _continuous_group_summary(obs_values.to_numpy(dtype=float), filtered_group_series, group_categories)
            node_fill_colors = _values_to_colors(node_color_values, continuous_color_map)
            legend_items = _continuous_colorbar_legend(node_color_values, continuous_color_map, f"mean {obs_key}")
            continuous_value_label = obs_key
            node_mode = "solid"
        else:
            obs_categories = _ordered_obs_categories(adata.obs[obs_key])
            pie_categories, pie_proportions, category_details = _categorical_group_composition(
                obs_values,
                filtered_group_series,
                group_categories,
                category_order=obs_categories,
            )
            color_map = _resolve_obs_color_mapping(
                adata,
                obs_key,
                obs_categories,
                discrete_palette or ["#636EFA"],
                prefer_fallback=discrete_palette is not None,
            )
            legend_items = _build_paga_legend_items(color_map)
            node_mode = "pie"

    for index, group in enumerate(group_categories):
        x_pos, y_pos = positions[index]
        node_data = {
            "id": node_ids[index],
            "label": str(group),
            "group": str(group),
            "cells": int(filtered_counts[index]),
            "size_px": float(node_sizes[index]),
            "background_color": CYTOSCAPE_EMPTY_COLOR,
        }
        n_pie_slots = MAX_RENDERED_CYTOSCAPE_PIE_SLICES if node_mode == "pie" else 0
        for slice_index in range(n_pie_slots):
            node_data[f"pie_{slice_index + 1}_size"] = 0.0
            node_data[f"pie_{slice_index + 1}_color"] = CYTOSCAPE_EMPTY_COLOR

        if node_mode == "solid":
            node_data["background_color"] = node_fill_colors[index]
            if continuous_value_label is not None:
                node_data["hover_text"] = _format_continuous_hover_text(
                    group,
                    int(filtered_counts[index]),
                    node_color_values[index],
                    continuous_value_label,
                )
        else:
            node_data["background_color"] = CYTOSCAPE_EMPTY_COLOR
            group_props = pie_proportions.get(group, {})
            visible_categories = [category for category in pie_categories if group_props.get(category, 0.0) > 0]
            if len(visible_categories) > MAX_RENDERED_CYTOSCAPE_PIE_SLICES:
                category_rank = {category: rank for rank, category in enumerate(pie_categories)}
                visible_categories = sorted(
                    visible_categories,
                    key=lambda category: (-group_props.get(category, 0.0), category_rank[category]),
                )
            for slice_index, category in enumerate(visible_categories[:MAX_RENDERED_CYTOSCAPE_PIE_SLICES]):
                node_data[f"pie_{slice_index + 1}_size"] = float(group_props.get(category, 0.0))
                node_data[f"pie_{slice_index + 1}_color"] = color_map.get(category, CYTOSCAPE_EMPTY_COLOR)

            detail = category_details.get(group, {})
            node_data["hover_text"] = _format_pie_hover_text(
                group,
                int(filtered_counts[index]),
                group_props,
                detail,
            )

        node_elements.append(
            {
                "data": node_data,
                "position": {"x": float(x_pos), "y": float(y_pos)},
            }
        )

    for i, j, weight in zip(connectivities.row, connectivities.col, connectivities.data):
        weight = float(weight)
        if i >= j or weight < edge_threshold:
            continue
        edge_elements.append(
            {
                "data": {
                    "id": f"paga-edge-{i}-{j}",
                    "source": node_ids[int(i)],
                    "target": node_ids[int(j)],
                    "weight": weight,
                    "width": float(_edge_width(weight, min_weight, max_weight, min_width=1.0, max_width=10.0)),
                }
            }
        )

    stylesheet = [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "width": "data(size_px)",
                "height": "data(size_px)",
                "background-color": "data(background_color)",
                "border-width": 1,
                "border-color": "white",
                "font-size": 11,
                "color": "#2F3E46",
                "text-valign": "bottom",
                "text-halign": "center",
                "text-margin-y": 10,
                "pie-size": "100%",
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "line-color": "#8A8A8A",
                "opacity": 0.6,
                "width": "data(width)",
            },
        },
    ]

    n_pie_slots = MAX_RENDERED_CYTOSCAPE_PIE_SLICES if node_mode == "pie" else 0
    for slice_index in range(n_pie_slots):
        stylesheet[0]["style"][f"pie-{slice_index + 1}-background-color"] = f"data(pie_{slice_index + 1}_color)"
        stylesheet[0]["style"][f"pie-{slice_index + 1}-background-size"] = f"data(pie_{slice_index + 1}_size)"

    cytoscape_component = cyto.Cytoscape(
        id=component_id,
        layout={"name": "preset", "fit": True, "padding": 40},
        elements=node_elements + edge_elements,
        stylesheet=stylesheet,
        minZoom=0.2,
        maxZoom=3.0,
        autoungrabify=False,
        responsive=True,
        style={"width": "100%", "height": "100%", "flex": "1 1 auto", "backgroundColor": "white"},
    )

    hover_detail = html.Div(
        [
            legend_items if legend_items is not None else None,
            html.Div(id=component_id.replace("-cytoscape-view", "-hover-detail")),
        ],
        style={
            "width": "220px",
            "minWidth": "220px",
            "height": "100%",
            "overflowY": "auto",
            "padding": "8px 10px",
            "fontSize": "12px",
            "color": "#2F3E46",
            "lineHeight": "1.35",
            "borderLeft": "1px solid #D8DDE6",
            "backgroundColor": "#FAFBFC",
            "boxSizing": "border-box",
            "whiteSpace": "normal",
            "overflowWrap": "break-word",
        },
    )
    graph_row = html.Div(
        [cytoscape_component, hover_detail],
        style={
            "minHeight": "0",
            "flex": "1 1 auto",
            "display": "flex",
            "width": "100%",
        },
    )

    return html.Div(
        [graph_row],
        style={"height": "100%", "width": "100%", "display": "flex", "flexDirection": "column"},
    )
