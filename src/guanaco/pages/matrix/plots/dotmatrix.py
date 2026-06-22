import numpy as np
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
import pandas as pd
from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression,
    apply_transformation,
    prewarm_gene_cache,
)
from guanaco.pages.matrix.plots.heatmap import ZSCORE_COLOR_CLIP
from guanaco.data.loader import obs_col

from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as _scipy_dendro
from scipy.spatial.distance import pdist


DOT_SIZE_LEGEND_SCALE = 1.6

# scipy's dendrogram lays leaves out on the icoord axis at evenly spaced ticks:
# leaf i sits at (_DENDRO_LEAF_BASE + i * _DENDRO_LEAF_STEP). _map_to_axis inverts
# this to recover the fractional leaf index from a raw icoord value.
_DENDRO_LEAF_BASE = 5.0
_DENDRO_LEAF_STEP = 10.0


def _colorbar(colorbar_title):
    """Shared colorbar styling for the dot/matrix plots (kept identical across both)."""
    return dict(
        title=colorbar_title,
        tickfont=dict(color='DarkSlateGrey', size=10),
        len=0.6,
        yanchor="middle",
        y=0.5,
        x=0.98,
    )


def _dendro_axis_layout(main_x_right, main_y_top, show_right_dendro, show_top_dendro, clamp_right_domain=False):
    """Build the optional dendrogram axis (x3/y3 right, x4/y4 top) layout entries.

    Returned only for the dendrograms actually shown so empty space is not reserved.
    ``clamp_right_domain`` keeps the right dendrogram inside the paper for the
    matrixplot, which uses a wider main domain than the dotplot.
    """
    extra = {}
    if show_right_dendro:
        right_upper = main_x_right + 0.04
        if clamp_right_domain:
            right_upper = min(right_upper, 0.99)
        extra['xaxis3'] = dict(domain=[main_x_right + 0.01, right_upper], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False)
        extra['yaxis3'] = dict(domain=[0.0, main_y_top], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False)
    if show_top_dendro:
        extra['xaxis4'] = dict(domain=[0.0, main_x_right], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False)
        extra['yaxis4'] = dict(domain=[main_y_top + 0.02, min(main_y_top + 0.14, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False)
    return extra


def _selected_row_indices(adata, group_series, selected_labels, selected_cells):
    use_selected_cells = selected_cells is not None and len(selected_cells) > 0
    if use_selected_cells:
        selected_pos = adata.obs_names.get_indexer(selected_cells)
        row_indices = selected_pos[selected_pos >= 0].astype(np.int64, copy=False)
    elif selected_labels:
        label_mask = group_series.isin(selected_labels).to_numpy()
        row_indices = np.flatnonzero(label_mask).astype(np.int64, copy=False)
    else:
        row_indices = np.arange(adata.n_obs, dtype=np.int64)
    return row_indices, use_selected_cells


def _groups_to_process(group_values, selected_labels, use_selected_cells):
    if selected_labels and not use_selected_cells:
        present_groups = set(pd.unique(group_values))
        return [group for group in selected_labels if group in present_groups]
    return list(pd.unique(group_values))


def _expression_frame(adata, valid_genes, row_indices, layer=None, transformation=None):
    n_cells = row_indices.size
    # Read all genes in one column slice up front so the per-gene loop below hits the
    # cache instead of issuing one slice per gene (slow on large/sparse/backed X).
    prewarm_gene_cache(adata, valid_genes, layer=layer, dtype=np.float32)
    gene_matrix = np.empty((n_cells, len(valid_genes)), dtype=np.float32)
    for j, gene in enumerate(valid_genes):
        expr = extract_gene_expression(adata, gene, layer=layer, use_cache=True, dtype=np.float32)
        gene_matrix[:, j] = expr[row_indices]

    expr_df = pd.DataFrame(gene_matrix, columns=valid_genes)
    if transformation and transformation not in ("None", "none"):
        expr_df = apply_transformation(expr_df, method=transformation, copy=True)
    return expr_df


def _aggregate_by_group(expr_df, group_values, groups_to_process):
    group_cat = pd.Categorical(group_values, categories=groups_to_process, ordered=True)
    aggregated_data = expr_df.groupby(group_cat, observed=True, sort=False).mean()
    fraction_expressing = (expr_df > 0).groupby(group_cat, observed=True, sort=False).mean()
    aggregated_data.index.name = None
    fraction_expressing.index.name = None
    return aggregated_data, fraction_expressing


def _normalize_standardization(standardization):
    if standardization in ("across_cells", "across_groups", "z_score", "zscore"):
        return "zscore"
    return standardization


def _standardize_aggregated_data(aggregated_data, expr_df, standardization):
    standardization = _normalize_standardization(standardization)
    diverging = standardization in ("zscore", "group")

    if standardization == "zscore":
        mu = expr_df.mean(axis=0)
        sd = expr_df.std(axis=0).replace(0, 1.0)
        aggregated_data = (aggregated_data - mu) / sd
    elif standardization in ("minmax", "var"):
        lo = aggregated_data.min(axis=0)
        rng = (aggregated_data.max(axis=0) - lo).replace(0, 1.0)
        aggregated_data = (aggregated_data - lo) / rng
    elif standardization == "group":
        row_mean = aggregated_data.mean(axis=1)
        row_std = aggregated_data.std(axis=1).replace(0, 1.0)
        aggregated_data = aggregated_data.sub(row_mean, axis=0).div(row_std, axis=0)

    return aggregated_data, standardization, diverging


def _base_groups_and_genes(aggregated_data, valid_genes, selected_labels):
    if selected_labels:
        present_groups = set(aggregated_data.index)
        base_groups = [label for label in selected_labels if label in present_groups]
    else:
        base_groups = list(aggregated_data.index)
    return base_groups, list(valid_genes)


def _compute_linkage(X, method, metric):
    try:
        if X.shape[0] < 2:
            return None
        X = np.nan_to_num(X, nan=0.0)
        distances = pdist(X, metric=metric)
        if distances.size == 0 or np.allclose(distances, 0):
            return None
        return linkage(distances, method=method)
    except Exception:
        return None


def _cluster_orders(aggregated_data, base_genes, base_groups, cluster, method, metric):
    cluster_df = aggregated_data[base_genes].loc[base_groups]
    row_labels_base = list(cluster_df.index)
    col_labels_base = list(cluster_df.columns)
    row_Z = _compute_linkage(cluster_df.values, method, metric) if cluster in ('row', 'both') else None
    col_Z = _compute_linkage(cluster_df.values.T, method, metric) if cluster in ('col', 'both') else None

    group_order = (
        [row_labels_base[i] for i in leaves_list(row_Z)]
        if cluster in ('row', 'both') and row_Z is not None
        else base_groups
    )
    gene_order = (
        [col_labels_base[i] for i in leaves_list(col_Z)]
        if cluster in ('col', 'both') and col_Z is not None
        else base_genes
    )
    return group_order, gene_order, row_Z, col_Z, row_labels_base, col_labels_base


def _color_range(aggregated_data, gene_order, standardization):
    values = aggregated_data[gene_order].to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    vmin = float(finite_values.min()) if finite_values.size else 0.0
    vmax = float(finite_values.max()) if finite_values.size else 1.0

    if standardization in ("minmax", "var"):
        return 0.0, 1.0, "Scaled (0-1)"
    if standardization == "zscore":
        return -ZSCORE_COLOR_CLIP, ZSCORE_COLOR_CLIP, "Z-score"
    if standardization == "group":
        magnitude = max(abs(vmin), abs(vmax)) or 1.0
        return -magnitude, magnitude, "Z-score (group)"
    return vmin, vmax, "Mean Expression"


def _axis_context(gene_order, group_order, cluster, transpose, row_Z, col_Z):
    if transpose:
        x_items = group_order
        y_items = gene_order
        show_right_dendro = cluster in ('col', 'both') and len(gene_order) > 1 and col_Z is not None
        show_top_dendro = cluster in ('row', 'both') and len(group_order) > 1 and row_Z is not None
    else:
        x_items = gene_order
        y_items = group_order
        show_right_dendro = cluster in ('row', 'both') and len(group_order) > 1 and row_Z is not None
        show_top_dendro = cluster in ('col', 'both') and len(gene_order) > 1 and col_Z is not None

    x_categoryarray = x_items
    y_categoryarray = list(reversed(y_items))
    return {
        'x_items': x_items,
        'y_items': y_items,
        'x_categoryarray': x_categoryarray,
        'y_categoryarray': y_categoryarray,
        'show_right_dendro': show_right_dendro,
        'show_top_dendro': show_top_dendro,
        'right_dendro_items': y_categoryarray,
        'top_dendro_items': x_categoryarray,
    }


def _dot_size_scale(max_fraction):
    scale_params = (
        (150, [0.1, 0.075, 0.05, 0.025]) if max_fraction < 0.1 else
        (90, [0.2, 0.15, 0.1, 0.05]) if max_fraction < 0.2 else
        (40, [0.4, 0.3, 0.2, 0.1]) if max_fraction < 0.4 else
        (35, [0.5, 0.4, 0.3, 0.2]) if max_fraction < 0.5 else
        (30, [0.6, 0.5, 0.4, 0.3]) if max_fraction < 0.6 else
        (25, [0.7, 0.6, 0.5, 0.4]) if max_fraction < 0.7 else
        (20, [0.8, 0.6, 0.3, 0.1]) if max_fraction < 0.8 else
        (17, [1.0, 0.75, 0.50, 0.25])
    )
    scale, size_legend_values = scale_params
    return scale * DOT_SIZE_LEGEND_SCALE, size_legend_values


def _dot_hovertemplate(groupby, transpose=False):
    gene_line = 'Gene: %{customdata[2]}<br>' if transpose else 'Gene: %{x}<br>'
    group_line = f'{groupby}: %{{x}}<br>' if transpose else f'{groupby}: %{{y}}<br>'
    return (
        gene_line
        + group_line
        + 'Expression: %{customdata[0]:.4f}<br>'
        + 'Fraction: %{customdata[1]:.4f}<extra></extra>'
    )


def _map_to_axis(vals, leaves_labels, target_items):
    n = len(target_items)
    pos_map = {g: (i + 0.5) / n for i, g in enumerate(target_items)}
    leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
    leaf_pos = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
    out = []
    for v in vals:
        leaf_index = (v - _DENDRO_LEAF_BASE) / _DENDRO_LEAF_STEP
        p = max(0.0, min(leaf_index, len(leaf_pos) - 1))
        i0 = int(np.floor(p))
        i1 = min(i0 + 1, len(leaf_pos) - 1)
        frac = p - i0
        out.append(leaf_pos[i0] * (1 - frac) + leaf_pos[i1] * frac)
    return out


def _dendrogram_source(transpose, row_Z, col_Z, row_labels_base, col_labels_base, location):
    if location == "right":
        return (col_Z, col_labels_base) if transpose else (row_Z, row_labels_base)
    return (row_Z, row_labels_base) if transpose else (col_Z, col_labels_base)


def _add_dendrogram_trace(fig, Z, labels, orientation, target_items, xaxis, yaxis):
    try:
        if Z is None:
            raise ValueError("No linkage available")
        dendro = _scipy_dendro(Z, no_plot=True, orientation=orientation, labels=list(labels))
        max_h = max([max(dc) for dc in dendro['dcoord']]) or 1.0
        for ico, dco in zip(dendro['icoord'], dendro['dcoord']):
            mapped = _map_to_axis(ico, dendro['ivl'], target_items)
            heights = [h / max_h for h in dco]
            xs, ys = (heights, mapped) if orientation == 'right' else (mapped, heights)
            fig.add_trace(go.Scatter(
                x=xs, y=ys, xaxis=xaxis, yaxis=yaxis,
                mode='lines', line=dict(color='black', width=1),
                hoverinfo='skip', showlegend=False
            ))
    except Exception:
        pass


def _add_dendrograms(
    fig,
    show_right_dendro,
    show_top_dendro,
    transpose,
    row_Z,
    col_Z,
    row_labels_base,
    col_labels_base,
    right_dendro_items,
    top_dendro_items,
):
    if show_right_dendro:
        Z, labels = _dendrogram_source(
            transpose,
            row_Z,
            col_Z,
            row_labels_base,
            col_labels_base,
            "right",
        )
        _add_dendrogram_trace(fig, Z, labels, 'right', right_dendro_items, 'x3', 'y3')

    if show_top_dendro:
        Z, labels = _dendrogram_source(
            transpose,
            row_Z,
            col_Z,
            row_labels_base,
            col_labels_base,
            "top",
        )
        _add_dendrogram_trace(fig, Z, labels, 'top', top_dendro_items, 'x4', 'y4')


def plot_dot_matrix(
    adata, genes, groupby, selected_labels,
    transformation=None, standardization=None, layer=None,
    color_map='Viridis', plot_type='dotplot',
    cluster='none', method='average', metric='correlation',
    transpose=False, selected_cells=None
):
    color_map = resolve_continuous_colorscale(color_map)

    valid_genes = [gene for gene in genes if gene in adata.var_names]
    if not valid_genes:
        raise PreventUpdate

    # Scanpy-like flow: cached gene vectors + vectorized grouped aggregations.
    group_series = obs_col(adata.obs, groupby)
    row_indices, use_selected_cells = _selected_row_indices(
        adata,
        group_series,
        selected_labels,
        selected_cells,
    )
    if row_indices.size == 0:
        raise PreventUpdate

    group_values = group_series.iloc[row_indices]
    groups_to_process = _groups_to_process(group_values, selected_labels, use_selected_cells)

    if not groups_to_process:
        raise PreventUpdate

    expr_df = _expression_frame(adata, valid_genes, row_indices, layer, transformation)
    aggregated_data, fraction_expressing = _aggregate_by_group(expr_df, group_values, groups_to_process)

    # Standardization (per-gene / column):
    #   "zscore" -> (x-mean)/std over cells, applied to the per-group means, so each
    #               dot is the group's mean per-cell z-score -- matches the heatmap.
    #   "minmax" -> Scanpy standard_scale='var': per-gene min-max over the group
    #               means -> [0, 1].
    # Old aliases ('across_cells'/'across_groups') collapse to z-score; the legacy
    # notebook 'var' (min-max) and 'group' (row z-score) values are still honoured.
    aggregated_data, standardization, diverging = _standardize_aggregated_data(
        aggregated_data,
        expr_df,
        standardization,
    )

    # Figure out base lists
    base_groups, base_genes = _base_groups_and_genes(aggregated_data, valid_genes, selected_labels)

    # Optional hierarchical clustering to order rows/columns like Scanpy
    group_order, gene_order, row_Z, col_Z, row_labels_base, col_labels_base = _cluster_orders(
        aggregated_data,
        base_genes,
        base_groups,
        cluster,
        method,
        metric,
    )
    vmin, vmax, colorbar_title = _color_range(aggregated_data, gene_order, standardization)

    # Keep the visible top-to-bottom / left-to-right order aligned with user input.
    axis_context = _axis_context(gene_order, group_order, cluster, transpose, row_Z, col_Z)
    x_items = axis_context['x_items']
    y_items = axis_context['y_items']
    x_categoryarray = axis_context['x_categoryarray']
    y_categoryarray = axis_context['y_categoryarray']
    show_right_dendro = axis_context['show_right_dendro']
    show_top_dendro = axis_context['show_top_dendro']
    right_dendro_items = axis_context['right_dendro_items']
    top_dendro_items = axis_context['top_dendro_items']

    if plot_type == 'dotplot':
        df_expression = aggregated_data.reset_index().melt(id_vars=['index'], value_vars=gene_order, var_name='gene', value_name='expression')
        df_expression.rename(columns={'index': groupby}, inplace=True)
        df_fraction = fraction_expressing.reset_index().melt(id_vars=['index'], value_vars=gene_order, var_name='gene', value_name='fraction')
        df_fraction.rename(columns={'index': groupby}, inplace=True)
        df_merged = pd.merge(df_expression, df_fraction, on=[groupby, 'gene'])

        # Scale calculation from original
        max_fraction = df_merged['fraction'].max()
        scale, size_legend_values = _dot_size_scale(max_fraction)

        marker_sizes = (df_merged['fraction'] * scale)
        size_legend_sizes = [(s * scale)  for s in size_legend_values]
        custom_data = np.stack([df_merged['expression'], df_merged['fraction']], axis=-1)

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_merged[groupby] if transpose else df_merged['gene'],
            y=df_merged['gene'] if transpose else df_merged[groupby],
            mode='markers',
            showlegend=False,
            marker=dict(
                size=marker_sizes.tolist() if hasattr(marker_sizes, "tolist") else marker_sizes,
                color=df_merged['expression'].astype(float).tolist(),
                colorscale=color_map,
                cmin=vmin,
                cmax=vmax,
                line=dict(color='black', width=0.5),
                colorbar=_colorbar(colorbar_title)
            ),
            customdata=custom_data,
            hovertemplate=_dot_hovertemplate(groupby, transpose)
        ))
        
        # Fix hovertemplate gene label for transposed case
        if transpose:
            fig.data[0].customdata = np.stack([df_merged['expression'], df_merged['fraction'], df_merged['gene']], axis=-1)

        # Build a custom size legend as a right-side inset aligned with the colorbar
        y_positions = np.linspace(0.8, 0.2, len(size_legend_values))

        # Main axes domains
        main_x_right = 0.68 if show_right_dendro else 0.72
        main_y_top = 0.86 if show_top_dendro else 1.0

        layout_kwargs = dict(
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       domain=[0.0, main_x_right], categoryorder='array', categoryarray=x_categoryarray),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       categoryorder='array', categoryarray=y_categoryarray, domain=[0.0, main_y_top]),
            # Size legend (right-middle)
            xaxis2=dict(domain=[0.74, 0.96], range=[0.1, 0.9], autorange=False, fixedrange=True,
                        showgrid=False, zeroline=False, showticklabels=False, uirevision='frac-legend'),
            yaxis2=dict(domain=[0.2, 0.8], range=[0.1, 1.0], autorange=False, fixedrange=True,
                        showgrid=False, zeroline=False, showticklabels=False, uirevision='frac-legend'),
            margin=dict(r=260, t=60 if show_top_dendro else 20),
            plot_bgcolor='white', paper_bgcolor='white'
        )

        # Add dendrogram axes only when needed (prevents reserving empty space)
        layout_kwargs.update(_dendro_axis_layout(main_x_right, main_y_top, show_right_dendro, show_top_dendro))

        fig.update_layout(**layout_kwargs)

        # Title for size legend
        fig.add_trace(go.Scatter(
            x=[0.6], y=[0.95], xaxis='x2', yaxis='y2',
            mode='text', text=["Frac. cells"],
            textposition='top center',
            cliponaxis=False,
            showlegend=False, hoverinfo='skip',
            textfont=dict(size=11, color='black')
        ))

        # Dots with percent labels to the right. The dots are fixed-pixel markers;
        # seating them at x=0.4 (rather than hard against the legend's left edge)
        # keeps a real gap from the main plot, and cliponaxis clips them to this
        # subplot so on a very narrow figure they can't bleed onto the data grid.
        for size, value, y in zip(size_legend_sizes, size_legend_values, y_positions):
            percent = f"{int(round(value * 100))}%"
            fig.add_trace(go.Scatter(
                x=[0.40], y=[y], xaxis='x2', yaxis='y2',
                mode='markers',
                marker=dict(size=size, color='grey', line=dict(color='black', width=0.5)),
                cliponaxis=True,
                showlegend=False, hoverinfo='skip'
            ))
            fig.add_trace(go.Scatter(
                x=[0.80], y=[y], xaxis='x2', yaxis='y2',
                mode='text', text=[percent],
                textposition='middle left',
                cliponaxis=False,
                showlegend=False, hoverinfo='skip',
                textfont=dict(size=10, color='black')
            ))

    else:  # matrixplot
        z_data = aggregated_data.loc[group_order, gene_order].values
        if transpose:
            z_data = z_data.T
            
        fig = go.Figure(data=go.Heatmap(
            z=z_data,
            x=x_items,
            y=y_items,
            colorscale=color_map,
            zmid=0 if diverging else None,
            zmin=vmin,
            zmax=vmax,
            colorbar=_colorbar(colorbar_title),
            hovertemplate='%{y}<br>%{x}<br>Expression: %{z:.2f}<extra></extra>'
        ))

        # Dendrogram layout (same style as dotplot, but without the frac. cells inset)
        main_x_right = 0.84 if show_right_dendro else 0.96  # Leave extra space for colorbar when no row dendrogram
        main_y_top = 0.84 if show_top_dendro else 1.0

        layout_kwargs = dict(
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       domain=[0.0, main_x_right], categoryorder='array', categoryarray=x_categoryarray),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       categoryorder='array', categoryarray=y_categoryarray, domain=[0.0, main_y_top]),
            margin=dict(r=120 if show_right_dendro else 80, t=70 if show_top_dendro else 20,
                        b=100, l=100),
            plot_bgcolor='white', paper_bgcolor='white'
        )

        layout_kwargs.update(_dendro_axis_layout(main_x_right, main_y_top, show_right_dendro,
                                                 show_top_dendro, clamp_right_domain=True))

        fig.update_layout(**layout_kwargs)

    # ============ Dendrograms (Scanpy-like) ============
    _add_dendrograms(
        fig,
        show_right_dendro,
        show_top_dendro,
        transpose,
        row_Z,
        col_Z,
        row_labels_base,
        col_labels_base,
        right_dendro_items,
        top_dendro_items,
    )

    return fig
