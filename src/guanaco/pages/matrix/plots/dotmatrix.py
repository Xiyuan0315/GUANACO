import numpy as np
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
import pandas as pd
from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation
from guanaco.pages.matrix.plots.heatmap import ZSCORE_COLOR_CLIP

from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as _scipy_dendro
from scipy.spatial.distance import pdist


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
    group_series = adata.obs[groupby]
    use_selected_cells = selected_cells is not None and len(selected_cells) > 0
    if use_selected_cells:
        selected_pos = adata.obs_names.get_indexer(selected_cells)
        row_indices = selected_pos[selected_pos >= 0].astype(np.int64, copy=False)
    elif selected_labels:
        label_mask = group_series.isin(selected_labels).to_numpy()
        row_indices = np.flatnonzero(label_mask).astype(np.int64, copy=False)
    else:
        row_indices = np.arange(adata.n_obs, dtype=np.int64)

    if row_indices.size == 0:
        raise PreventUpdate

    group_values = group_series.iloc[row_indices]
    if selected_labels and not use_selected_cells:
        present_groups = set(pd.unique(group_values))
        groups_to_process = [g for g in selected_labels if g in present_groups]
    else:
        groups_to_process = list(pd.unique(group_values))

    if not groups_to_process:
        raise PreventUpdate

    n_cells = row_indices.size
    gene_matrix = np.empty((n_cells, len(valid_genes)), dtype=np.float32)
    for j, gene in enumerate(valid_genes):
        expr = extract_gene_expression(adata, gene, layer=layer, use_cache=True, dtype=np.float32)
        gene_matrix[:, j] = expr[row_indices]
    expr_df = pd.DataFrame(gene_matrix, columns=valid_genes)

    if transformation and transformation not in ("None", "none"):
        expr_df = apply_transformation(expr_df, method=transformation, copy=True)

    group_cat = pd.Categorical(group_values, categories=groups_to_process, ordered=True)
    aggregated_data = expr_df.groupby(group_cat, observed=True, sort=False).mean()
    fraction_expressing = (expr_df > 0).groupby(group_cat, observed=True, sort=False).mean()
    aggregated_data.index.name = None
    fraction_expressing.index.name = None

    # Standardization (per-gene / column):
    #   "zscore" -> (x-mean)/std over cells, applied to the per-group means, so each
    #               dot is the group's mean per-cell z-score -- matches the heatmap.
    #   "minmax" -> Scanpy standard_scale='var': per-gene min-max over the group
    #               means -> [0, 1].
    # Old aliases ('across_cells'/'across_groups') collapse to z-score; the legacy
    # notebook 'var' (min-max) and 'group' (row z-score) values are still honoured.
    if standardization in ("across_cells", "across_groups", "z_score", "zscore"):
        standardization = "zscore"
    diverging = standardization in ("zscore", "group")  # zero-centred colour scale
    if standardization == "zscore":
        mu = expr_df.mean(axis=0)
        sd = expr_df.std(axis=0).replace(0, 1.0)
        aggregated_data = (aggregated_data - mu) / sd
    elif standardization in ("minmax", "var"):
        lo = aggregated_data.min(axis=0)
        rng = (aggregated_data.max(axis=0) - lo).replace(0, 1.0)
        aggregated_data = (aggregated_data - lo) / rng
    elif standardization == "group":  # legacy: z-score per group (row)
        aggregated_data = (aggregated_data.sub(aggregated_data.mean(axis=1), axis=0)
                                        .div(aggregated_data.std(axis=1), axis=0))

    # Figure out base lists
    if selected_labels:
        present_groups = set(aggregated_data.index)
        base_groups = [label for label in selected_labels if label in present_groups]
    else:
        base_groups = list(aggregated_data.index)
    base_genes = [g for g in valid_genes]

    # Optional hierarchical clustering to order rows/columns like Scanpy
    def _compute_linkage(X):

        try:
            if X.shape[0] < 2:
                return None
            X = np.nan_to_num(X, nan=0.0)
            D = pdist(X, metric=metric)
            if D.size == 0 or np.allclose(D, 0):
                return None
            return linkage(D, method=method)
        except Exception:
            return None

    cluster_df = aggregated_data[base_genes].loc[base_groups]
    row_labels_base = list(cluster_df.index)
    col_labels_base = list(cluster_df.columns)
    row_Z = _compute_linkage(cluster_df.values) if cluster in ('row', 'both') else None
    col_Z = _compute_linkage(cluster_df.values.T) if cluster in ('col', 'both') else None

    if cluster in ('row', 'both'):
        group_order = [row_labels_base[i] for i in leaves_list(row_Z)] if row_Z is not None else base_groups
    else:
        group_order = base_groups

    if cluster in ('col', 'both'):
        gene_order = [col_labels_base[i] for i in leaves_list(col_Z)] if col_Z is not None else base_genes
    else:
        gene_order = base_genes

    vmin = float(aggregated_data[gene_order].min().min())
    vmax = float(aggregated_data[gene_order].max().max())
    if standardization in ("minmax", "var"):
        # Already in [0, 1] (Scanpy standard_scale), sequential scale.
        vmin, vmax = 0.0, 1.0
        colorbar_title = "Scaled (0–1)"
    elif standardization == "zscore":
        # Cap the symmetric range so outlier z-scores don't wash out the colours
        # (same ±cap as the heatmap, so a z-score reads identically in both).
        vmin, vmax = -ZSCORE_COLOR_CLIP, ZSCORE_COLOR_CLIP
        colorbar_title = "Z-score"
    elif standardization == "group":  # legacy row z-score: data-driven symmetric range
        m = max(abs(vmin), abs(vmax)) or 1.0
        vmin, vmax = -m, m
        colorbar_title = "Z-score (group)"
    else:
        colorbar_title = "Mean Expression"

    # Keep the visible top-to-bottom / left-to-right order aligned with user input.
    if transpose:
        x_items = group_order
        y_items = gene_order
        x_categoryarray = x_items
        y_categoryarray = list(reversed(y_items))
        # Dendrogram flags
        show_right_dendro = (cluster in ('col', 'both') and len(gene_order) > 1)
        show_top_dendro = (cluster in ('row', 'both') and len(group_order) > 1)
        right_dendro_items = y_categoryarray
        top_dendro_items = x_categoryarray
    else:
        x_items = gene_order
        y_items = group_order
        x_categoryarray = x_items
        y_categoryarray = list(reversed(y_items))
        # Dendrogram flags
        show_right_dendro = (cluster in ('row', 'both') and len(group_order) > 1)
        show_top_dendro = (cluster in ('col', 'both') and len(gene_order) > 1)
        right_dendro_items = y_categoryarray
        top_dendro_items = x_categoryarray

    if plot_type == 'dotplot':
        df_expression = aggregated_data.reset_index().melt(id_vars=['index'], value_vars=gene_order, var_name='gene', value_name='expression')
        df_expression.rename(columns={'index': groupby}, inplace=True)
        df_fraction = fraction_expressing.reset_index().melt(id_vars=['index'], value_vars=gene_order, var_name='gene', value_name='fraction')
        df_fraction.rename(columns={'index': groupby}, inplace=True)
        df_merged = pd.merge(df_expression, df_fraction, on=[groupby, 'gene'])

        # Scale calculation from original
        max_fraction = df_merged['fraction'].max()
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
        # Increase scale for better visibility
        scale *= 1.6

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
                colorbar=dict(
                    title=colorbar_title,
                    tickfont=dict(color='DarkSlateGrey', size=10),
                    len=0.6,
                    yanchor="middle",
                    y=0.5,
                    x=0.98
                )
            ),
            customdata=custom_data,
            hovertemplate=(
                'Gene: %{customdata[2]}<br>' if transpose else 'Gene: %{x}<br>'
                f'{groupby}: %{{x}}<br>' if transpose else f'{groupby}: %{{y}}<br>'
                'Expression: %{customdata[0]:.4f}<br>'
                'Fraction: %{customdata[1]:.4f}<extra></extra>'
            )
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
        if show_right_dendro:
            layout_kwargs.update({
                'xaxis3': dict(domain=[main_x_right + 0.01, main_x_right + 0.04], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis3': dict(domain=[0.0, main_y_top], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })
        if show_top_dendro:
            layout_kwargs.update({
                'xaxis4': dict(domain=[0.0, main_x_right], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis4': dict(domain=[main_y_top + 0.02, min(main_y_top + 0.14, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })

        fig.update_layout(**layout_kwargs)

        # Title for size legend
        fig.add_trace(go.Scatter(
            x=[0.5], y=[0.95], xaxis='x2', yaxis='y2',
            mode='text', text=["Frac. cells"],
            textposition='top center',
            cliponaxis=False,
            showlegend=False, hoverinfo='skip',
            textfont=dict(size=11, color='black')
        ))

        # Dots with percent labels to the right
        for size, value, y in zip(size_legend_sizes, size_legend_values, y_positions):
            percent = f"{int(round(value * 100))}%"
            fig.add_trace(go.Scatter(
                x=[0.20], y=[y], xaxis='x2', yaxis='y2',
                mode='markers',
                marker=dict(size=size, color='grey', line=dict(color='black', width=0.5)),
                cliponaxis=False,
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
            colorbar=dict(
                title=colorbar_title,
                tickfont=dict(color='DarkSlateGrey', size=10),
                len=0.6,
                yanchor="middle",
                y=0.5,
                x=0.98
            ),
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

        if show_right_dendro:
            layout_kwargs.update({
                'xaxis3': dict(domain=[main_x_right + 0.01, min(main_x_right + 0.04, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis3': dict(domain=[0.0, main_y_top], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })
        if show_top_dendro:
            layout_kwargs.update({
                'xaxis4': dict(domain=[0.0, main_x_right], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis4': dict(domain=[main_y_top + 0.02, min(main_y_top + 0.14, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })

        fig.update_layout(**layout_kwargs)

    # ============ Dendrograms (Scanpy-like) ============
    def _map_to_axis(vals, leaves_labels, target_items):
        n = len(target_items)
        pos_map = {g: (i + 0.5) / n for i, g in enumerate(target_items)}
        leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
        leaf_pos = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
        out = []
        for v in vals:
            p = max(0.0, min((v - 5.0) / 10.0, len(leaf_pos) - 1))
            i0 = int(np.floor(p))
            i1 = min(i0 + 1, len(leaf_pos) - 1)
            frac = p - i0
            out.append(leaf_pos[i0] * (1 - frac) + leaf_pos[i1] * frac)
        return out

    # Right dendrogram (x3, y3)
    if show_right_dendro:
        try:
            if transpose:
                Z = col_Z
                labels = col_labels_base
            else:
                Z = row_Z
                labels = row_labels_base
            if Z is None:
                raise ValueError("No linkage available")
            dendro = _scipy_dendro(Z, no_plot=True, orientation='right', labels=list(labels))
            max_h = max([max(dc) for dc in dendro['dcoord']]) or 1.0
            for ico, dco in zip(dendro['icoord'], dendro['dcoord']):
                ys = _map_to_axis(ico, dendro['ivl'], right_dendro_items)
                xs = [h / max_h for h in dco]
                fig.add_trace(go.Scatter(
                    x=xs, y=ys, xaxis='x3', yaxis='y3',
                    mode='lines', line=dict(color='black', width=1),
                    hoverinfo='skip', showlegend=False
                ))
        except Exception:
            pass

    # Top dendrogram (x4, y4)
    if show_top_dendro:
        try:
            if transpose:
                Z = row_Z
                labels = row_labels_base
            else:
                Z = col_Z
                labels = col_labels_base
            if Z is None:
                raise ValueError("No linkage available")
            dendro = _scipy_dendro(Z, no_plot=True, orientation='top', labels=list(labels))
            max_h = max([max(dc) for dc in dendro['dcoord']]) or 1.0
            for ico, dco in zip(dendro['icoord'], dendro['dcoord']):
                xs = _map_to_axis(ico, dendro['ivl'], top_dendro_items)
                ys = [h / max_h for h in dco]
                fig.add_trace(go.Scatter(
                    x=xs, y=ys, xaxis='x4', yaxis='y4',
                    mode='lines', line=dict(color='black', width=1),
                    hoverinfo='skip', showlegend=False
                ))
        except Exception:
            pass

    return fig
