import numpy as np
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
import pandas as pd
from .gene_extraction_utils import extract_gene_expression, apply_transformation
try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    _SCIPY_AVAILABLE = True
except Exception:
    _SCIPY_AVAILABLE = False

def plot_dot_matrix(
    adata, genes, groupby, selected_labels,
    aggregation='mean', transformation=None, standardization=None,
    vmin=None, vmax=None, expression_threshold=0,
    color_map='Viridis', plot_type='dotplot',
    cluster='none', method='average', metric='correlation'
):
    """Optimized dotplot function that processes data group by group to avoid loading all cells at once."""
    valid_genes = [gene for gene in genes if gene in adata.var_names]
    if not valid_genes:
        raise PreventUpdate

    # Determine which groups to process
    if selected_labels:
        groups_to_process = selected_labels
    else:
        groups_to_process = adata.obs[groupby].unique()
    
    # Initialize results dictionaries
    aggregated_results = {}
    fraction_results = {}
    gene_expressions = {}
    for gene in valid_genes:
        # Extract from original adata (not views) to ensure cache hits
        gene_expressions[gene] = extract_gene_expression(adata, gene)
    
    # Process each group separately using the cached gene data
    for group in groups_to_process:
        # Get indices for this group only
        group_mask = adata.obs[groupby] == group
        group_indices = np.where(group_mask)[0]
        
        if len(group_indices) == 0:
            continue
            
        group_expression_list = []
        
        for gene in valid_genes:
            gene_expr = gene_expressions[gene][group_indices]
            
            if transformation:
                gene_expr = apply_transformation(gene_expr, transformation, copy=False)
            agg_value = np.mean(gene_expr)
            
            fraction = np.mean(gene_expr > expression_threshold)
            
            group_expression_list.append({
                'gene': gene,
                'aggregated': agg_value,
                'fraction': fraction
            })
        
        # Store results for this group
        aggregated_results[group] = {gene_data['gene']: gene_data['aggregated'] for gene_data in group_expression_list}
        fraction_results[group] = {gene_data['gene']: gene_data['fraction'] for gene_data in group_expression_list}
    
    # Convert to DataFrames
    aggregated_data = pd.DataFrame(aggregated_results).T
    fraction_expressing = pd.DataFrame(fraction_results).T
    
    # Apply standardization if needed
    if standardization == 'var':
        # Standardize per gene (each column)
        aggregated_data = (aggregated_data - aggregated_data.min()) / (aggregated_data.max() - aggregated_data.min())
    elif standardization == 'group':
        # Standardize per group (each row)
        aggregated_data = (aggregated_data.sub(aggregated_data.mean(axis=1), axis=0)
                                        .div(aggregated_data.std(axis=1), axis=0))
    
    # Figure out base lists
    if selected_labels:
        base_groups = [label for label in aggregated_data.index if label in selected_labels]
    else:
        base_groups = list(aggregated_data.index)
    base_genes = [g for g in valid_genes]

    # Optional hierarchical clustering to order rows/columns like Scanpy
    def _cluster_order(df_vals, axis='row'):
        if not _SCIPY_AVAILABLE:
            return None
        try:
            X = df_vals.values if axis == 'row' else df_vals.values.T
            if X.shape[0] < 2:
                return list(df_vals.index if axis == 'row' else df_vals.columns)
            X = np.nan_to_num(X, nan=0.0)
            D = pdist(X, metric=metric)
            if np.allclose(D, 0):
                return list(df_vals.index if axis == 'row' else df_vals.columns)
            Z = linkage(D, method=method)
            order = leaves_list(Z)
            labels = df_vals.index if axis == 'row' else df_vals.columns
            return [labels[i] for i in order]
        except Exception:
            return None

    if cluster in ('row', 'both'):
        groups = _cluster_order(aggregated_data[base_genes].loc[base_groups], axis='row') or base_groups
    else:
        # Keep previous reverse default when not clustering
        groups = base_groups[::-1]

    if cluster in ('col', 'both'):
        ordered_genes = _cluster_order(aggregated_data[base_genes].loc[base_groups], axis='col') or base_genes
        valid_genes = ordered_genes
    # else keep valid_genes as is
    
    vmin = vmin if vmin is not None else float(aggregated_data[valid_genes].min().min())
    vmax = vmax if vmax is not None else float(aggregated_data[valid_genes].max().max())

    if plot_type == 'dotplot':
        df_expression = aggregated_data.reset_index().melt(id_vars=['index'], value_vars=valid_genes, var_name='gene', value_name='expression')
        df_expression.rename(columns={'index': groupby}, inplace=True)
        df_fraction = fraction_expressing.reset_index().melt(id_vars=['index'], value_vars=valid_genes, var_name='gene', value_name='fraction')
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
            x=df_merged['gene'],
            y=df_merged[groupby],
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
                    title=f'{aggregation.capitalize()} Expression ({transformation})' if transformation and transformation != 'None' else f'{aggregation.capitalize()} Expression',
                    tickfont=dict(color='DarkSlateGrey', size=10),
                    len=0.6,
                    yanchor="middle",
                    y=0.5,
                    x=0.98
                )
            ),
            customdata=custom_data,
            hovertemplate=(
                'Gene: %{x}<br>'
                f'{groupby}: %{{y}}<br>'
                'Expression: %{customdata[0]:.4f}<br>'
                'Fraction: %{customdata[1]:.4f}<extra></extra>'
            )
        ))
        # Build a custom size legend as a right-side inset aligned with the colorbar
        y_positions = np.linspace(0.8, 0.2, len(size_legend_values))

        # Determine if dendrograms are shown to allocate space dynamically
        show_row_dendro = (_SCIPY_AVAILABLE and cluster in ('row', 'both') and len(groups) > 1)
        show_col_dendro = (_SCIPY_AVAILABLE and cluster in ('col', 'both') and len(valid_genes) > 1)

        # Main axes domains
        main_x_right = 0.68 if show_row_dendro else 0.72
        main_y_top = 0.86 if show_col_dendro else 1.0

        layout_kwargs = dict(
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       domain=[0.0, main_x_right], categoryorder='array', categoryarray=valid_genes),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       categoryorder='array', categoryarray=groups, domain=[0.0, main_y_top]),
            # Size legend (right-middle)
            xaxis2=dict(domain=[0.74, 0.96], range=[0.1, 0.9], autorange=False, fixedrange=True,
                        showgrid=False, zeroline=False, showticklabels=False, uirevision='frac-legend'),
            yaxis2=dict(domain=[0.2, 0.8], range=[0.1, 1.0], autorange=False, fixedrange=True,
                        showgrid=False, zeroline=False, showticklabels=False, uirevision='frac-legend'),
            margin=dict(r=260, t=60 if show_col_dendro else 20),
            plot_bgcolor='white', paper_bgcolor='white'
        )

        # Add dendrogram axes only when needed (prevents reserving empty space)
        if show_row_dendro:
            layout_kwargs.update({
                'xaxis3': dict(domain=[main_x_right + 0.01, main_x_right + 0.04], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis3': dict(domain=[0.0, main_y_top], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })
        if show_col_dendro:
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
            textposition='top left',
            cliponaxis=False,
            showlegend=False, hoverinfo='skip',
            textfont=dict(size=11, color='black')
        ))

        # Dots with percent labels to the right
        for size, value, y in zip(size_legend_sizes, size_legend_values, y_positions):
            percent = f"{int(round(value * 100))}%"
            fig.add_trace(go.Scatter(
                x=[0.25], y=[y], xaxis='x2', yaxis='y2',
                mode='markers',
                marker=dict(size=size, color='grey', line=dict(color='black', width=0.5)),
                cliponaxis=False,
                showlegend=False, hoverinfo='skip'
            ))
            fig.add_trace(go.Scatter(
                x=[0.62], y=[y], xaxis='x2', yaxis='y2',
                mode='text', text=[percent],
                textposition='middle left',
                cliponaxis=False,
                showlegend=False, hoverinfo='skip',
                textfont=dict(size=10, color='black')
            ))

        # ============ Dendrograms (Scanpy-like) ============
        if _SCIPY_AVAILABLE:
            # Helper to map scipy dendrogram coords to axis coords (rows)
            def _map_row_y(yy_vals, leaves_labels):
                # Display positions for current groups order
                n = len(groups)
                pos_map = {g: (i + 0.5) / n for i, g in enumerate(groups)}
                # Dendrogram leaf order labels
                leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
                leaf_pos_y = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
                out = []
                for yy in yy_vals:
                    p = max(0.0, min((yy - 5.0) / 10.0, len(leaf_pos_y) - 1))
                    i0 = int(np.floor(p))
                    i1 = min(i0 + 1, len(leaf_pos_y) - 1)
                    frac = p - i0
                    yy_mapped = leaf_pos_y[i0] * (1 - frac) + leaf_pos_y[i1] * frac
                    out.append(yy_mapped)
                return out

            # Helper to map scipy dendrogram coords to axis coords (cols)
            def _map_col_x(xx_vals, leaves_labels):
                n = len(valid_genes)
                pos_map = {g: (i + 0.5) / n for i, g in enumerate(valid_genes)}
                leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
                leaf_pos_x = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
                out = []
                for xx in xx_vals:
                    p = max(0.0, min((xx - 5.0) / 10.0, len(leaf_pos_x) - 1))
                    i0 = int(np.floor(p))
                    i1 = min(i0 + 1, len(leaf_pos_x) - 1)
                    frac = p - i0
                    xx_mapped = leaf_pos_x[i0] * (1 - frac) + leaf_pos_x[i1] * frac
                    out.append(xx_mapped)
                return out

            from scipy.cluster.hierarchy import dendrogram as _scipy_dendro

            # Row dendrogram
            if show_row_dendro:
                try:
                    row_mat = aggregated_data[valid_genes].loc[groups].values
                    row_mat = np.nan_to_num(row_mat, nan=0.0)
                    Zr = linkage(pdist(row_mat, metric=metric), method=method)
                    drow = _scipy_dendro(Zr, no_plot=True, orientation='right', labels=list(groups))
                    max_h = max([max(dc) for dc in drow['dcoord']]) or 1.0
                    for ico, dco in zip(drow['icoord'], drow['dcoord']):
                        ys = _map_row_y(ico, drow['ivl'])  # ivl are leaf labels in order
                        xs = [h / max_h for h in dco]
                        fig.add_trace(go.Scatter(
                            x=xs, y=ys, xaxis='x3', yaxis='y3',
                            mode='lines', line=dict(color='black', width=1),
                            hoverinfo='skip', showlegend=False
                        ))
                except Exception:
                    pass

            # Column dendrogram
            if show_col_dendro:
                try:
                    col_mat = aggregated_data[valid_genes].loc[groups].values.T
                    col_mat = np.nan_to_num(col_mat, nan=0.0)
                    Zc = linkage(pdist(col_mat, metric=metric), method=method)
                    dcol = _scipy_dendro(Zc, no_plot=True, orientation='top', labels=list(valid_genes))
                    max_hc = max([max(dc) for dc in dcol['dcoord']]) or 1.0
                    for ico, dco in zip(dcol['icoord'], dcol['dcoord']):
                        xs = _map_col_x(ico, dcol['ivl'])
                        ys = [h / max_hc for h in dco]
                        fig.add_trace(go.Scatter(
                            x=xs, y=ys, xaxis='x4', yaxis='y4',
                            mode='lines', line=dict(color='black', width=1),
                            hoverinfo='skip', showlegend=False
                        ))
                except Exception:
                    pass

    else:  # matrixplot
        fig = go.Figure(data=go.Heatmap(
            z=aggregated_data.loc[groups, valid_genes].values,
            x=valid_genes,
            y=groups,
            colorscale=color_map,
            zmid=None,
            zmin=vmin,
            zmax=vmax,
            colorbar=dict(
                title=f'{aggregation.capitalize()} Expression ({transformation})' if transformation and transformation != 'None' else f'{aggregation.capitalize()} Expression',
                tickfont=dict(color='DarkSlateGrey', size=10),
                len=0.6,
                yanchor="middle",
                y=0.5,
                x=0.98
            ),
            hovertemplate='%{y}<br>%{x}<br>Expression: %{z:.2f}<extra></extra>'
        ))

        # Dendrogram layout (same style as dotplot, but without the frac. cells inset)
        show_row_dendro = (_SCIPY_AVAILABLE and cluster in ('row', 'both') and len(groups) > 1)
        show_col_dendro = (_SCIPY_AVAILABLE and cluster in ('col', 'both') and len(valid_genes) > 1)

        main_x_right = 0.84 if show_row_dendro else 0.96  # Leave extra space for colorbar when no row dendrogram
        main_y_top = 0.84 if show_col_dendro else 1.0

        layout_kwargs = dict(
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       domain=[0.0, main_x_right], categoryorder='array', categoryarray=valid_genes),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False,
                       categoryorder='array', categoryarray=groups, domain=[0.0, main_y_top]),
            margin=dict(r=120 if show_row_dendro else 80, t=70 if show_col_dendro else 20,
                        b=100, l=100),
            plot_bgcolor='white', paper_bgcolor='white'
        )

        if show_row_dendro:
            layout_kwargs.update({
                'xaxis3': dict(domain=[main_x_right + 0.01, min(main_x_right + 0.04, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis3': dict(domain=[0.0, main_y_top], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })
        if show_col_dendro:
            layout_kwargs.update({
                'xaxis4': dict(domain=[0.0, main_x_right], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
                'yaxis4': dict(domain=[main_y_top + 0.02, min(main_y_top + 0.14, 0.99)], range=[0, 1],
                               showgrid=False, zeroline=False, showticklabels=False),
            })

        fig.update_layout(**layout_kwargs)

        # Draw dendrograms
        if _SCIPY_AVAILABLE:
            from scipy.cluster.hierarchy import dendrogram as _scipy_dendro

            def _map_row_y(yy_vals, leaves_labels):
                n = len(groups)
                pos_map = {g: (i + 0.5) / n for i, g in enumerate(groups)}
                leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
                leaf_pos_y = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
                out = []
                for yy in yy_vals:
                    p = max(0.0, min((yy - 5.0) / 10.0, len(leaf_pos_y) - 1))
                    i0 = int(np.floor(p))
                    i1 = min(i0 + 1, len(leaf_pos_y) - 1)
                    frac = p - i0
                    out.append(leaf_pos_y[i0] * (1 - frac) + leaf_pos_y[i1] * frac)
                return out

            def _map_col_x(xx_vals, leaves_labels):
                n = len(valid_genes)
                pos_map = {g: (i + 0.5) / n for i, g in enumerate(valid_genes)}
                leaf_pos_labels = [leaves_labels[i] for i in range(len(leaves_labels))]
                leaf_pos_x = [pos_map.get(lbl, 0.0) for lbl in leaf_pos_labels]
                out = []
                for xx in xx_vals:
                    p = max(0.0, min((xx - 5.0) / 10.0, len(leaf_pos_x) - 1))
                    i0 = int(np.floor(p))
                    i1 = min(i0 + 1, len(leaf_pos_x) - 1)
                    frac = p - i0
                    out.append(leaf_pos_x[i0] * (1 - frac) + leaf_pos_x[i1] * frac)
                return out

            if show_row_dendro:
                try:
                    row_mat = aggregated_data[valid_genes].loc[groups].values
                    row_mat = np.nan_to_num(row_mat, nan=0.0)
                    Zr = linkage(pdist(row_mat, metric=metric), method=method)
                    drow = _scipy_dendro(Zr, no_plot=True, orientation='right', labels=list(groups))
                    max_h = max([max(dc) for dc in drow['dcoord']]) or 1.0
                    for ico, dco in zip(drow['icoord'], drow['dcoord']):
                        ys = _map_row_y(ico, drow['ivl'])
                        xs = [h / max_h for h in dco]
                        fig.add_trace(go.Scatter(
                            x=xs, y=ys, xaxis='x3', yaxis='y3',
                            mode='lines', line=dict(color='black', width=1),
                            hoverinfo='skip', showlegend=False
                        ))
                except Exception:
                    pass

            if show_col_dendro:
                try:
                    col_mat = aggregated_data[valid_genes].loc[groups].values.T
                    col_mat = np.nan_to_num(col_mat, nan=0.0)
                    Zc = linkage(pdist(col_mat, metric=metric), method=method)
                    dcol = _scipy_dendro(Zc, no_plot=True, orientation='top', labels=list(valid_genes))
                    max_hc = max([max(dc) for dc in dcol['dcoord']]) or 1.0
                    for ico, dco in zip(dcol['icoord'], dcol['dcoord']):
                        xs = _map_col_x(ico, dcol['ivl'])
                        ys = [h / max_hc for h in dco]
                        fig.add_trace(go.Scatter(
                            x=xs, y=ys, xaxis='x4', yaxis='y4',
                            mode='lines', line=dict(color='black', width=1),
                            hoverinfo='skip', showlegend=False
                        ))
                except Exception:
                    pass

    return fig
