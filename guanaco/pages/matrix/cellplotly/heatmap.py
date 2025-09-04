import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from guanaco.pages.matrix.cellplotly.gene_extraction_utils import (
    extract_gene_expression, extract_multiple_genes, apply_transformation, bin_cells_for_heatmap
)

# ==========================================================
# Helper functions

def _validate_genes(adata, genes):
    return [g for g in genes if g in adata.var_names]


def _map_transformation_args(transformation, log, z_score):
    if transformation:
        if transformation == 'log':
            log = True
        elif transformation in ['z_score', 'zscore']:
            z_score = True
    return log, z_score


def _is_continuous_annotation(adata, annotation, threshold=50):
    if annotation not in adata.obs.columns:
        return False
    dtype = str(adata.obs[annotation].dtype)
    if any(t in dtype for t in ['float', 'int']):
        return adata.obs[annotation].nunique() >= threshold
    return False


def _filter_cells_and_obs(adata, groupby1, labels):
    is_backed = hasattr(adata, 'isbacked') and adata.isbacked
    original_adata = adata
    filtered_obs = adata.obs
    filtered_obs_names = adata.obs_names
    cell_indices_array = None
    if labels:
        mask = adata.obs[groupby1].isin(labels)
        if is_backed:
            cell_indices_array = np.where(mask)[0]
            filtered_obs = adata.obs.iloc[cell_indices_array]
            filtered_obs_names = adata.obs_names[cell_indices_array]
        else:
            adata = adata[mask].copy()
            filtered_obs = adata.obs
            filtered_obs_names = adata.obs_names
    return {
        'adata': adata,
        'original_adata': original_adata,
        'filtered_obs': filtered_obs,
        'filtered_obs_names': filtered_obs_names,
        'cell_indices_array': cell_indices_array,
        'is_backed': is_backed,
    }


def _extract_gene_df_after_filter(ctx, valid_genes):
    adata = ctx['adata']
    original_adata = ctx['original_adata']
    cell_indices_array = ctx['cell_indices_array']
    filtered_obs_names = ctx['filtered_obs_names']
    is_backed = ctx['is_backed']
    if valid_genes is None or len(valid_genes) == 0:
        return pd.DataFrame()
    if ctx['cell_indices_array'] is not None and is_backed:
        cols = []
        for g in valid_genes:
            expr = extract_gene_expression(original_adata, g)
            cols.append(pd.Series(expr[cell_indices_array], name=g, index=filtered_obs_names))
        gene_df = pd.concat(cols, axis=1)
    else:
        gene_df = extract_multiple_genes(adata, valid_genes)
    gene_df.insert(0, 'CellID', filtered_obs_names)
    return gene_df


def _apply_transformations(gene_df, valid_genes, log, z_score):
    if log:
        for g in valid_genes:
            gene_df[g] = apply_transformation(gene_df[g], method='log1p')
    if z_score:
        for g in valid_genes:
            gene_df[g] = apply_transformation(gene_df[g], method='z_score')
    return gene_df


def _make_expression_heatmap(matrix, genes, color_map, log, z_score, colorbar_len):
    if z_score:
        z_max = max(abs(matrix.min()), abs(matrix.max()))
        zmin, zmax, zmid = -z_max, z_max, 0
    else:
        zmin, zmax, zmid = matrix.min(), matrix.max(), None
    return go.Heatmap(
        z=matrix,
        x=list(range(matrix.shape[1])),
        y=list(range(len(genes))),
        colorscale=color_map,
        hoverinfo='skip',
        zmin=zmin,
        zmax=zmax,
        zmid=zmid,
        colorbar=dict(
            title=f'Expression(log)' if log else f'Expression(z-score)' if z_score else 'Expression',
            len=colorbar_len,
            y=1,
            yanchor='top'
        ),
    )


def _make_annotation_bar(groups, group_sizes, color_map, width=1):
    traces = []
    for i, label in enumerate(groups):
        traces.append(go.Bar(
            x=[group_sizes[i]],
            marker_color=(color_map.get(label, 'grey') if isinstance(color_map, dict) else 'grey'),
            name=f"{label}",
            hovertemplate=f'<b>{label}</b><br>Count: {group_sizes[i]}<extra></extra>',
            orientation='h',
            showlegend=False,
            width=width,
        ))
    return traces


def _add_boundaries(fig, group_sizes, row=1, col=1, width=1, n_genes=None):
    positions = np.cumsum(group_sizes[:-1]).tolist()
    y0 = -0.5
    y1 = (n_genes - 0.5) if n_genes is not None else y0
    for pos in positions:
        fig.add_shape(
            type="line",
            x0=pos, y0=y0, x1=pos, y1=y1,
            xref=f"x{col}", yref=f"y{row}",
            line=dict(color="rgba(0, 0, 0, 0.8)", width=width)
        )


def _default_color_maps(adata_obs, original_adata, groupby1, groupby2, groupby1_label_color_map, groupby2_label_color_map):
    from guanaco.data_loader import color_config
    tab20_colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
                    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
                    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
                    '#17becf', '#9edae5']
    if groupby1_label_color_map is None:
        unique_labels_primary = sorted(adata_obs[groupby1].unique()) if adata_obs is not None else sorted(original_adata.obs[groupby1].unique())
        groupby1_label_color_map = {label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels_primary)}
    if groupby2 and groupby2_label_color_map is None:
        unique_labels_secondary = sorted(adata_obs[groupby2].unique()) if adata_obs is not None else sorted(original_adata.obs[groupby2].unique())
        groupby2_label_color_map = dict(zip(unique_labels_secondary, tab20_colors[:len(unique_labels_secondary)]))
    return groupby1_label_color_map, groupby2_label_color_map


def _legend_annotations(groupby1, label_list1, groupby1_label_color_map, has_secondary, groupby2, label2_dict, groupby2_label_color_map):
    legend_annotations = []
    legend_start_x = 1.01
    current_y = 0.5
    legend_annotations.append(dict(
        x=legend_start_x, y=current_y,
        xref='paper', yref='paper', text=f"<b>{groupby1}</b>", showarrow=False,
        font=dict(size=12, color='black'), xanchor='left', yanchor='top'))
    current_y -= 0.08
    for label in [l for l in label_list1 if pd.notna(l)]:
        color = groupby1_label_color_map.get(label, 'grey')
        legend_annotations.append(dict(
            x=legend_start_x, y=current_y, xref='paper', yref='paper',
            text=f"<span style='color:{color}'>■</span> {label}",
            showarrow=False, font=dict(size=12), xanchor='left', yanchor='middle'))
        current_y -= 0.04
    if has_secondary:
        current_y -= 0.04
        legend_annotations.append(dict(
            x=legend_start_x, y=current_y,
            xref='paper', yref='paper', text=f"<b>{groupby2}</b>", showarrow=False,
            font=dict(size=12, color='black'), xanchor='left', yanchor='top'))
        current_y -= 0.08
        for label in sorted([v for v in set(label2_dict.values()) if pd.notna(v)]):
            color = groupby2_label_color_map.get(label, 'grey') if isinstance(groupby2_label_color_map, dict) else 'grey'
            legend_annotations.append(dict(
                x=legend_start_x, y=current_y, xref='paper', yref='paper',
                text=f"<span style='color:{color}'>■</span> {label}",
                showarrow=False, font=dict(size=12), xanchor='left', yanchor='middle'))
            current_y -= 0.04
    return legend_annotations


# ==========================================================
# Unified Heatmap (categorical annotations)
# ==========================================================

def plot_unified_heatmap(
    adata, genes, groupby1, groupby2=None, labels=None, log=False, z_score=False,
    boundary=False, color_map='Viridis', groupby1_label_color_map=None,
    groupby2_label_color_map=None, max_cells=50000, n_bins=10000, transformation=None, adata_obs=None
):
    log, z_score = _map_transformation_args(transformation, log, z_score)
    if groupby2 and _is_continuous_annotation(adata, groupby2):
        return plot_heatmap2_continuous(
            adata, genes, groupby1, groupby2, labels, log, z_score,
            color_map, groupby1_label_color_map, max_cells, n_bins, adata_obs
        )

    valid_genes = _validate_genes(adata, genes)
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(text='No valid genes found in the dataset', xref='paper', yref='paper', x=0.5, y=0.5, showarrow=False, font=dict(size=14))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    ctx = _filter_cells_and_obs(adata, groupby1, labels)
    original_adata = ctx['original_adata']
    filtered_obs = ctx['filtered_obs']
    filtered_obs_names = ctx['filtered_obs_names']

    gene_df = _extract_gene_df_after_filter(ctx, valid_genes)
    gene_df = _apply_transformations(gene_df, valid_genes, log, z_score)

    annotation_columns = [groupby1]
    if groupby2 and groupby2 != 'None' and groupby2 != groupby1:
        annotation_columns.append(groupby2)
    label_df = pd.DataFrame(filtered_obs[annotation_columns])
    label_df.insert(0, 'CellID', filtered_obs_names)
    heatmap_df = pd.merge(gene_df, label_df, on='CellID')

    # Binning for large datasets (apply even when secondary categorical annotation is present)
    if len(heatmap_df) > max_cells:
        if groupby2 and groupby2 != 'None' and groupby2 != groupby1:
            heatmap_df['__combined__'] = heatmap_df[groupby1].astype(str) + '_' + heatmap_df[groupby2].astype(str)
            binned = bin_cells_for_heatmap(heatmap_df, valid_genes, '__combined__', n_bins, continuous_key=None)
            # Reconstruct groupby1 and groupby2 from combined key
            def _split_combined(val):
                s = str(val)
                if '_' in s:
                    a, b = s.rsplit('_', 1)
                    return a, b
                return s, s
            if len(binned) > 0:
                gb1_vals, gb2_vals = zip(*[_split_combined(v) for v in binned['__combined__'].values])
                binned[groupby1] = list(gb1_vals)
                binned[groupby2] = list(gb2_vals)
            else:
                binned[groupby1] = []
                binned[groupby2] = []
            heatmap_df = binned.drop(columns=['__combined__'])
        else:
            heatmap_df = bin_cells_for_heatmap(heatmap_df, valid_genes, groupby1, n_bins, continuous_key=None)

    if labels:
        heatmap_df[groupby1] = pd.Categorical(heatmap_df[groupby1], categories=labels, ordered=True)
    sort_columns = [groupby1]
    if len(annotation_columns) > 1:
        sort_columns.append(groupby2)
    sorted_heatmap_df = heatmap_df.sort_values(sort_columns)
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T

    label_list1 = sorted_heatmap_df[groupby1].unique().tolist()
    value_list1 = [sorted_heatmap_df[sorted_heatmap_df[groupby1] == item].shape[0] for item in label_list1]

    has_secondary = len(annotation_columns) > 1
    groupby1_label_color_map, groupby2_label_color_map = _default_color_maps(adata_obs, original_adata, groupby1, groupby2 if has_secondary else None, groupby1_label_color_map, groupby2_label_color_map)

    if has_secondary:
        sorted_heatmap_df['combined'] = sorted_heatmap_df[groupby1].astype(str) + '_' + sorted_heatmap_df[groupby2].astype(str)
        label_list2 = sorted_heatmap_df['combined'].unique().tolist()
        label2_dict = {item: item.split('_')[-1] for item in label_list2}
        value_list2 = [sorted_heatmap_df[sorted_heatmap_df['combined'] == item].shape[0] for item in label_list2]
    else:
        label_list2, label2_dict, value_list2 = [], {}, []

    y_position = list(range(len(valid_genes)))
    heatmap_height = 40 * len(valid_genes)
    bar_chart_height1 = 30
    bar_chart_height2 = 30 if has_secondary else 0
    if has_secondary:
        total_height = [heatmap_height, bar_chart_height1, bar_chart_height2]
        rows = 3
    else:
        total_height = [heatmap_height, bar_chart_height1]
        rows = 2
    total_x_range = sum(value_list1)

    fig = make_subplots(rows=rows, cols=1, row_heights=total_height, shared_xaxes=True, vertical_spacing=0.02)

    colorbar_len = 0.4 if has_secondary else 0.5
    fig.add_trace(_make_expression_heatmap(heatmap_gene_matrix, valid_genes, color_map, log, z_score, colorbar_len), row=1, col=1)

    if boundary is not False:
        boundary_width = boundary if isinstance(boundary, (int, float)) else 1
        _add_boundaries(fig, value_list1, row=1, col=1, width=boundary_width, n_genes=len(valid_genes))

    for tr in _make_annotation_bar(label_list1, value_list1, groupby1_label_color_map, width=1):
        fig.add_trace(tr, row=2, col=1)

    if has_secondary:
        for i, label in enumerate(label_list2):
            secondary_label = label2_dict[label]
            color = groupby2_label_color_map.get(secondary_label, 'grey') if isinstance(groupby2_label_color_map, dict) else 'grey'
            fig.add_trace(go.Bar(
                x=[value_list2[i]],
                marker_color=color,
                name=f"{secondary_label}",
                hovertemplate=f'<b>{secondary_label}</b><br>Count: {value_list2[i]}<extra></extra>',
                orientation='h',
                showlegend=False,
                width=2,
            ), row=3, col=1)

    legend_annotations = _legend_annotations(groupby1, label_list1, groupby1_label_color_map, has_secondary, groupby2, label2_dict, groupby2_label_color_map if has_secondary else {})

    layout_updates = {
        'barmode': 'stack',
        'showlegend': False,
        'plot_bgcolor': 'white',
        'paper_bgcolor': 'white',
        'annotations': legend_annotations,
        'hovermode': 'closest',
        'height': max(450, sum(total_height)),
        'margin': dict(t=50, b=150 if not has_secondary else 150, l=50, r=200),
        'xaxis': dict(range=[0, total_x_range], constrain='domain', constraintoward='center'),
        'yaxis': dict(
            tickmode='array', tickvals=y_position, ticktext=valid_genes, showgrid=False, zeroline=False,
            tickfont=dict(size=14), constrain='domain', constraintoward='middle', range=[-0.5, len(valid_genes)-0.5]
        ),
        'yaxis2': dict(
            tickmode='array', tickvals=[0], ticktext=[f"<b>{groupby1}</b>"], showgrid=False, zeroline=False,
            tickfont=dict(size=14), constrain='domain', constraintoward='middle'
        )
    }
    if has_secondary:
        layout_updates['xaxis2'] = dict(range=[0, total_x_range], visible=False, constrain='domain', constraintoward='center')
        layout_updates['xaxis3'] = dict(range=[0, total_x_range], visible=False, constrain='domain', constraintoward='center')
        layout_updates['yaxis3'] = dict(
            tickmode='array', tickvals=[0], ticktext=[f"<b>{groupby2}</b>"], showgrid=False, zeroline=False,
            tickfont=dict(size=14), constrain='domain', constraintoward='middle'
        )
    else:
        layout_updates['xaxis2'] = dict(range=[0, total_x_range], visible=False, constrain='domain', constraintoward='center')

    fig.update_layout(**layout_updates)
    fig.update_xaxes(tickangle=90, visible=False)
    return fig


# ==========================================================
# Heatmap with continuous secondary annotation
# ==========================================================

def plot_heatmap2_continuous(
    adata, genes, groupby1, continuous_key, labels=None, log=False, z_score=False,
    color_map='Viridis', groupby1_label_color_map=None, max_cells=50000, n_bins=10000, adata_obs=None
):
    valid_genes = _validate_genes(adata, genes)
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(text='No valid genes found in the dataset', xref='paper', yref='paper', x=0.5, y=0.5, showarrow=False, font=dict(size=14))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    ctx = _filter_cells_and_obs(adata, groupby1, labels)
    filtered_obs = ctx['filtered_obs']
    filtered_obs_names = ctx['filtered_obs_names']

    if labels and ctx['is_backed']:
        gene_df_list = []
        for gene in valid_genes:
            gene_expr = extract_gene_expression(adata, gene)
            gene_expr_filtered = gene_expr[ctx['cell_indices_array']]
            gene_df_list.append(pd.Series(gene_expr_filtered, name=gene, index=filtered_obs_names))
        gene_df = pd.concat(gene_df_list, axis=1)
    else:
        gene_df = extract_multiple_genes(ctx['adata'], valid_genes)

    gene_df.insert(0, 'CellID', filtered_obs_names)
    label_df = pd.DataFrame(filtered_obs[[groupby1, continuous_key]])
    label_df.insert(0, 'CellID', filtered_obs_names)
    heatmap_df = pd.merge(gene_df, label_df, on='CellID')

    gene_df = _apply_transformations(heatmap_df[['CellID'] + valid_genes].copy(), valid_genes, log, z_score)
    heatmap_df[valid_genes] = gene_df[valid_genes]

    sorted_heatmap_df = heatmap_df.sort_values(continuous_key)
    use_binning = len(sorted_heatmap_df) > max_cells
    if use_binning:
        sorted_heatmap_df = bin_cells_for_heatmap(sorted_heatmap_df, valid_genes, groupby1, n_bins, continuous_key=continuous_key)
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T

    unique_labels = sorted(adata_obs[groupby1].unique()) if adata_obs is not None else sorted(adata.obs[groupby1].unique())
    if groupby1_label_color_map is None:
        from guanaco.data_loader import color_config
        groupby1_label_color_map = {label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels)}

    heatmap_height = 40 * len(valid_genes)
    continuous_bar_height = 30
    bar_chart_height = 30

    fig = make_subplots(rows=3, cols=1, row_heights=[heatmap_height, continuous_bar_height, bar_chart_height], shared_xaxes=True, vertical_spacing=0.01)

    fig.add_trace(_make_expression_heatmap(heatmap_gene_matrix, valid_genes, color_map, log, z_score, colorbar_len=0.4), row=1, col=1)

    fig.add_trace(go.Heatmap(
        z=sorted_heatmap_df[continuous_key].values.reshape(1, -1),
        x=list(range(len(sorted_heatmap_df))),
        y=[continuous_key],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=f'{continuous_key}: %{{z:.4f}}<extra></extra>'
    ), row=2, col=1)

    category_to_num = {cat: i for i, cat in enumerate(sorted_heatmap_df[groupby1].unique())}
    category_values = sorted_heatmap_df[groupby1].map(category_to_num).values.reshape(1, -1)
    unique_categories = sorted(sorted_heatmap_df[groupby1].unique())
    n_categories = len(unique_categories)
    colorscale = []
    for i, cat in enumerate(unique_categories):
        color = groupby1_label_color_map.get(cat, 'grey') if isinstance(groupby1_label_color_map, dict) else 'grey'
        colorscale.append([i/n_categories, color])
        colorscale.append([(i+1)/n_categories, color])
    fig.add_trace(go.Heatmap(
        z=category_values,
        x=list(range(len(sorted_heatmap_df))),
        y=[groupby1],
        colorscale=colorscale,
        showscale=False,
        text=[[cat for cat in sorted_heatmap_df[groupby1]]],
        hovertemplate='%{text}<extra></extra>',
        zmin=-0.5,
        zmax=n_categories-0.5
    ), row=3, col=1)

    legend_annotations = []
    legend_start_x = 1.01
    legend_start_y = 0.5
    legend_annotations.append(dict(
        x=legend_start_x, y=legend_start_y,
        xref='paper', yref='paper', text=f"<b>{groupby1}</b>", showarrow=False,
        font=dict(size=12, color='black'), xanchor='left', yanchor='top'
    ))
    current_y = legend_start_y - 0.08
    for label in unique_categories:
        color = groupby1_label_color_map.get(label, 'grey') if isinstance(groupby1_label_color_map, dict) else 'grey'
        legend_annotations.append(dict(
            x=legend_start_x, y=current_y,
            xref='paper', yref='paper', text=f"<span style='color:{color}'>■</span> {label}",
            showarrow=False, font=dict(size=12), xanchor='left', yanchor='middle'
        ))
        current_y -= 0.04

    fig.update_layout(
        plot_bgcolor='white', paper_bgcolor='white',
        annotations=legend_annotations, hovermode='closest',
        height=max(450, sum([heatmap_height, continuous_bar_height, bar_chart_height])),
        margin=dict(t=50, b=50, l=50, r=150),
        xaxis=dict(visible=False), xaxis2=dict(visible=False), xaxis3=dict(visible=False),
        yaxis=dict(
            tickmode='array', tickvals=list(range(len(valid_genes))), ticktext=valid_genes,
            showgrid=False, zeroline=False, tickfont=dict(size=12), range=[-0.5, len(valid_genes)-0.5]
        ),
        yaxis2=dict(
            tickmode='array', tickvals=[0], ticktext=[continuous_key], showgrid=False, zeroline=False,
            tickfont=dict(size=12)
        ),
        yaxis3=dict(
            tickmode='array', tickvals=[0], ticktext=[groupby1], showgrid=False, zeroline=False,
            tickfont=dict(size=12)
        ),
    )
    return fig
