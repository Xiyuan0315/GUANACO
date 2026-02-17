import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression, apply_transformation, bin_cells_for_heatmap
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


def _resolve_effective_bins(n_rows, requested_bins, max_bins=4000):
    if n_rows <= 0:
        return 0
    return max(1, min(int(requested_bins), int(max_bins), int(n_rows)))


def _is_continuous_annotation(adata, annotation, threshold=50):
    if annotation not in adata.obs.columns:
        return False
    dtype = str(adata.obs[annotation].dtype)
    if any(t in dtype for t in ['float', 'int']):
        return adata.obs[annotation].nunique() >= threshold
    return False


def _filter_cells_and_obs(adata, groupby1, labels, data_already_filtered=False):
    is_backed = hasattr(adata, 'isbacked') and adata.isbacked
    original_adata = adata
    filtered_obs = adata.obs
    filtered_obs_names = adata.obs_names
    cell_indices_array = None
    if labels and not data_already_filtered:
        mask = adata.obs[groupby1].isin(labels)
        cell_indices_array = np.where(mask.to_numpy())[0]
        filtered_obs = adata.obs.iloc[cell_indices_array]
        filtered_obs_names = adata.obs_names[cell_indices_array]
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
    if valid_genes is None or len(valid_genes) == 0:
        return pd.DataFrame()
    source_adata = original_adata if cell_indices_array is not None else adata
    row_indices = np.asarray(cell_indices_array, dtype=np.int64) if cell_indices_array is not None else None
    n_cells = len(filtered_obs_names)
    gene_matrix = np.empty((n_cells, len(valid_genes)), dtype=np.float32)
    for j, gene in enumerate(valid_genes):
        expr = np.asarray(extract_gene_expression(source_adata, gene, use_cache=True), dtype=np.float32)
        gene_matrix[:, j] = expr[row_indices] if row_indices is not None else expr
    gene_df = pd.DataFrame(gene_matrix, columns=valid_genes, index=filtered_obs_names)
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
    color_map = _resolve_continuous_colorscale(color_map)
    if z_score:
        z_max = max(abs(matrix.min()), abs(matrix.max()))
        zmin, zmax, zmid = -z_max, z_max, 0
    else:
        zmin, zmax, zmid = matrix.min(), matrix.max(), None
    return go.Heatmap(
        z=matrix,
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


def _resolve_continuous_colorscale(color_map):
    if not isinstance(color_map, str):
        return color_map
    if not color_map.startswith("cc:"):
        return color_map

    cc_name = color_map.split(":", 1)[1].strip()
    if not cc_name:
        return "Viridis"

    try:
        import colorcet as cc
    except Exception:
        return "Viridis"

    palette = cc.palette.get(cc_name)
    if palette is None:
        return "Viridis"

    colors = list(palette)
    if len(colors) < 2:
        return "Viridis"

    denom = len(colors) - 1
    return [[i / denom, c] for i, c in enumerate(colors)]


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

    default_color = color_config
    try:
        from matplotlib import colormaps
        from matplotlib.colors import to_hex

        cmap = colormaps.get_cmap("tab20")
        tab20_colors = [to_hex(cmap(i), keep_alpha=False) for i in range(cmap.N)]
    except Exception:
        tab20_colors = default_color
    if groupby1_label_color_map is None:
        unique_labels_primary = sorted(adata_obs[groupby1].unique()) if adata_obs is not None else sorted(original_adata.obs[groupby1].unique())
        groupby1_label_color_map = {label: default_color[i % len(default_color)] for i, label in enumerate(unique_labels_primary)}
    if groupby2 and groupby2_label_color_map is None:
        unique_labels_secondary = sorted(adata_obs[groupby2].unique()) if adata_obs is not None else sorted(original_adata.obs[groupby2].unique())
        groupby2_label_color_map = {label: tab20_colors[i % len(tab20_colors)] for i, label in enumerate(unique_labels_secondary)}
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
    groupby2_label_color_map=None, max_cells=10000, n_bins=4000, transformation=None,
    adata_obs=None, data_already_filtered=False
):
    log, z_score = _map_transformation_args(transformation, log, z_score)
    if groupby2 and _is_continuous_annotation(adata, groupby2):
        return plot_heatmap2_continuous(
            adata, genes, groupby1, groupby2, labels, log, z_score,
            color_map, groupby1_label_color_map, max_cells, n_bins, adata_obs,
            data_already_filtered=data_already_filtered,
        )

    valid_genes = _validate_genes(adata, genes)
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(text='No valid genes found in the dataset', xref='paper', yref='paper', x=0.5, y=0.5, showarrow=False, font=dict(size=14))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    ctx = _filter_cells_and_obs(adata, groupby1, labels, data_already_filtered=data_already_filtered)
    original_adata = ctx['original_adata']
    filtered_obs = ctx['filtered_obs']
    filtered_obs_names = ctx['filtered_obs_names']

    gene_df = _extract_gene_df_after_filter(ctx, valid_genes)
    gene_df = _apply_transformations(gene_df, valid_genes, log, z_score)

    annotation_columns = [groupby1]
    if groupby2 and groupby2 != 'None' and groupby2 != groupby1:
        annotation_columns.append(groupby2)
    heatmap_df = gene_df.copy()
    for annotation in annotation_columns:
        heatmap_df[annotation] = filtered_obs[annotation].to_numpy()

    # Binning for large datasets (apply even when secondary categorical annotation is present)
    if len(heatmap_df) > max_cells:
        effective_bins = _resolve_effective_bins(len(heatmap_df), n_bins)
        if groupby2 and groupby2 != 'None' and groupby2 != groupby1:
            pair_key = '__pair_key__'
            heatmap_df[pair_key] = list(zip(heatmap_df[groupby1], heatmap_df[groupby2]))
            binned = bin_cells_for_heatmap(heatmap_df, valid_genes, pair_key, effective_bins, continuous_key=None)

            def _split_pair(val):
                if isinstance(val, tuple) and len(val) >= 2:
                    return val[0], val[1]
                return val, val

            if len(binned) > 0:
                gb1_vals, gb2_vals = zip(*[_split_pair(v) for v in binned[pair_key].values])
                binned[groupby1] = list(gb1_vals)
                binned[groupby2] = list(gb2_vals)
            else:
                binned[groupby1] = []
                binned[groupby2] = []
            heatmap_df = binned.drop(columns=[pair_key])
        else:
            heatmap_df = bin_cells_for_heatmap(heatmap_df, valid_genes, groupby1, effective_bins, continuous_key=None)

    if labels:
        heatmap_df[groupby1] = pd.Categorical(heatmap_df[groupby1], categories=labels, ordered=True)
    sort_columns = [groupby1]
    if len(annotation_columns) > 1:
        sort_columns.append(groupby2)
    sorted_heatmap_df = heatmap_df.sort_values(sort_columns)
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].to_numpy(dtype=np.float32, copy=False).T

    label_list1 = sorted_heatmap_df[groupby1].unique().tolist()
    value_list1 = (
        sorted_heatmap_df[groupby1]
        .value_counts(sort=False)
        .reindex(label_list1, fill_value=0)
        .tolist()
    )

    has_secondary = len(annotation_columns) > 1
    groupby1_label_color_map, groupby2_label_color_map = _default_color_maps(adata_obs, original_adata, groupby1, groupby2 if has_secondary else None, groupby1_label_color_map, groupby2_label_color_map)

    if has_secondary:
        sorted_heatmap_df['combined'] = list(zip(sorted_heatmap_df[groupby1], sorted_heatmap_df[groupby2]))
        label_list2 = sorted_heatmap_df['combined'].unique().tolist()
        label2_dict = {item: item[1] if isinstance(item, tuple) and len(item) >= 2 else item for item in label_list2}
        value_list2 = (
            sorted_heatmap_df['combined']
            .value_counts(sort=False)
            .reindex(label_list2, fill_value=0)
            .tolist()
        )
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
    color_map='Viridis', groupby1_label_color_map=None, max_cells=10000, n_bins=4000,
    adata_obs=None, data_already_filtered=False
):
    valid_genes = _validate_genes(adata, genes)
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(text='No valid genes found in the dataset', xref='paper', yref='paper', x=0.5, y=0.5, showarrow=False, font=dict(size=14))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    ctx = _filter_cells_and_obs(adata, groupby1, labels, data_already_filtered=data_already_filtered)
    filtered_obs = ctx['filtered_obs']
    filtered_obs_names = ctx['filtered_obs_names']

    gene_df = _extract_gene_df_after_filter(ctx, valid_genes)
    gene_df = _apply_transformations(gene_df, valid_genes, log, z_score)
    heatmap_df = gene_df.copy()
    heatmap_df[groupby1] = filtered_obs[groupby1].to_numpy()
    heatmap_df[continuous_key] = filtered_obs[continuous_key].to_numpy()

    sorted_heatmap_df = heatmap_df.sort_values(continuous_key)
    use_binning = len(sorted_heatmap_df) > max_cells
    if use_binning:
        effective_bins = _resolve_effective_bins(len(sorted_heatmap_df), n_bins)
        sorted_heatmap_df = bin_cells_for_heatmap(
            sorted_heatmap_df,
            valid_genes,
            groupby1,
            effective_bins,
            continuous_key=continuous_key,
        )
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].to_numpy(dtype=np.float32, copy=False).T

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
        y=[continuous_key],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=f'{continuous_key}: %{{z:.4f}}<extra></extra>'
    ), row=2, col=1)
    
    unique_categories = sorted(sorted_heatmap_df[groupby1].unique())

    category_to_num = {cat: i for i, cat in enumerate(unique_categories)}
    category_values = sorted_heatmap_df[groupby1].map(category_to_num).values.reshape(1, -1)

    n_categories = len(unique_categories)
    colorscale = []
    for i, cat in enumerate(unique_categories):
        color = groupby1_label_color_map.get(cat, 'grey') if isinstance(groupby1_label_color_map, dict) else 'grey'
        colorscale.append([i / n_categories, color])
        colorscale.append([(i + 1) / n_categories, color])
    fig.add_trace(go.Heatmap(
        z=category_values,
        y=[groupby1],
        colorscale=colorscale,
        showscale=False,
        hoverinfo='skip',
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
