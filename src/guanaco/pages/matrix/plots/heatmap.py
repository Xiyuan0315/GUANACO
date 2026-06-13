import hashlib
from collections import OrderedDict

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression, apply_transformation, bin_cells_for_heatmap
)

# Standardization is one of: None (raw), "minmax", "zscore".
#   "minmax" -> Scanpy standard_scale: per-gene (x-min)/(max-min) -> [0, 1].
#   "zscore" -> per-gene (x-mean)/std, colour range clipped to ±ZSCORE_COLOR_CLIP.
# z-scores are unbounded -- a single low-variance or very sparse gene can reach ±20
# and wash out the whole map -- so the symmetric colour range is capped at ±2.5 (the
# Seurat/Scanpy vmin/vmax convention). The SAME cap is used in the heatmap and dotplot
# so a given z-score maps to the same colour in both.
ZSCORE_COLOR_CLIP = 2.5


def _canonical_standardization(standardization):
    """Map UI / legacy values onto the canonical None | 'minmax' | 'zscore'."""
    if standardization in (None, "None", "none"):
        return None
    if standardization in ("minmax", "var"):  # 'var' = legacy notebook min-max
        return "minmax"
    if standardization in ("zscore", "z_score", "across_cells", "across_groups"):
        return "zscore"
    return None


# Cache of standardized + binned per-gene vectors (each ~n_bins long, ~16 KB).
# Keyed by (bin-plan id, gene, layer, standardization). The bin plan depends only
# on grouping/labels/cells/layer/standardization — NOT on which genes are shown —
# so when the user adds/removes genes one at a time, every unchanged gene is a hit
# and only the new gene is computed. This keeps the heatmap from rebuilding (and
# holding) the full cells x genes matrix on every edit.
_BINNED_GENE_CACHE = OrderedDict()
_BINNED_GENE_CACHE_MAX = 4000  # ~64 MB ceiling at 4000 bins x float32


def _binned_cache_get(key):
    val = _BINNED_GENE_CACHE.get(key)
    if val is not None:
        _BINNED_GENE_CACHE.move_to_end(key)
    return val


def _binned_cache_put(key, val):
    _BINNED_GENE_CACHE[key] = val
    _BINNED_GENE_CACHE.move_to_end(key)
    while len(_BINNED_GENE_CACHE) > _BINNED_GENE_CACHE_MAX:
        _BINNED_GENE_CACHE.popitem(last=False)

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


def _extract_gene_df_after_filter(ctx, valid_genes, layer=None):
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
        expr = np.asarray(extract_gene_expression(source_adata, gene, layer=layer, use_cache=True), dtype=np.float32)
        gene_matrix[:, j] = expr[row_indices] if row_indices is not None else expr
    gene_df = pd.DataFrame(gene_matrix, columns=valid_genes, index=filtered_obs_names)
    return gene_df


def _apply_transformations(gene_df, valid_genes, log, standardization):
    """Apply optional log1p then per-gene standardization over all cells.

    standardization:
        "minmax" -> each gene min-max scaled to [0, 1] (Scanpy standard_scale='var')
        "zscore" -> each gene z-scored (x-mean)/std
    """
    if log:
        gene_df[valid_genes] = apply_transformation(gene_df[valid_genes], method='log1p', copy=True)
    if standardization == "minmax":
        block = gene_df[valid_genes]
        lo = block.min(axis=0)
        rng = (block.max(axis=0) - lo).replace(0, 1.0)
        gene_df[valid_genes] = (block - lo) / rng
    elif standardization == "zscore":
        block = gene_df[valid_genes]
        mu = block.mean(axis=0)
        sd = block.std(axis=0).replace(0, 1.0)
        gene_df[valid_genes] = (block - mu) / sd
    return gene_df


def _make_expression_heatmap(matrix, genes, color_map, standardization, colorbar_len, x=None):
    color_map = resolve_continuous_colorscale(color_map)
    if standardization == "minmax":
        # Scanpy standard_scale: already in [0, 1], sequential scale.
        zmin, zmax, zmid = 0.0, 1.0, None
        title = "Scaled (0–1)"
    elif standardization == "zscore":
        # Cap the symmetric range so outlier z-scores don't wash out the map.
        zmin, zmax, zmid = -ZSCORE_COLOR_CLIP, ZSCORE_COLOR_CLIP, 0
        title = "Z-score"
    else:
        zmin, zmax, zmid = matrix.min(), matrix.max(), None
        title = "Expression"
    kwargs = {}
    if x is not None:
        # Cell-weighted column centres so each bin's width reflects how many
        # original cells it represents (keeps the annotation bar aligned).
        kwargs["x"] = x
    return go.Heatmap(
        z=matrix,
        colorscale=color_map,
        hoverinfo='skip',
        zmin=zmin,
        zmax=zmax,
        zmid=zmid,
        colorbar=dict(
            title=title,
            len=colorbar_len,
            y=1,
            yanchor='top'
        ),
        **kwargs,
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


def _default_color_maps(adata_obs, original_adata, groupby1, groupby2, groupby1_label_color_map, groupby2_label_color_map, color_config=None):
    if color_config is None:
        from guanaco.data.registry import color_config

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
    for label in [lbl for lbl in label_list1 if pd.notna(lbl)]:
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

def _streaming_primary_matrix(ctx, valid_genes, groupby1, labels, max_cells, n_bins,
                              standardization, layer, filtered_obs, cache_sig=None):
    """Build the (genes x columns) heatmap matrix one gene at a time.

    Computes the bin assignment once (from group labels only), then for each gene
    extracts its column, standardizes, and reduces it to per-bin means via a cached
    lookup. Peak memory is one full gene column plus the binned output, instead of
    a dense cells x genes matrix. Returns the same column ordering (by group) the
    full-matrix path produces. Only the primary-categorical, non-log case.
    """
    adata_src = ctx['original_adata'] if ctx['cell_indices_array'] is not None else ctx['adata']
    row_indices = ctx['cell_indices_array']
    row_indices = None if row_indices is None else np.asarray(row_indices, dtype=np.int64)

    group_series = filtered_obs[groupby1]
    groups_arr = group_series.to_numpy()
    n_filtered = len(groups_arr)

    present = set(pd.unique(groups_arr).tolist())
    if labels:
        order = [g for g in labels if g in present]
    else:
        order = sorted(present, key=str)

    # --- bin assignment (matches bin_cells_for_heatmap group-wise semantics) ---
    binning = n_filtered > max_cells
    row_to_bin = np.empty(n_filtered, dtype=np.int64)
    bin_sizes = []
    bin_group = []
    nb = 0
    for g in order:
        pos = np.flatnonzero(groups_arr == g)
        m = pos.size
        if m == 0:
            continue
        if (not binning) or m <= 10:
            row_to_bin[pos] = nb + np.arange(m)
            bin_sizes.extend([1] * m)
            bin_group.extend([g] * m)
            nb += m
        else:
            gbins = min(max(1, int(n_bins * m / n_filtered)), m)
            edges = np.unique(np.linspace(0, m, num=gbins + 1, dtype=int))
            starts, ends = edges[:-1], edges[1:]
            for b, (s, e) in enumerate(zip(starts, ends)):
                row_to_bin[pos[s:e]] = nb + b
                bin_sizes.append(int(e - s))
                bin_group.append(g)
            nb += len(starts)
    n_out = nb
    bin_sizes = np.asarray(bin_sizes, dtype=np.float64)
    bin_group = np.asarray(bin_group, dtype=object)

    # Stable id for the bin plan: unchanged across gene add/remove (same adata,
    # grouping, labels, cell set), so per-gene binned vectors stay cacheable.
    # Prefer a caller-supplied content signature -- it stays stable for the same
    # cell selection, whereas a fresh adata[selected_cells] view changes id() on
    # every render and would defeat the cache during an active lasso/filter.
    if cache_sig is not None:
        plan_id = str(cache_sig)
    else:
        backed = bool(getattr(adata_src, "isbacked", False) and getattr(adata_src, "filename", None))
        src_id = ("backed", str(adata_src.filename)) if backed else ("mem", id(adata_src))
        plan_id = hashlib.md5(
            repr((src_id, tuple(adata_src.shape), groupby1,
                  tuple(labels) if labels else None, int(max_cells), int(n_bins))).encode()
        ).hexdigest()

    matrix = np.empty((len(valid_genes), n_out), dtype=np.float32)
    for j, gene in enumerate(valid_genes):
        key = (plan_id, gene, layer, standardization)
        cached = _binned_cache_get(key)
        if cached is None or cached.shape[0] != n_out:
            col = np.asarray(
                extract_gene_expression(adata_src, gene, layer=layer, use_cache=True),
                dtype=np.float32,
            )
            if row_indices is not None:
                col = col[row_indices]
            if standardization == "minmax":
                lo = col.min()
                rng = col.max() - lo
                col = (col - lo) / (rng if rng else 1.0)
            elif standardization == "zscore":
                # Single-pass: compute mean and std together to avoid scanning twice.
                mu = col.mean()
                sd = np.sqrt(((col - mu) ** 2).mean())
                col = (col - mu) / (sd if sd else 1.0)
            sums = np.bincount(row_to_bin, weights=col, minlength=n_out)
            cached = (sums / bin_sizes).astype(np.float32)
            _binned_cache_put(key, cached)
        matrix[j, :] = cached

    label_list1 = list(order)
    value_list1 = [int(bin_sizes[bin_group == g].sum()) for g in order]
    return matrix, label_list1, value_list1, bin_sizes, binning


def plot_unified_heatmap(
    adata, genes, groupby1, groupby2=None, labels=None, log=False, z_score=False,
    boundary=False, color_map='Viridis', groupby1_label_color_map=None,
    groupby2_label_color_map=None, max_cells=10000, n_bins=4000, transformation=None,
    standardization=None, layer=None, cache_sig=None,
    adata_obs=None, data_already_filtered=False, color_config=None
):
    log, z_score = _map_transformation_args(transformation, log, z_score)
    # Legacy z_score (notebook API) maps onto z-score standardization.
    if standardization is None and z_score:
        standardization = "zscore"
    standardization = _canonical_standardization(standardization)
    if groupby2 and _is_continuous_annotation(adata, groupby2):
        return plot_heatmap2_continuous(
            adata, genes, groupby1, groupby2, labels, log, standardization,
            color_map, groupby1_label_color_map, max_cells, n_bins, adata_obs,
            data_already_filtered=data_already_filtered,
            color_config=color_config, layer=layer,
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

    has_secondary = bool(groupby2 and groupby2 != 'None' and groupby2 != groupby1)
    # Streaming, per-gene-cached path for the common case (single categorical
    # annotation, no log). It avoids materializing the full cells x genes matrix
    # and makes add/remove-a-gene reuse every unchanged gene. The secondary-
    # annotation and log paths keep the original full-matrix implementation.
    use_streaming = (not has_secondary) and (not log)

    if use_streaming:
        heatmap_gene_matrix, label_list1, value_list1, bin_sizes, is_binned = _streaming_primary_matrix(
            ctx, valid_genes, groupby1, labels, max_cells, n_bins, standardization, layer, filtered_obs,
            cache_sig=cache_sig,
        )
        groupby1_label_color_map, _ = _default_color_maps(
            adata_obs, original_adata, groupby1, None,
            groupby1_label_color_map, None, color_config=color_config,
        )
        groupby2_label_color_map = None
        label_list2, label2_dict, value_list2 = [], {}, []
    else:
        gene_df = _extract_gene_df_after_filter(ctx, valid_genes, layer=layer)
        gene_df = _apply_transformations(gene_df, valid_genes, log, standardization)

        annotation_columns = [groupby1]
        if has_secondary:
            annotation_columns.append(groupby2)
        heatmap_df = gene_df
        for annotation in annotation_columns:
            heatmap_df[annotation] = filtered_obs[annotation].to_numpy()

        # Binning for large datasets (apply even when secondary categorical annotation is present)
        if len(heatmap_df) > max_cells:
            effective_bins = _resolve_effective_bins(len(heatmap_df), n_bins)
            if has_secondary:
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

        # Per-column weight = number of original cells it represents. When binning
        # was applied each column is a bin (>=1 cell); otherwise every column is one
        # cell. Annotation-bar widths and group boundaries use these true cell counts
        # rather than column counts, and the heatmap x is weighted to match.
        is_binned = '__bin_size__' in sorted_heatmap_df.columns
        if is_binned:
            bin_sizes = sorted_heatmap_df['__bin_size__'].to_numpy(dtype=np.float64)
        else:
            bin_sizes = np.ones(len(sorted_heatmap_df), dtype=np.float64)

        label_list1 = sorted_heatmap_df[groupby1].unique().tolist()
        gb1_arr = sorted_heatmap_df[groupby1].to_numpy()
        value_list1 = [int(bin_sizes[gb1_arr == lbl].sum()) for lbl in label_list1]

        groupby1_label_color_map, groupby2_label_color_map = _default_color_maps(adata_obs, original_adata, groupby1, groupby2 if has_secondary else None, groupby1_label_color_map, groupby2_label_color_map, color_config=color_config)

        if has_secondary:
            sorted_heatmap_df['combined'] = list(zip(sorted_heatmap_df[groupby1], sorted_heatmap_df[groupby2]))
            label_list2 = sorted_heatmap_df['combined'].unique().tolist()
            label2_dict = {item: item[1] if isinstance(item, tuple) and len(item) >= 2 else item for item in label_list2}
            combined_vals = sorted_heatmap_df['combined'].tolist()
            value_list2 = [
                int(sum(bs for bs, c in zip(bin_sizes, combined_vals) if c == lbl))
                for lbl in label_list2
            ]
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
    if is_binned:
        # Place each column at its cumulative cell-count centre so column widths
        # (and hence each group's horizontal extent) reflect true cell numbers.
        edges = np.concatenate([[0.0], np.cumsum(bin_sizes)])
        x_centers = (edges[:-1] + edges[1:]) / 2.0
    else:
        x_centers = None
    fig.add_trace(_make_expression_heatmap(heatmap_gene_matrix, valid_genes, color_map, standardization, colorbar_len, x=x_centers), row=1, col=1)

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
    adata, genes, groupby1, continuous_key, labels=None, log=False, standardization=None,
    color_map='Viridis', groupby1_label_color_map=None, max_cells=10000, n_bins=4000,
    adata_obs=None, data_already_filtered=False, color_config=None, layer=None
):
    standardization = _canonical_standardization(standardization)
    valid_genes = _validate_genes(adata, genes)
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(text='No valid genes found in the dataset', xref='paper', yref='paper', x=0.5, y=0.5, showarrow=False, font=dict(size=14))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    ctx = _filter_cells_and_obs(adata, groupby1, labels, data_already_filtered=data_already_filtered)
    filtered_obs = ctx['filtered_obs']

    gene_df = _extract_gene_df_after_filter(ctx, valid_genes, layer=layer)
    gene_df = _apply_transformations(gene_df, valid_genes, log, standardization)
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
        if color_config is None:
            from guanaco.data.registry import color_config
        groupby1_label_color_map = {label: color_config[i % len(color_config)] for i, label in enumerate(unique_labels)}

    heatmap_height = 40 * len(valid_genes)
    continuous_bar_height = 30
    bar_chart_height = 30

    fig = make_subplots(rows=3, cols=1, row_heights=[heatmap_height, continuous_bar_height, bar_chart_height], shared_xaxes=True, vertical_spacing=0.01)

    fig.add_trace(_make_expression_heatmap(heatmap_gene_matrix, valid_genes, color_map, standardization, colorbar_len=0.4), row=1, col=1)

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
