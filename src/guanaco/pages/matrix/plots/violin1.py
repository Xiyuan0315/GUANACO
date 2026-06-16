from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression, apply_transformation, prewarm_gene_cache
)
import hashlib
import time

# Global cache for violin plot data
_violin_data_cache = {}

# A violin is a KDE of a distribution, so a uniform random subsample reproduces the
# same shape while keeping the figure payload and client-side KDE cost bounded.
# Points are never drawn individually here (points=False); the sample only feeds the
# KDE. We bound the points emitted *per gene* (summed across that gene's group
# violins), distributed across groups in proportion to size with a per-group floor.
# This caps the payload regardless of cell count -- the dominant cost of switching to
# the violin tab -- so a 10k- and a 1M-cell dataset ship the same small figure.
MAX_POINTS_PER_GENE = 3_000   # total KDE sample per subplot (one row per gene)
MIN_POINTS_PER_GROUP = 50     # floor so small groups still form a visible violin
MAX_VIOLIN_POINTS = 5_000     # hard per-group ceiling

default_color = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"]

# Neutral statistical overlay (box / IQR), so it reads as a summary annotation
# rather than competing with the category-colored violin fill.
OVERLAY_LINE_COLOR = "#2C3E50"
OVERLAY_FILL_COLOR = "rgba(255,255,255,0.5)"
# Alpha applied to the category color of the violin fill.
VIOLIN_FILL_ALPHA = 0.6


def _to_rgba(color, alpha):
    """Return ``color`` as an ``rgba(...)`` string with ``alpha`` (hex/rgb in)."""
    if isinstance(color, str):
        c = color.strip()
        if c.startswith("#") and len(c) == 7:
            r, g, b = int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)
            return f"rgba({r},{g},{b},{alpha})"
        if c.lower().startswith("rgb(") and c.endswith(")"):
            parts = [p.strip() for p in c[c.find("(") + 1:-1].split(",")[:3]]
            if len(parts) == 3:
                return f"rgba({parts[0]},{parts[1]},{parts[2]},{alpha})"
    return color


def _build_group_sample(label_masks, labels, rng):
    """Row positions to sample for each group's violin, shared across all genes.

    Bounds the total points per gene to ``MAX_POINTS_PER_GENE``, allocated across
    groups in proportion to size (with a per-group floor, capped at the group size
    and the hard per-group ceiling). One plan is reused for every gene so all
    subplots sample a consistent set of cells. Returns ``label -> row positions``.
    """
    positions = {lab: np.flatnonzero(label_masks[lab]) for lab in labels if lab in label_masks}
    total = sum(pos.size for pos in positions.values())
    if total <= MAX_POINTS_PER_GENE:
        return positions
    plan = {}
    for lab, pos in positions.items():
        n = pos.size
        alloc = int(round(MAX_POINTS_PER_GENE * n / total))
        alloc = min(max(alloc, MIN_POINTS_PER_GROUP), n, MAX_VIOLIN_POINTS)
        plan[lab] = np.sort(rng.choice(pos, alloc, replace=False)) if alloc < n else pos
    return plan


def _label_masks(obs_values, labels):
    obs_values_array = obs_values.to_numpy() if hasattr(obs_values, 'to_numpy') else np.array(obs_values)
    masks = {}
    for label in labels:
        mask = obs_values_array == label
        if np.any(mask):
            masks[label] = mask
    return masks


def _row_positions_for_labels(adata, groupby, labels, data_already_filtered):
    if labels and not data_already_filtered:
        return np.where(adata.obs[groupby].isin(labels).to_numpy())[0]
    return None


def _extract_gene_frame(adata, valid_genes, row_pos, layer=None):
    prewarm_gene_cache(adata, valid_genes, layer=layer)
    gene_data = {}
    for gene in valid_genes:
        col = extract_gene_expression(adata, gene, layer=layer)
        gene_data[gene] = col[row_pos] if row_pos is not None else col
    return pd.DataFrame(gene_data)


def _apply_gene_transformations(gene_df, valid_genes, transformation):
    if not transformation:
        return gene_df
    for gene in valid_genes:
        gene_df[gene] = apply_transformation(gene_df[gene], method=transformation)
    return gene_df


def _cache_violin_data(cache_key, cached_data):
    if len(_violin_data_cache) >= 50:
        oldest_key = next(iter(_violin_data_cache))
        del _violin_data_cache[oldest_key]
    _violin_data_cache[cache_key] = cached_data


def _create_cache_key(genes, labels, groupby, transformation, adata_id):
    """Create a unique cache key for violin data."""
    key_data = {
        # Preserve user-specified order to avoid incorrect cache hits.
        'genes': tuple(genes),
        'labels': tuple(labels) if labels else None,
        'groupby': groupby,
        'transformation': transformation,
        'adata_id': adata_id,
    }
    return hashlib.md5(str(key_data).encode()).hexdigest()


def _get_adata_id(adata):
    """Get a unique identifier for the adata object."""
    if hasattr(adata, 'isbacked') and adata.isbacked and hasattr(adata, 'filename'):
        return adata.filename
    return f"{id(adata)}_{adata.shape}"


def _label_color_map(adata, groupby, groupby_label_color_map, adata_obs=None):
    unique_labels = (
        sorted(adata_obs[groupby].unique())
        if adata_obs is not None
        else sorted(adata.obs[groupby].unique())
    )
    if groupby_label_color_map is not None:
        return groupby_label_color_map
    return {label: default_color[i % len(default_color)] for i, label in enumerate(unique_labels)}


def _subplot_geometry(valid_genes):
    num_genes = len(valid_genes)
    fig_height = max(400, num_genes * 140)
    # Keep a consistent ~28px gap between rows so adjacent y-axis tick labels (the
    # 0.0 of one row and the top tick of the next) don't overlap, for any number
    # of genes. vertical_spacing is a fraction of height, capped to Plotly's max.
    target_gap_px = 28
    vertical_spacing = min(target_gap_px / fig_height, 0.9 / max(num_genes - 1, 1))
    return num_genes, fig_height, vertical_spacing


def _violin_spanmode(group_expr):
    return 'soft' if np.var(group_expr) < 1e-10 else 'hard'


def _expression_range(gene_values):
    expr_min = gene_values.min()
    expr_max = gene_values.max()
    if expr_max == expr_min:
        pad = 1e-3 if expr_max == 0 else abs(expr_max) * 0.05
        return [expr_min - pad, expr_max + pad]
    return [expr_min, expr_max]


def _add_gene_row_axes(fig, gene, row, expr_range, labels):
    """Per-row axes for one gene: the gene name as the left y tick, numbers on the right.

    The two labels sit on opposite sides (left/right) so they're each at a natural
    axis edge and can't overlap. The gene name is a real y-axis *tick label* with
    ``automargin`` -- exactly how the heatmap shows row names: Plotly measures the
    actual rendered text (in the real page font) and grows the left margin to fit
    it, so long genome loci like "GL000194.1:100880-101887" are never clipped,
    regardless of font or DPI. The numeric expression scale moves to a right-hand
    secondary axis -- it's the only magnitude cue, since ``scalemode='width'``
    normalizes each violin's width.
    """
    mid = (expr_range[0] + expr_range[1]) / 2
    fig.update_yaxes(
        range=expr_range,
        tickmode="array", tickvals=[mid], ticktext=[str(gene)],
        tickfont=dict(size=12), automargin=True,
        showgrid=False, zeroline=False, showline=True, linewidth=2, linecolor="black",
        row=row, col=1, secondary_y=False,
    )
    # A secondary axis with no trace assigned isn't drawn, so anchor it with an
    # invisible point to make the right-hand numeric scale render.
    fig.add_trace(
        go.Scatter(
            x=[labels[0]], y=[expr_range[0]], mode="markers",
            marker=dict(opacity=0), showlegend=False, hoverinfo="skip",
        ),
        row=row, col=1, secondary_y=True,
    )
    fig.update_yaxes(
        range=expr_range, tickformat=".1f", tickfont=dict(size=10),
        showgrid=False, zeroline=False, showline=True, linewidth=2, linecolor="black",
        row=row, col=1, secondary_y=True,
    )


def _add_violin_trace(fig, group_expr, label, row, show_box, color):
    fig.add_trace(
        go.Violin(
            y=group_expr,
            x=[label] * len(group_expr),
            width=0.8,
            # Show Box Plot -> the violin's inner box (Q1-Q3 + median)
            # styled as a neutral, narrow, light overlay (Median + IQR).
            box=dict(
                visible=show_box,
                width=0.12,
                fillcolor=OVERLAY_FILL_COLOR,
                line=dict(color=OVERLAY_LINE_COLOR, width=1),
            ),
            points=False,
            meanline_visible=True,
            name=label,
            spanmode=_violin_spanmode(group_expr),
            fillcolor=_to_rgba(color, VIOLIN_FILL_ALPHA),
            line_color='DarkSlateGrey',
            hoveron='violins',
            hoverinfo='y',
            jitter=0.05,
            scalemode='width',
            scalegroup=label,
            marker=dict(size=2),
        ),
        row=row,
        col=1,
    )


def _extract_and_cache_violin_data(
    adata,
    genes,
    labels,
    groupby,
    transformation,
    data_already_filtered=False,
    layer=None,
):
    """Extract and cache violin plot data for reuse."""
    adata_id = _get_adata_id(adata)
    cache_key = _create_cache_key(
        genes,
        labels,
        groupby,
        transformation,
        f"{adata_id}_prefilter={data_already_filtered}_layer={layer}",
    )
    
    # Check if data is already cached
    if cache_key in _violin_data_cache:
        return _violin_data_cache[cache_key]
    
    # Extract data
    start_time = time.time()
    
    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        return None

    row_pos = _row_positions_for_labels(adata, groupby, labels, data_already_filtered)
    if adata.n_obs == 0 or (row_pos is not None and row_pos.size == 0):
        return None

    # Per-gene extraction through the shared cache-backed single-column path: each
    # gene's full column is read once, cached, and reused across plots; we then keep
    # the selected rows. Consumers index gene_df / obs_values positionally, so it is
    # row order -- not the index -- that must line up, and both come from row_pos.
    gene_df = _extract_gene_frame(adata, valid_genes, row_pos, layer=layer)
    obs_values = adata.obs.iloc[row_pos][groupby] if row_pos is not None else adata.obs[groupby]
    
    gene_df = _apply_gene_transformations(gene_df, valid_genes, transformation)
    
    # Create cached data structure
    cached_data = {
        'gene_df': gene_df,
        'obs_values': obs_values,
        'valid_genes': valid_genes,
        'extraction_time': time.time() - start_time
    }
    
    _cache_violin_data(cache_key, cached_data)
    return cached_data

def plot_violin1(
    adata,
    genes,
    groupby,
    labels=None,
    transformation=None,
    layer=None,
    show_box=False,
    groupby_label_color_map=None,
    adata_obs=None,
    data_already_filtered=False,
):

    # if no label or no genes, stop updating.
    if len(genes) == 0 or labels is None or len(labels) == 0:
        raise PreventUpdate
    # Extract data
    cached_data = _extract_and_cache_violin_data(
        adata, genes, labels, groupby, transformation, data_already_filtered=data_already_filtered, layer=layer
    )
    if cached_data is None:
        return go.Figure()
    
    gene_df = cached_data['gene_df']
    obs_values = cached_data['obs_values']
    valid_genes = cached_data['valid_genes']
    
    groupby_label_color_map = _label_color_map(adata, groupby, groupby_label_color_map, adata_obs=adata_obs)

    # Create subplots. secondary_y on every row gives each gene a left axis (the
    # gene-name tick) and a right axis (the numeric scale); see _add_gene_row_axes.
    num_genes, fig_height, vertical_spacing = _subplot_geometry(valid_genes)
    fig = make_subplots(
        rows=num_genes, cols=1, shared_xaxes=True, vertical_spacing=vertical_spacing,
        specs=[[{"secondary_y": True}] for _ in range(num_genes)],
    )

    # Pre-calculate boolean masks for each label to avoid repeated groupby operations
    # Use numpy array for faster boolean indexing
    label_masks = _label_masks(obs_values, labels)

    # Fixed seed so each violin's subsample (and thus its shape) is stable across
    # re-renders. The plan is built once and reused for every gene.
    rng = np.random.default_rng(0)
    group_sample = _build_group_sample(label_masks, labels, rng)

    # Create violin plots
    for i, gene in enumerate(valid_genes):
        # Access gene values directly as numpy array (full data drives the y-range).
        gene_values = gene_df[gene].values

        for label in labels:
            if label in group_sample:
                # Payload-bounded KDE sample for this group (see _build_group_sample).
                group_expr = gene_values[group_sample[label]]

                _add_violin_trace(
                    fig,
                    group_expr,
                    label,
                    i + 1,
                    show_box,
                    groupby_label_color_map[label],
                )

        _add_gene_row_axes(fig, gene, i + 1, _expression_range(gene_values), labels)
    # Layout
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=10),
        showlegend=False,
        height=fig_height,
        # Small base margins; automargin on the gene-name (left) axes grows the
        # left side to fit the longest name on its own.
        margin=dict(l=10, r=20, t=20, b=40),
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')

    return fig
