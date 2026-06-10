from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from guanaco.utils.gene_extraction_utils import (
    extract_gene_expression, extract_multiple_genes, apply_transformation
)
import hashlib
import time

# Global cache for violin plot data
_violin_data_cache = {}

# Cap the points sent per violin. A violin is a KDE of the distribution, so a
# uniform random sample of this many points reproduces the same shape while
# keeping the figure payload and client-side KDE cost bounded. Only groups larger
# than this are subsampled (small/rare groups are left intact).
MAX_VIOLIN_POINTS = 5_000

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
    else:
        return f"{id(adata)}_{adata.shape}"

def _extract_and_cache_violin_data(adata, genes, labels, groupby, transformation, data_already_filtered=False, layer=None):
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
    
    # Filter cells based on selected labels, unless caller already passed filtered data.
    if labels and not data_already_filtered:
        cell_indices = adata.obs[groupby].isin(labels)
        if hasattr(adata, 'isbacked') and adata.isbacked:
            cell_indices_array = np.where(cell_indices)[0]
        else:
            cell_indices_array = None
    else:
        cell_indices = None
        cell_indices_array = None
    
    # Extract multiple genes at once for efficiency
    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        return None
    
    # Use batch extraction for better performance
    if len(valid_genes) > 1:
        try:
            if labels and not (hasattr(adata, 'isbacked') and adata.isbacked) and not data_already_filtered:
                filtered_adata = adata[cell_indices]
                gene_df = extract_multiple_genes(filtered_adata, valid_genes, layer=layer)
                obs_values = filtered_adata.obs[groupby]
            elif hasattr(adata, 'isbacked') and adata.isbacked and cell_indices_array is not None:
                # Backed mode + filtered labels: avoid extracting all cells.
                gene_data = {}
                for gene in valid_genes:
                    gene_expr = extract_gene_expression(adata, gene, layer=layer)
                    gene_data[gene] = gene_expr[cell_indices_array]
                gene_df = pd.DataFrame(gene_data)
                obs_values = adata.obs.iloc[cell_indices_array][groupby]
            else:
                gene_df = extract_multiple_genes(adata, valid_genes, layer=layer)
                if cell_indices_array is not None:
                    gene_df = gene_df.iloc[cell_indices_array]
                    obs_values = adata.obs.iloc[cell_indices_array][groupby]
                else:
                    obs_values = adata.obs[groupby]
        except Exception:
            # Fallback to individual extraction
            gene_data = {}
            for gene in valid_genes:
                gene_expr = extract_gene_expression(adata, gene, layer=layer)
                if cell_indices_array is not None:
                    gene_expr = gene_expr[cell_indices_array]
                elif cell_indices is not None:
                    gene_expr = gene_expr[cell_indices]
                gene_data[gene] = gene_expr
            
            gene_df = pd.DataFrame(gene_data)
            if cell_indices is not None:
                if hasattr(adata, 'isbacked') and adata.isbacked:
                    obs_values = adata.obs.iloc[cell_indices_array][groupby]
                else:
                    obs_values = adata.obs[cell_indices][groupby]
            else:
                obs_values = adata.obs[groupby]
    else:
        # Single gene extraction
        gene = valid_genes[0]
        gene_expr = extract_gene_expression(adata, gene, layer=layer)
        if cell_indices_array is not None:
            gene_expr = gene_expr[cell_indices_array]
        elif cell_indices is not None:
            gene_expr = gene_expr[cell_indices]
        
        gene_df = pd.DataFrame({gene: gene_expr})
        if cell_indices is not None:
            if hasattr(adata, 'isbacked') and adata.isbacked:
                obs_values = adata.obs.iloc[cell_indices_array][groupby]
            else:
                obs_values = adata.obs[cell_indices][groupby]
        else:
            obs_values = adata.obs[groupby]
    
    # Apply transformations
    if transformation:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method=transformation)
    
    # Create cached data structure
    cached_data = {
        'gene_df': gene_df,
        'obs_values': obs_values,
        'valid_genes': valid_genes,
        'extraction_time': time.time() - start_time
    }
    
    # Cache with size limit (keep last 50 datasets)
    if len(_violin_data_cache) > 50:
        # Remove oldest entries
        oldest_key = next(iter(_violin_data_cache))
        del _violin_data_cache[oldest_key]
    
    _violin_data_cache[cache_key] = cached_data
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
    
    # Set up colors
    unique_labels = sorted(adata_obs[groupby].unique()) if adata_obs is not None else sorted(adata.obs[groupby].unique())
    if groupby_label_color_map is None:
        groupby_label_color_map = {label: default_color[i % len(default_color)] for i, label in enumerate(unique_labels)}

    # Create subplots
    num_genes = len(valid_genes)
    fig_height = max(400, num_genes * 140)
    # Keep a consistent ~28px gap between rows so adjacent y-axis tick labels (the
    # 0.0 of one row and the top tick of the next) don't overlap, for any number
    # of genes. vertical_spacing is a fraction of height, capped to Plotly's max.
    target_gap_px = 28
    vertical_spacing = min(target_gap_px / fig_height, 0.9 / max(num_genes - 1, 1))
    fig = make_subplots(rows=num_genes, cols=1, shared_xaxes=True, vertical_spacing=vertical_spacing)

    # Reserve a fixed left label column sized to the longest gene name so every
    # subplot's y-axis starts at the same x and the (horizontal) labels never
    # overlap. tick_pad leaves room for the y-axis tick numbers between the label
    # and the axis; the label is right-aligned just to the left of them.
    tick_pad = 40
    max_label_len = max((len(str(g)) for g in valid_genes), default=1)
    left_margin = int(min(360, tick_pad + max_label_len * 7 + 15))

    # Pre-calculate boolean masks for each label to avoid repeated groupby operations
    # Use numpy array for faster boolean indexing
    obs_values_array = obs_values.to_numpy() if hasattr(obs_values, 'to_numpy') else np.array(obs_values)
    label_masks = {}
    for label in labels:
        mask = (obs_values_array == label)
        if np.any(mask):
            label_masks[label] = mask

    # Fixed seed so each violin's subsample (and thus its shape) is stable across
    # re-renders.
    rng = np.random.default_rng(0)

    # Create violin plots
    for i, gene in enumerate(valid_genes):
        # Access gene values directly as numpy array
        gene_values = gene_df[gene].values
        expr_max = gene_values.max()

        for label in labels:
            if label in label_masks:
                mask = label_masks[label]
                group_expr = gene_values[mask]

                # Subsample large groups: a KDE over this many points matches the
                # full distribution but keeps payload/compute bounded. scalemode is
                # 'width', so the reduced count doesn't change the violin's width.
                if group_expr.shape[0] > MAX_VIOLIN_POINTS:
                    group_expr = rng.choice(group_expr, MAX_VIOLIN_POINTS, replace=False)

                # Calculate variance to determine bandwidth (seaborn-like behavior)
                variance = np.var(group_expr)
                if variance < 1e-10:  # Near-zero variance
                    spanmode = 'soft'
                else:
                    # Normal bandwidth for data with variance
                    spanmode = 'hard'
                
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
                        spanmode=spanmode,
                        fillcolor=_to_rgba(groupby_label_color_map[label], VIOLIN_FILL_ALPHA),
                        line_color='DarkSlateGrey',
                        hoveron='violins',
                        hoverinfo='y',
                        jitter=0.05,
                        scalemode='width',
                        scalegroup=label,
                        marker=dict(size=2),
                    ),
                    row=i + 1, col=1
                )

        expr_min = gene_values.min()
        if expr_max == expr_min:
            pad = 1e-3 if expr_max == 0 else abs(expr_max) * 0.05
            y_range = [expr_min - pad, expr_max + pad]
        else:
            y_range = [expr_min, expr_max]
        fig.update_yaxes(
            range=y_range,
            tickformat=".1f",
            row=i + 1, col=1
        )
        # Horizontal gene label in the reserved left column, right-aligned just to
        # the left of the y-axis tick numbers, vertically centered on the subplot.
        yaxis_domain_ref = "y domain" if i == 0 else f"y{i + 1} domain"
        fig.add_annotation(
            text=str(gene),
            xref="paper", x=0, xanchor="right", xshift=-tick_pad,
            yref=yaxis_domain_ref, y=0.5, yanchor="middle",
            showarrow=False,
            font=dict(size=12),
        )
    # Layout
    fig.update_layout(
    plot_bgcolor='white', paper_bgcolor='white', font=dict(size=10),
    showlegend=False,
    height=fig_height,
    margin=dict(l=left_margin, r=20, t=20, b=40)
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    
    return fig
