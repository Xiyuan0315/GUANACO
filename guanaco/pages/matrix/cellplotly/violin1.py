from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from guanaco.pages.matrix.cellplotly.gene_extraction_utils import (
    extract_gene_expression, extract_multiple_genes, apply_transformation
)
import hashlib
import time

# Global cache for violin plot data
_violin_data_cache = {}
_violin_plot_cache = {}

default_color = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"]

def _create_cache_key(genes, labels, groupby, transformation, adata_id, adata_obs_id=None):
    """Create a unique cache key for violin data."""
    key_data = {
        'genes': tuple(sorted(genes)),
        'labels': tuple(sorted(labels)) if labels else None,
        'groupby': groupby,
        'transformation': transformation,
        'adata_id': adata_id,
        'adata_obs_id': adata_obs_id
    }
    return hashlib.md5(str(key_data).encode()).hexdigest()

def _get_adata_id(adata):
    """Get a unique identifier for the adata object."""
    if hasattr(adata, 'isbacked') and adata.isbacked and hasattr(adata, 'filename'):
        return adata.filename
    else:
        return f"{id(adata)}_{adata.shape}"

def _extract_and_cache_violin_data(adata, genes, labels, groupby, transformation):
    """Extract and cache violin plot data for reuse."""
    adata_id = _get_adata_id(adata)
    cache_key = _create_cache_key(genes, labels, groupby, transformation, adata_id)
    
    # Check if data is already cached
    if cache_key in _violin_data_cache:
        return _violin_data_cache[cache_key]
    
    # Extract data
    start_time = time.time()
    
    # Filter cells based on selected labels
    if labels:
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
            if labels and not (hasattr(adata, 'isbacked') and adata.isbacked):
                filtered_adata = adata[cell_indices]
                gene_df = extract_multiple_genes(filtered_adata, valid_genes)
                obs_values = filtered_adata.obs[groupby]
            else:
                gene_df = extract_multiple_genes(adata, valid_genes)
                if cell_indices_array is not None:
                    gene_df = gene_df.iloc[cell_indices_array]
                    obs_values = adata.obs.iloc[cell_indices_array][groupby]
                else:
                    obs_values = adata.obs[groupby]
        except Exception:
            # Fallback to individual extraction
            gene_data = {}
            for gene in valid_genes:
                gene_expr = extract_gene_expression(adata, gene)
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
        gene_expr = extract_gene_expression(adata, gene)
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

def plot_violin1(adata, genes, groupby, labels=None, transformation=None, show_box=False, show_points=False, groupby_label_color_map=None, adata_obs=None):
   
    # if no label or no genes, stop updating.
    if len(genes) == 0 or labels is None or len(labels) == 0:
        raise PreventUpdate
    # Check plot cache
    adata_id = _get_adata_id(adata)
    adata_obs_id = id(adata_obs) if adata_obs is not None else None
    plot_cache_key = _create_cache_key(genes, labels, groupby, transformation, adata_id, adata_obs_id) + f"_{show_box}_{show_points}"
    
    if plot_cache_key in _violin_plot_cache:
        cached_plot = _violin_plot_cache[plot_cache_key]
        if groupby_label_color_map:
            cached_plot = _update_plot_colors(cached_plot, groupby_label_color_map)
        return cached_plot
    
    # Extract data
    cached_data = _extract_and_cache_violin_data(adata, genes, labels, groupby, transformation)
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
    fig = make_subplots(rows=num_genes, cols=1, shared_xaxes=True, vertical_spacing=0.02)

    # Pre-calculate boolean masks for each label to avoid repeated groupby operations
    # Use numpy array for faster boolean indexing
    obs_values_array = obs_values.to_numpy() if hasattr(obs_values, 'to_numpy') else np.array(obs_values)
    label_masks = {}
    for label in labels:
        mask = (obs_values_array == label)
        if np.any(mask):
            label_masks[label] = mask

    # Create violin plots
    for i, gene in enumerate(valid_genes):
        # Access gene values directly as numpy array
        gene_values = gene_df[gene].values
        expr_max = gene_values.max()

        points_mode = 'all' if show_points else False
        
        for label in labels:
            if label in label_masks:
                mask = label_masks[label]
                group_expr = gene_values[mask]
                
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
                        box_visible=show_box,
                        points=points_mode,
                        meanline_visible=True,
                        showlegend=(i == 0),
                        name=label,
                        spanmode=spanmode,
                        fillcolor=groupby_label_color_map[label],
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
        
        fig.update_yaxes(
            range=[0, expr_max],
            title_text=gene,
            title_standoff=10,
            tickformat=".1f",
            row=i + 1, col=1
        )
    # Layout
    fig_height = max(400, num_genes * 120)
    fig.update_layout(
    plot_bgcolor='white', paper_bgcolor='white', font=dict(size=10),
    showlegend=False,
    height=fig_height,
    margin=dict(l=80, r=20, t=20, b=40)
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    
    # Cache plot
    if len(_violin_plot_cache) > 20:
        oldest_key = next(iter(_violin_plot_cache))
        del _violin_plot_cache[oldest_key]
    
    _violin_plot_cache[plot_cache_key] = fig
    return fig

def _update_plot_colors(fig, new_color_map):
    # Create a copy to avoid modifying cached plot
    new_fig = go.Figure(fig)
    
    for trace in new_fig.data:
        if hasattr(trace, 'name') and trace.name in new_color_map:
            trace.fillcolor = new_color_map[trace.name]
    
    return new_fig

def clear_violin_cache():
    global _violin_data_cache, _violin_plot_cache
    _violin_data_cache.clear()
    _violin_plot_cache.clear()

def get_violin_cache_info():
    return {
        'data_cache_size': len(_violin_data_cache),
        'plot_cache_size': len(_violin_plot_cache),
        'data_cache_keys': list(_violin_data_cache.keys())[:5],  # Show first 5
        'plot_cache_keys': list(_violin_plot_cache.keys())[:5]   # Show first 5
    }
