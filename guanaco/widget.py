"""
Convenience wrappers that return Plotly FigureWidget objects for GUANACO plots.

Usage (after installing the package):

    from guanaco.widget import embedding_continuous_widget
    fw = embedding_continuous_widget(adata, embedding_key='X_umap', gene='TNF')
    fw  # in Jupyter, this renders an interactive widget

The wrappers lazily import the underlying plotting functions from
guanaco.pages.matrix.cellplotly.* and convert the returned go.Figure to
go.FigureWidget for use in notebooks.
"""

from __future__ import annotations

import plotly.graph_objects as go


def _to_widget(fig: "go.Figure") -> "go.FigureWidget":
    """Convert a Plotly Figure to FigureWidget."""
    return go.FigureWidget(fig)


# -------- Embedding plots --------
def embedding_continuous_widget(
    adata,
    embedding_key: str,
    gene: str,
    x_axis: str | None = None,
    y_axis: str | None = None,
    transformation: str | None = None,
    order: str | None = None,
    color_map: str = 'Viridis',
    marker_size: int = 5,
    opacity: float = 1.0,
    annotation: str | None = None,
    axis_show: bool = True,
):
    """Return a FigureWidget for a continuous embedding (e.g., gene expression)."""
    from guanaco.pages.matrix.cellplotly.embedding import plot_continuous_embedding

    fig = plot_continuous_embedding(
        adata=adata,
        embedding_key=embedding_key,
        color=gene,
        x_axis=x_axis,
        y_axis=y_axis,
        transformation=transformation,
        order=order,
        color_map=color_map,
        marker_size=marker_size,
        opacity=opacity,
        annotation=annotation,
        axis_show=axis_show,
    )
    return _to_widget(fig)


def embedding_categorical_widget(
    adata,
    embedding_key: str,
    color: str,
    x_axis: str | None = None,
    y_axis: str | None = None,
    color_map=None,
    marker_size: int = 5,
    opacity: float = 1.0,
    legend_show: str = 'on legend',
    axis_show: bool = True,
):
    """Return a FigureWidget for a categorical embedding with fixed colors."""
    from guanaco.pages.matrix.cellplotly.embedding import plot_categorical_embedding_with_fixed_colors

    fig = plot_categorical_embedding_with_fixed_colors(
        adata=adata,
        adata_full=adata,
        embedding_key=embedding_key,
        color=color,
        x_axis=x_axis,
        y_axis=y_axis,
        color_map=color_map,
        marker_size=marker_size,
        opacity=opacity,
        legend_show=legend_show,
        axis_show=axis_show,
    )
    return _to_widget(fig)


def embedding_coexpression_widget(
    adata,
    embedding_key: str,
    gene1: str,
    gene2: str,
    threshold1 = 0.5,
    threshold2 = 0.5,
    x_axis: str | None = None,
    y_axis: str | None = None,
    transformation: str | None = None,
    marker_size: int = 5,
    opacity: float = 1.0,
    legend_show: str = 'on legend',
):
    """Return a FigureWidget for two-gene coexpression categories on embeddings."""
    from guanaco.pages.matrix.cellplotly.embedding import plot_coexpression_embedding

    fig = plot_coexpression_embedding(
        adata=adata,
        embedding_key=embedding_key,
        gene1=gene1,
        gene2=gene2,
        threshold1=threshold1,
        threshold2=threshold2,
        x_axis=x_axis,
        y_axis=y_axis,
        transformation=transformation,
        marker_size=marker_size,
        opacity=opacity,
        legend_show=legend_show,
    )
    return _to_widget(fig)


def embedding_continuous_annotation_widget(
    adata,
    embedding_key: str,
    annotation: str,
    x_axis: str | None = None,
    y_axis: str | None = None,
    color_map: str = 'Viridis',
    marker_size: int = 5,
    opacity: float = 1.0,
    axis_show: bool = True,
):
    """Return a FigureWidget for continuous annotations on embeddings."""
    from guanaco.pages.matrix.cellplotly.embedding import plot_continuous_annotation

    fig = plot_continuous_annotation(
        adata=adata,
        embedding_key=embedding_key,
        annotation=annotation,
        x_axis=x_axis,
        y_axis=y_axis,
        color_map=color_map,
        marker_size=marker_size,
        opacity=opacity,
        axis_show=axis_show,
    )
    return _to_widget(fig)


# -------- Heatmap --------
def heatmap_widget(
    adata,
    genes,
    groupby1: str,
    groupby2: str | None = None,
    labels=None,
    log: bool = False,
    z_score: bool = False,
    boundary: bool = False,
    color_map='Viridis',
    groupby1_label_color_map=None,
    groupby2_label_color_map=None,
    max_cells: int = 50000,
    n_bins: int = 10000,
    transformation=None,
    adata_obs=None,
):
    """Return a FigureWidget for the unified heatmap plot."""
    from guanaco.pages.matrix.cellplotly.heatmap import plot_unified_heatmap

    fig = plot_unified_heatmap(
        adata=adata,
        genes=genes,
        groupby1=groupby1,
        groupby2=groupby2,
        labels=labels,
        log=log,
        z_score=z_score,
        boundary=boundary,
        color_map=color_map,
        groupby1_label_color_map=groupby1_label_color_map,
        groupby2_label_color_map=groupby2_label_color_map,
        max_cells=max_cells,
        n_bins=n_bins,
        transformation=transformation,
        adata_obs=adata_obs,
    )
    return _to_widget(fig)


# -------- Violin plots --------
def violin1_widget(
    adata,
    genes,
    groupby: str,
    labels = None,
    transformation=None,
    show_box: bool = False,
    show_points: bool = False,
    groupby_label_color_map=None,
    adata_obs=None,
):
    """Return a FigureWidget for the first violin plot variant."""
    from guanaco.pages.matrix.cellplotly.violin1 import plot_violin1

    fig = plot_violin1(
        adata=adata,
        genes=genes,
        labels=labels,
        groupby=groupby,
        transformation=transformation,
        show_box=show_box,
        show_points=show_points,
        groupby_label_color_map=groupby_label_color_map,
        adata_obs=adata_obs,
    )
    return _to_widget(fig)


def violin2_widget(
    adata,
    key: str,
    meta1: str,
    meta2: str,
    mode: str,
    transformation = None,
    labels = None,
    color_map = None
):
    """Return a FigureWidget for the second violin plot variant."""
    from guanaco.pages.matrix.cellplotly.violin2 import plot_violin2_new

    fig = plot_violin2_new(
        adata, key, meta1, meta2, mode,
        transformation= transformation,
        show_box=True,
        show_points=True,
        test_method='auto',
        labels=labels,
        color_map=color_map)
    
    return _to_widget(fig)


# -------- Dot matrix --------
def dotmatrix_widget(
    adata,
    genes,
    groupby: str,
    labels=None,
    color_map='Viridis',
    transformation=None,
    plot_type='dotplot',
    cluster='none', 
    method='average', 
    metric='correlation'
):
    """Return a FigureWidget for the dot matrix plot."""
    from guanaco.pages.matrix.cellplotly.dotmatrix import plot_dot_matrix
    fig = plot_dot_matrix(
        adata=adata,
        genes=genes,
        groupby=groupby,
        selected_labels=labels,
        color_map=color_map,
        transformation=transformation,
        plot_type=plot_type,
        cluster=cluster,
        method=method,
        metric=metric,
    )
    return _to_widget(fig)


# -------- Stacked bar --------
def stacked_bar_widget(
    adata,
    x_meta: str,
    y_meta: str,
    norm = 'prop',
    color_map=None,
    y_order=None,
    x_order=None,
):
    """Return a FigureWidget for the stacked bar plot (or histogram when x_meta == y_meta)."""
    from guanaco.pages.matrix.cellplotly.stacked_bar import plot_stacked_bar

    fig = plot_stacked_bar(
        x_meta=x_meta,
        y_meta=y_meta,
        norm=norm,
        adata=adata,
        color_map=color_map,
        y_order=y_order,
        x_order=x_order,
    )
    return _to_widget(fig)


# -------- Pseudotime --------
def pseudotime_widget(
    adata,
    genes,
    pseudotime_key: str = 'pseudotime',
    groupby: str | None = None,
    min_expr: float = 0.5,
    transformation: str = 'none',
    color_map=None,
):
    """Return a FigureWidget for gene expression along pseudotime."""
    from guanaco.pages.matrix.cellplotly.pseudotime import plot_genes_in_pseudotime

    fig = plot_genes_in_pseudotime(
        adata=adata,
        genes=genes,
        pseudotime_key=pseudotime_key,
        groupby=groupby,
        min_expr=min_expr,
        transformation=transformation,
        color_map=color_map,
    )
    return _to_widget(fig)


__all__ = [
    # Embedding
    'embedding_continuous_widget',
    'embedding_categorical_widget',
    'embedding_coexpression_widget',
    'embedding_continuous_annotation_widget',
    # Heatmap
    'heatmap_widget',
    # Violins
    'violin1_widget',
    'violin2_widget',
    # Dot matrix
    'dotmatrix_widget',
    # Stacked bar
    'stacked_bar_widget',
    # Pseudotime
    'pseudotime_widget',
]

