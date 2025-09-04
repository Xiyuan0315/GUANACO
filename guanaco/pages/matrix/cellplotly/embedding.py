import plotly.graph_objs as go
import pandas as pd
import numpy as np
from .gene_extraction_utils import extract_gene_expression, apply_transformation
import plotly.express as px

def plot_continuous_embedding(
    adata, embedding_key, color, x_axis=None, y_axis=None,
    transformation=None, order=None, color_map='Viridis',
    marker_size=5, opacity=1, annotation=None, axis_show=True
):
    """
    Plot a continuous feature (e.g., gene expression) on a 2D embedding.
    """
    embedding_prefixes = {
        'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
        'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
    }
    embedding_prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    embedding_data = adata.obsm[embedding_key]

    # Set column names for the embedding
    num_dimensions = embedding_data.shape[1]
    embedding_columns = [f'{embedding_prefix}{i + 1}' for i in range(num_dimensions)]
    embedding_df = pd.DataFrame(embedding_data, columns=embedding_columns)

    # Default x and y axis
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])

    # Extract gene expression using optimized method
    gene_expression = extract_gene_expression(adata, color)
    
    # Apply transformation if specified
    if transformation:
        gene_expression = apply_transformation(gene_expression, transformation, copy=False)

    embedding_df[color] = gene_expression

    # Add original cell indices BEFORE sorting/shuffling
    embedding_df['_original_idx'] = np.arange(len(embedding_df))
    
    # Sort order
    if order == 'max':
        embedding_df_sorted = embedding_df.sort_values(by=color)
    elif order == 'min':
        embedding_df_sorted = embedding_df.sort_values(by=color, ascending=False)
    elif order == 'random':
        embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
    else:
        embedding_df_sorted = embedding_df

    # Add annotation if provided
    if annotation:
        embedding_df_sorted[annotation] = adata.obs[annotation].values

    # Use the original indices for selection tracking (not the sorted indices)
    embedding_df_sorted['_cell_idx'] = embedding_df_sorted['_original_idx']

    # Create scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=embedding_df_sorted[x_axis],
        y=embedding_df_sorted[y_axis],
        mode='markers',
        marker=dict(
            color=embedding_df_sorted[color],
            colorscale=color_map,
            cmin=embedding_df_sorted[color].min(),
            cmax=embedding_df_sorted[color].max(),
            size=marker_size,
            opacity=opacity,
            colorbar=dict(
                title=transformation if transformation else color,
                len = 0.8
            )
        ),
        customdata=embedding_df_sorted['_cell_idx'],  # Add cell indices for selection
        hoverinfo='skip',  # Disable hover
        selectedpoints=None,  # Enable selection
        selected=dict(marker=dict(opacity=1)),  # Keep selected points fully visible
        unselected=dict(marker=dict(opacity=0.2))  # Dim unselected points
    ))

    fig.update_layout(
        autosize=True,
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(
            text=f'<b>{color}</b>',
            x=0.5,
            y=0.98,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title=x_axis,
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False,
            zeroline=False
        ),
        margin=dict(t=40, r=20, l=50, b=50)
    )

    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )

    return fig


def plot_categorical_embedding_with_fixed_colors(
    adata, adata_full, gene, embedding_key, color,
    x_axis=None, y_axis=None,
    color_map=None, marker_size=5, opacity=1,
    legend_show='on legend', axis_show=True
):

    # Get all unique labels from the full dataset to ensure consistent colors
    all_unique_labels = sorted(adata_full.obs[color].unique())
    
    # Create color mapping for ALL categories (not just filtered ones)
    color_map = color_map or px.colors.qualitative.Plotly
    label_to_color_dict = {
        label: color_map[i % len(color_map)]
        for i, label in enumerate(all_unique_labels)
    }
    
    embedding_prefixes = {
        "X_umap": "UMAP", "X_pca": "PCA", "X_tsne": "t-SNE",
        "X_diffmap": "DiffMap", "X_phate": "PHATE", "X_draw_graph_fa": "FA"
    }
    on_data = legend_show == 'on data'

    # Prepare embedding coordinates
    embedding_data = adata.obsm[embedding_key]
    prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    dims = [f"{prefix}{i+1}" for i in range(embedding_data.shape[1])]
    x_axis = x_axis or dims[0]
    y_axis = y_axis or (dims[1] if len(dims) > 1 else dims[0])

    # Prepare DataFrame
    df = pd.DataFrame(embedding_data, columns=dims)
    df[color] = adata.obs[color].values
    
    # Only extract gene expression if gene is provided
    if gene is not None and gene in adata.var_names:
        from guanaco.pages.matrix.cellplotly.gene_extraction_utils import extract_gene_expression
        df[gene] = extract_gene_expression(adata, gene)

    # Get unique labels in the filtered data
    unique_labels_filtered = sorted(df[color].unique())

    fig = go.Figure()
    
    # First, add a grey background trace for all cells
    fig.add_trace(go.Scattergl(
        x=df[x_axis],
        y=df[y_axis],
        mode='markers',
        marker=dict(
            size=marker_size,
            color='lightgrey',
            opacity=opacity * 0.3,
        ),
        name='Background',
        hoverinfo='skip',
        showlegend=False,
        visible=True
    ))

    # Add one trace per category (only for categories present in filtered data)
    for label in unique_labels_filtered:
        mask = df[color] == label
        fig.add_trace(go.Scattergl(
            x=df.loc[mask, x_axis],
            y=df.loc[mask, y_axis],
            mode='markers',
            marker=dict(
                size=marker_size,
                color=label_to_color_dict[label],  # Use color from full dataset mapping
                opacity=opacity,
            ),
            name=str(label),
            customdata=df.loc[mask, color] if gene is None else np.stack([df.loc[mask, color], df.loc[mask, gene]], axis=-1),
            hoverinfo='skip',  # Disable hover
            showlegend=not on_data,
            legendgroup=str(label),
        ))

    # Add labels at cluster medians if requested
    if on_data:
        for label in unique_labels_filtered:
            mask = df[color] == label
            median_x = df.loc[mask, x_axis].median()
            median_y = df.loc[mask, y_axis].median()
            fig.add_annotation(
                x=median_x, y=median_y,
                text=f"<b>{label}</b>",
                showarrow=False,
                font=dict(size=12, color='black'),
                xanchor='center', yanchor='middle',
                opacity=0.9,
            )

    # Layout settings
    fig.update_layout(
        autosize=True,
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(text=f"<b>{color}</b>", x=0.5, y=0.98, xanchor='center', yanchor='top'),
        xaxis=dict(
            title=x_axis,
            showgrid=False, zeroline=False,
            tickfont=dict(color="rgba(0,0,0,0)" if not axis_show else "black")
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False, zeroline=False,
            tickfont=dict(color="rgba(0,0,0,0)" if not axis_show else "black")
        ),
        legend=dict(
            orientation='v',
            itemsizing='constant',
            x=1.02, y=0.5,
            bgcolor='rgba(0,0,0,0)',
            itemclick='toggle',
            itemdoubleclick='toggleothers',
            font=dict(size=12)
        ) if not on_data else None,
        margin=dict(t=40, r=100, l=50, b=50)
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

    return fig


def plot_coexpression_embedding(
    adata, embedding_key, gene1, gene2,
    x_axis=None, y_axis=None,
    threshold1=0.5, threshold2=0.5,
    transformation=None,
    color_map=None,
    marker_size=5, opacity=1,
    legend_show='right', axis_show=True
):
    """
    Plot co-expression of two genes on a 2D embedding.
    Cells are categorized into 4 groups based on expression thresholds.
    """
    embedding_prefixes = {
        'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
        'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
    }
    embedding_prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    embedding_data = adata.obsm[embedding_key]

    # Set column names for the embedding
    num_dimensions = embedding_data.shape[1]
    embedding_columns = [f'{embedding_prefix}{i + 1}' for i in range(num_dimensions)]
    embedding_df = pd.DataFrame(embedding_data, columns=embedding_columns)

    # Default x and y axis
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])

    # Extract gene expression for both genes
    gene1_expr = extract_gene_expression(adata, gene1)
    gene2_expr = extract_gene_expression(adata, gene2)
    
    # Apply transformation if specified
    if transformation:
        gene1_expr = apply_transformation(gene1_expr, transformation, copy=True)
        gene2_expr = apply_transformation(gene2_expr, transformation, copy=True)

    # Create categories based on thresholds
    # Use threshold values directly (they are already actual expression values from the callback)
    gene1_threshold = threshold1
    gene2_threshold = threshold2
    
    # Categorize cells
    gene1_expressed = gene1_expr > gene1_threshold
    gene2_expressed = gene2_expr > gene2_threshold
    
    # Create category labels
    categories = np.zeros(len(adata), dtype='object')
    categories[:] = 'Neither'
    categories[gene1_expressed & ~gene2_expressed] = f'{gene1} only'
    categories[~gene1_expressed & gene2_expressed] = f'{gene2} only'
    categories[gene1_expressed & gene2_expressed] = 'Co-expressed'
    
    embedding_df['category'] = categories
    embedding_df['gene1_expr'] = gene1_expr
    embedding_df['gene2_expr'] = gene2_expr

    # Define color mapping for 4 categories
    if color_map is None:
        color_map = {
            'Neither': '#E8E8E8',  # Light gray
            f'{gene1} only': '#648fff',  # Blue (color-blind friendly)
            f'{gene2} only': '#ffb000',  # Orange (color-blind friendly)
            'Co-expressed': '#dc267f'  # Magenta (color-blind friendly)
        }
    
    # Ensure all 4 categories appear in order
    category_order = ['Neither', f'{gene1} only', f'{gene2} only', 'Co-expressed']
    
    fig = go.Figure()

    # Add traces for each category
    for category in category_order:
        mask = embedding_df['category'] == category
        if mask.any():
            fig.add_trace(go.Scattergl(
                x=embedding_df.loc[mask, x_axis],
                y=embedding_df.loc[mask, y_axis],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=color_map.get(category, '#808080'),
                    opacity=opacity,
                ),
                name=category,
                customdata=np.arange(len(adata))[mask],  # Add cell indices for selection
                showlegend=(legend_show == 'right'),
                hoverinfo='skip',
                selectedpoints=None,  # Enable selection
                selected=dict(marker=dict(opacity=1)),  # Keep selected points fully visible
                unselected=dict(marker=dict(opacity=0.2))  # Dim unselected points
            ))

    # Add category labels on plot if requested
    if legend_show == 'on data':
        for category in category_order:
            mask = embedding_df['category'] == category
            if mask.any():
                median_x = embedding_df.loc[mask, x_axis].median()
                median_y = embedding_df.loc[mask, y_axis].median()
                fig.add_annotation(
                    x=median_x, y=median_y,
                    text=f'<b>{category}</b>',
                    showarrow=False,
                    font=dict(size=10, color='black'),
                    xanchor='center', yanchor='middle',
                    opacity=0.9,
                )

    # Update layout
    fig.update_layout(
        autosize=True,
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(
            text=f'<b>Co-expression: {gene1} & {gene2}</b>',
            x=0.5,
            y=0.98,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title=x_axis,
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False,
            zeroline=False
        ),
        legend=dict(
            orientation='v',
            itemsizing='constant',
            x=1.02, y=0.5,
            bgcolor='rgba(0,0,0,0)',
            itemclick='toggle',
            itemdoubleclick='toggleothers',
            font=dict(size=10)
        ) if legend_show == 'right' else None,
        margin=dict(t=40, r=100, l=50, b=50)
    )

    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )

    return fig


def plot_continuous_annotation(
    adata, embedding_key, annotation, x_axis=None, y_axis=None,
    transformation=None, order=None, color_map='Viridis',
    marker_size=5, opacity=1, axis_show=True
):
    """
    Plot a continuous annotation (from obs) on a 2D embedding.
    Modified to work with obs columns instead of gene expression.
    """

    
    embedding_prefixes = {
        'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
        'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
    }
    embedding_prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    embedding_data = adata.obsm[embedding_key]

    # Set column names for the embedding
    num_dimensions = embedding_data.shape[1]
    embedding_columns = [f'{embedding_prefix}{i + 1}' for i in range(num_dimensions)]
    embedding_df = pd.DataFrame(embedding_data, columns=embedding_columns)

    # Default x and y axis
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])

    # Extract annotation values (from obs instead of expression)
    annotation_values = adata.obs[annotation].values
    
    # Only apply transformations if explicitly requested and data is numeric
    # For annotation data, we usually want to see the raw values
    if transformation and annotation_values.dtype in ['float32', 'float64', 'int32', 'int64']:
        if transformation == 'log':
            # Handle negative values for log transformation
            min_val = annotation_values.min()
            if min_val <= 0:
                annotation_values = annotation_values - min_val + 1
            annotation_values = np.log1p(annotation_values)
        elif transformation == 'z_score':
            annotation_values = (annotation_values - np.mean(annotation_values)) / np.std(annotation_values)

    embedding_df[annotation] = annotation_values
    # Add cell indices for selection tracking
    embedding_df['_cell_idx'] = np.arange(len(embedding_df))

    # Sort order
    if order == 'max':
        embedding_df_sorted = embedding_df.sort_values(by=annotation)
    elif order == 'min':
        embedding_df_sorted = embedding_df.sort_values(by=annotation, ascending=False)
    elif order == 'random':
        embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
    else:
        embedding_df_sorted = embedding_df

    # Create scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=embedding_df_sorted[x_axis],
        y=embedding_df_sorted[y_axis],
        mode='markers',
        marker=dict(
            color=embedding_df_sorted[annotation],
            colorscale=color_map,
            cmin=embedding_df_sorted[annotation].min(),
            cmax=embedding_df_sorted[annotation].max(),
            size=marker_size,
            opacity=opacity,
            colorbar=dict(
                title=f"{annotation}<br>{transformation if transformation else ''}",
                len=0.8
            )
        ),
        customdata=np.stack([embedding_df_sorted[annotation], embedding_df_sorted['_cell_idx']], axis=-1),  # Add customdata with cell index
        hoverinfo='skip',  # Disable hover
        selectedpoints=None,  # Enable selection
        selected=dict(marker=dict(opacity=1)),  # Keep selected points fully visible
        unselected=dict(marker=dict(opacity=0.2))  # Dim unselected points
    ))

    fig.update_layout(
        autosize=True,
        plot_bgcolor='white',
        paper_bgcolor='white',
        title=dict(
            text=f'<b>{annotation}</b>',
            x=0.5,
            y=0.98,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title=x_axis,
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title=y_axis,
            showgrid=False,
            zeroline=False
        ),
        margin=dict(t=40, r=20, l=50, b=50)
    )

    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
    )

    return fig
