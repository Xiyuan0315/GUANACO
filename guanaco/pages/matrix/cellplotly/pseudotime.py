import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import plotly.express as px
from scipy.interpolate import UnivariateSpline

try:
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.pipeline import make_pipeline
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


def plot_genes_in_pseudotime(adata, genes, pseudotime_key='pseudotime', groupby=None, 
                           min_expr=0.5, transformation='none', color_map=None):
    """Plot gene expression along pseudotime with smoothed curves."""
    
    # Filter out genes that don't exist in the dataset
    valid_genes = [gene for gene in genes if gene in adata.var_names]
    if not valid_genes:
        # Return empty figure if no valid genes
        fig = go.Figure()
        fig.add_annotation(
            text="No valid genes found in the dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=400
        )
        return fig
    
    # Check if pseudotime exists
    if pseudotime_key not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text=f"'{pseudotime_key}' not found in adata.obs",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=400
        )
        return fig
    
    # Extract data
    if hasattr(adata, 'isbacked') and adata.isbacked:
        # Extract gene indices for valid genes
        gene_indices = [adata.var_names.get_loc(gene) for gene in valid_genes]
        # Extract expression data directly from the backed file
        if hasattr(adata.X, 'toarray'):
            gene_expression_matrix = adata.X[:, gene_indices].toarray()
        else:
            gene_expression_matrix = adata.X[:, gene_indices]
    else:
        # Non-backed AnnData
        adata_selected = adata[:, valid_genes]
        gene_expression_matrix = adata_selected.X.toarray()
    
    # Apply transformation
    if transformation == 'log':
        gene_expression_matrix = np.log1p(gene_expression_matrix)
    elif transformation == 'z_score':
        gene_expression_matrix = (gene_expression_matrix - gene_expression_matrix.mean(axis=0)) / gene_expression_matrix.std(axis=0)
    
    # Create DataFrame
    expr_df = pd.DataFrame(gene_expression_matrix, columns=valid_genes)
    expr_df['pseudotime'] = adata.obs[pseudotime_key].values
    expr_df['cell_id'] = adata.obs_names
    
    if groupby and groupby in adata.obs.columns:
        expr_df[groupby] = adata.obs[groupby].values
    
    # Filter cells based on minimum expression
    if min_expr > 0:
        # Keep cells where at least one gene passes the threshold
        mask = (gene_expression_matrix >= min_expr).any(axis=1)
        expr_df = expr_df[mask]
    
    # Remove cells with NaN pseudotime
    expr_df = expr_df.dropna(subset=['pseudotime'])
    
    if expr_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No cells pass the filtering criteria",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=400
        )
        return fig
    
    # Sort by pseudotime
    expr_df = expr_df.sort_values('pseudotime')
    
    # Setup colors
    if groupby and groupby in expr_df.columns:
        unique_groups = sorted(expr_df[groupby].unique())
        if color_map is None:
            colors = px.colors.qualitative.Plotly
            color_map = dict(zip(unique_groups, colors[:len(unique_groups)]))
    else:
        unique_groups = ['All cells']
        expr_df['_group'] = 'All cells'
        groupby = '_group'
        color_map = {'All cells': '#1f77b4'}
    
    # Create subplots
    n_genes = len(valid_genes)
    row_height = 300  # Height per gene plot
    
    fig = make_subplots(
        rows=n_genes, cols=1,
        subplot_titles=valid_genes,
        shared_xaxes=True,
        vertical_spacing=0.05,
        row_heights=[row_height] * n_genes
    )
    
    # Plot each gene
    for gene_idx, gene in enumerate(valid_genes):
        row = gene_idx + 1
        
        # Add scatter points for each group
        for group in unique_groups:
            group_data = expr_df[expr_df[groupby] == group]
            
            if not group_data.empty:
                fig.add_trace(
                    go.Scatter(
                        x=group_data['pseudotime'],
                        y=group_data[gene],
                        mode='markers',
                        marker=dict(
                            color=color_map[group],
                            size=3,
                            opacity=0.6
                        ),
                        name=str(group),
                        legendgroup=str(group),
                        showlegend=(gene_idx == 0),  # Only show legend for first gene
                        hoverinfo='skip'
                    ),
                    row=row, col=1
                )
        
        # Add smoothed curve for all cells combined
        if len(expr_df) > 10:  # Need enough points for smoothing
            # Prepare data for smoothing
            x_smooth = expr_df['pseudotime'].values
            y_smooth = expr_df[gene].values
            
            # Create smooth curve using polynomial regression
            try:
                if SKLEARN_AVAILABLE:
                    # Use polynomial regression for smooth curve
                    poly_degree = min(5, len(expr_df) // 10)  # Adaptive degree
                    model = make_pipeline(
                        PolynomialFeatures(poly_degree),
                        Ridge(alpha=1.0)
                    )
                    model.fit(x_smooth.reshape(-1, 1), y_smooth)
                    
                    # Generate smooth curve
                    x_range = np.linspace(x_smooth.min(), x_smooth.max(), 200)
                    y_pred = model.predict(x_range.reshape(-1, 1))
                else:
                    # Fallback to simple spline if sklearn not available
                    # Remove duplicates and sort
                    df_smooth = pd.DataFrame({'x': x_smooth, 'y': y_smooth})
                    df_smooth = df_smooth.groupby('x').mean().reset_index()
                    
                    if len(df_smooth) > 3:
                        spline = UnivariateSpline(df_smooth['x'], df_smooth['y'], s=len(df_smooth))
                        x_range = np.linspace(x_smooth.min(), x_smooth.max(), 200)
                        y_pred = spline(x_range)
                    else:
                        x_range = df_smooth['x'].values
                        y_pred = df_smooth['y'].values
                
                # Add smoothed curve
                fig.add_trace(
                    go.Scatter(
                        x=x_range,
                        y=y_pred,
                        mode='lines',
                        line=dict(
                            color='black',
                            width=3
                        ),
                        name='Smoothed',
                        legendgroup='smoothed',
                        showlegend=(gene_idx == 0),
                        hoverinfo='skip'
                    ),
                    row=row, col=1
                )
                
            except Exception as e:
                print(f"Warning: Could not fit smooth curve for {gene}: {e}")
        
        # Update y-axis label
        y_title = f'{gene} Expression'
        if transformation == 'log':
            y_title += ' (log)'
        elif transformation == 'z_score':
            y_title += ' (z-score)'
        
        fig.update_yaxes(title_text=y_title, row=row, col=1)
    
    # Update layout
    fig.update_xaxes(title_text='Pseudotime', row=n_genes, col=1)
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=row_height * n_genes,
        margin=dict(t=50, b=50, l=80, r=150),
        legend=dict(
            orientation='v',
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='rgba(0, 0, 0, 0.2)',
            borderwidth=1
        ),
        hovermode='closest'
    )
    
    # Update all axes
    fig.update_xaxes(
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor='black',
        zeroline=False
    )
    
    fig.update_yaxes(
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor='black',
        zeroline=False
    )
    
    return fig