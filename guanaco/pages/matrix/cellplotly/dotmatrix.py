import numpy as np
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
import pandas as pd
from .gene_extraction_utils import extract_gene_expression, apply_transformation

def plot_dot_matrix(adata, genes, groupby, selected_labels, aggregation='mean', transformation=None, standardization=None, vmin=None, vmax=None, expression_threshold=0, color_map='Viridis', plot_type='dotplot'):
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
    
    if selected_labels:
        # Keep the order of selected_labels but only include those actually present
        groups = [label for label in reversed(selected_labels) if label in aggregated_data.index]
    else:
        groups = list(reversed(aggregated_data.index))
    
    vmin = vmin if vmin is not None else aggregated_data[valid_genes].min().min()
    vmax = vmax if vmax is not None else aggregated_data[valid_genes].max().max()

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
            (17, [1.0, 0.75, 0.5, 0.25])
        )
        scale, size_legend_values = scale_params

        marker_sizes = df_merged['fraction'] * scale
        size_legend_sizes = [s * scale for s in size_legend_values]
        custom_data = np.stack([df_merged['expression'], df_merged['fraction']], axis=-1)

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_merged['gene'],
            y=df_merged[groupby],
            mode='markers',
            showlegend=False,
            marker=dict(
                size=marker_sizes,
                color=df_merged['expression'],
                colorscale=color_map,
                cmin=vmin,
                cmax=vmax,
                colorbar=dict(
                    title=f'{aggregation.capitalize()} Expression ({transformation})' if transformation and transformation != 'None' else f'{aggregation.capitalize()} Expression',
                    tickfont=dict(color='DarkSlateGrey', size=10),
                    len=0.5,
                    yanchor="bottom",
                    y=0
                )
            ),
            line=dict(color='black', width=0.5),
            customdata=custom_data,
            hovertemplate=(
                'Gene: %{x}<br>'
                f'{groupby}: %{{y}}<br>'
                'Expression: %{customdata[0]:.4f}<br>'
                'Fraction: %{customdata[1]:.4f}<extra></extra>'
            )
        ))

        # Add fraction legend title
        fig.add_trace(go.Scatter(
            x=[None], y=[""],
            mode="text",
            text=["Fraction of Cells"],
            textposition="top center",
            textfont=dict(size=10, color="black"),
            showlegend=False
        ))

        # Add size legend
        for size, value in zip(size_legend_sizes, size_legend_values):
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=size, color="grey", opacity=1),
                name=f"{'Percent Expressed <br>' if size == size_legend_sizes[0] else ''}{value:.2f}",
                legendgroup="size_legend",
                showlegend=True
            ))

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False),
            yaxis=dict(showline=True, linewidth=2, linecolor='black', showgrid=False, categoryorder='array', categoryarray=groups)
        )

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
                    len=0.5,
                    yanchor="bottom",
                    y=0
                ),
            hovertemplate='%{y}<br>%{x}<br>Expression: %{z:.2f}<extra></extra>'
        ))

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(b=100, l=100, t=50, r=50)
        )

    return fig
