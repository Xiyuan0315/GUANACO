import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from guanaco.data_loader import color_config

def plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=None, y_order=None, x_order=None):
    """Plot stacked bar chart."""
    # Check if x_meta and y_meta are the same - if so, create a histogram
    if x_meta == y_meta:
        # Create a simple count dataframe for histogram
        count_df = adata.obs[x_meta].value_counts().reset_index()
        count_df.columns = [x_meta, 'count']
        count_df[x_meta] = count_df[x_meta].astype(str)
        
        # For histogram, proportion doesn't make sense (always 1.0)
        if norm == 'prop':
            # Calculate proportion of total cells
            count_df['prop'] = count_df['count'] / count_df['count'].sum()
            y_value = 'prop'
            y_label = 'Proportion of Total Cells'
        else:
            y_value = 'count'
            y_label = 'Cell Count'
        
        # Apply filtering based on y_order (from Select Labels) since x_meta == y_meta
        if y_order is not None and len(y_order) > 0:
            y_order_str = [str(y) for y in y_order]
            count_df = count_df[count_df[x_meta].isin(y_order_str)]
            count_df[x_meta] = pd.Categorical(count_df[x_meta], categories=y_order_str, ordered=True)
            count_df = count_df.sort_values(x_meta)
        # Also check x_order in case user reordered in the UI
        elif x_order is not None and len(x_order) > 0:
            x_order_str = [str(x) for x in x_order]
            count_df = count_df[count_df[x_meta].isin(x_order_str)]
            count_df[x_meta] = pd.Categorical(count_df[x_meta], categories=x_order_str, ordered=True)
            count_df = count_df.sort_values(x_meta)
        
        # Create a simple bar plot (histogram)
        fig = px.bar(
            count_df,
            x=x_meta,
            y=y_value,
            labels={x_meta: f'{x_meta}', y_value: y_label},
            color_discrete_sequence=['#1f77b4']  # Single color for histogram
        )
        
        # Add value labels on bars
        fig.update_traces(
            texttemplate='%{y:.2f}' if norm == 'prop' else '%{y}',
            textposition='outside'
        )
        
    else:
        # Original stacked bar logic
        # Create count dataframe
        count_df = adata.obs.groupby([x_meta, y_meta]).size().reset_index(name='count')
        count_df[x_meta] = count_df[x_meta].astype(str)
        count_df[y_meta] = count_df[y_meta].astype(str)

        # Calculate proportions if needed
        if norm == 'prop':
            count_df['prop'] = count_df.groupby(x_meta)['count'].transform(lambda x: x / x.sum())
            y_value = 'prop'
            y_label = 'Proportion'
        else:
            y_value = 'count'
            y_label = 'Cell Count'
    
        # Set up color mapping
        if isinstance(color_map, dict):
            color_discrete_map = color_map
        elif color_map is None:
            categories = sorted(count_df[y_meta].unique())
            predefined_colors = color_config
            color_discrete_map = {cat: predefined_colors[i % len(predefined_colors)] for i, cat in enumerate(categories)}
        else:
            categories = sorted(count_df[y_meta].unique())
            color_discrete_map = {cat: color_map[i % len(color_map)] for i, cat in enumerate(categories)}

        # Apply x_order if specified (this comes from the draggable dropdown)
        if x_order is not None and len(x_order) > 0:
            # Convert x_order to strings for consistency
            x_order_str = [str(x) for x in x_order]
            # Filter the dataframe to only include x groups in x_order
            count_df = count_df[count_df[x_meta].isin(x_order_str)]
        
        # Apply y_order if specified (this comes from Select Labels)
        if y_order is not None and len(y_order) > 0:
            # Convert y_order to strings for consistency
            y_order_str = [str(y) for y in y_order]
            # Filter the dataframe to only include categories in y_order
            count_df = count_df[count_df[y_meta].isin(y_order_str)]
            # Sort the dataframe to match the order in y_order
            count_df[y_meta] = pd.Categorical(count_df[y_meta], categories=y_order_str, ordered=True)
            count_df = count_df.sort_values([x_meta, y_meta])
            category_orders = {y_meta: y_order_str}
        else:
            # If no order specified, use all categories in sorted order
            category_orders = None
        
        # Add x-axis ordering to category_orders
        if x_order is not None and len(x_order) > 0:
            if category_orders is None:
                category_orders = {}
            category_orders[x_meta] = x_order_str

        # Create the plot
        fig = px.bar(
            count_df,
            x=x_meta,
            y=y_value,
            color=y_meta,
            labels={x_meta: f'{x_meta}', y_value: y_label, y_meta: f'{y_meta}'},
            barmode='stack',
            color_discrete_map=color_discrete_map,
            category_orders=category_orders  # This ensures the stacking order follows y_order
        )

    # Update layout
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=False, 
            title_font=dict(size=18)
        ),
        yaxis=dict(showgrid=False, title_font=dict(size=18)),
        legend=dict(
            title=y_meta,
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02
        ),
        margin=dict(r=150)
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

    return fig