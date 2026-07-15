import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from guanaco.data.loader import obs_col


def _resolve_color_discrete_map(color_map, categories):
    """Map each category to a color using the same palette as the other plots.

    ``color_map`` may be a ready dict (category -> color, as resolved by the
    callback from the discrete-colormap dropdown), ``None`` (fall back to
    Plotly's qualitative palette), or a list of colors to cycle through.
    """
    if isinstance(color_map, dict):
        return color_map
    palette = px.colors.qualitative.Plotly if color_map is None else color_map
    return {cat: palette[i % len(palette)] for i, cat in enumerate(categories)}


def plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=None, y_order=None, x_order=None):
    """Plot stacked bar chart."""
    # Check if x_meta and y_meta are the same - if so, create a histogram
    if x_meta == y_meta:
        # Create a simple count dataframe for histogram
        count_df = obs_col(adata.obs, x_meta).value_counts().reset_index()
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
        
        # Color each bar by its category using the same discrete palette as the
        # embedding/other plots, so a category keeps a consistent color across
        # plots instead of collapsing to a single flat blue.
        bar_categories = [str(c) for c in count_df[x_meta].tolist()]
        color_discrete_map = _resolve_color_discrete_map(color_map, bar_categories)

        # Create a simple bar plot (histogram)
        fig = px.bar(
            count_df,
            x=x_meta,
            y=y_value,
            color=x_meta,
            labels={x_meta: f'{x_meta}', y_value: y_label},
            color_discrete_map=color_discrete_map,
        )

    else:
        # Original stacked bar logic
        # Create count dataframe
        count_df = pd.DataFrame({x_meta: obs_col(adata.obs, x_meta), y_meta: obs_col(adata.obs, y_meta)}).groupby([x_meta, y_meta]).size().reset_index(name='count')
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
    
        # Set up color mapping (consistent with the histogram branch above)
        color_discrete_map = _resolve_color_discrete_map(
            color_map, sorted(count_df[y_meta].unique())
        )

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


def plot_composition_hierarchy(
    parent_meta,
    child_meta,
    adata,
    *,
    parent_color_map=None,
    child_color_map=None,
    parent_order=None,
):
    """Plot a two-level composition hierarchy as an icicle chart."""
    def string_values(key):
        values = obs_col(adata.obs, key).astype(object)
        return values.where(values.notna(), "Missing").astype(str)

    frame = pd.DataFrame(
        {
            parent_meta: string_values(parent_meta),
            child_meta: string_values(child_meta),
        }
    )
    total_count = len(frame)

    if total_count == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No cells are available for this hierarchy.",
            x=0.5,
            y=0.5,
            showarrow=False,
        )
        return fig

    observed_parents = frame[parent_meta].drop_duplicates().tolist()
    requested_order = [str(value) for value in parent_order or []]
    ordered_parents = [
        value for value in requested_order if value in observed_parents
    ] + [value for value in observed_parents if value not in requested_order]

    parent_colors = _resolve_color_discrete_map(
        parent_color_map, ordered_parents
    )
    child_categories = frame[child_meta].drop_duplicates().tolist()
    child_colors = _resolve_color_discrete_map(child_color_map, child_categories)

    ids = ["all-cells"]
    labels = ["All cells"]
    parents = [""]
    counts = [total_count]
    parent_counts = [total_count]
    colors = ["#f1f3f5"]

    same_level = parent_meta == child_meta
    for parent_value in ordered_parents:
        parent_frame = frame[frame[parent_meta] == parent_value]
        parent_count = len(parent_frame)
        parent_id = f"{parent_meta}::{parent_value}"
        ids.append(parent_id)
        labels.append(parent_value)
        parents.append("all-cells")
        counts.append(parent_count)
        parent_counts.append(total_count)
        colors.append(parent_colors.get(parent_value, "#adb5bd"))

        if same_level:
            continue

        child_counts = parent_frame[child_meta].value_counts(sort=False)
        for child_value, child_count in child_counts.items():
            ids.append(f"{parent_id}/{child_meta}::{child_value}")
            labels.append(str(child_value))
            parents.append(parent_id)
            counts.append(int(child_count))
            parent_counts.append(parent_count)
            colors.append(child_colors.get(str(child_value), "#adb5bd"))

    customdata = [
        [
            count,
            count / total_count,
            count / containing_parent if containing_parent else 0,
        ]
        for count, containing_parent in zip(counts, parent_counts, strict=True)
    ]

    fig = go.Figure(
        go.Icicle(
            ids=ids,
            labels=labels,
            parents=parents,
            values=counts,
            branchvalues="total",
            sort=False,
            marker={"colors": colors, "line": {"color": "white", "width": 2}},
            customdata=customdata,
            textinfo="label+percent parent",
            hovertemplate=(
                "<b>%{label}</b><br>"
                "Cells: %{value:,}<br>"
                "Of all cells: %{customdata[1]:.1%}<br>"
                "Of parent: %{customdata[2]:.1%}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin={"t": 20, "r": 20, "b": 20, "l": 20},
    )
    return fig


def plot_composition_differential_abundance(
    results,
    *,
    alpha=0.05,
):
    """Render pairwise ALR effects as a population-by-comparison heatmap."""
    populations = results["population"].drop_duplicates().tolist()
    comparisons = results["comparison"].drop_duplicates().tolist()
    indexed = results.set_index(["population", "comparison"])

    effects = []
    significance = []
    customdata = []
    for population in populations:
        effect_row = []
        significance_row = []
        custom_row = []
        for comparison in comparisons:
            row = indexed.loc[(population, comparison)]
            effect_row.append(float(row["effect"]))
            significance_row.append("*" if row["q_value"] < alpha else "")
            custom_row.append(
                [
                    row["group_a"],
                    row["group_b"],
                    row["ci_low"],
                    row["ci_high"],
                    row["p_value"],
                    row["q_value"],
                    row["n_a"],
                    row["n_b"],
                    row["mean_proportion_a"],
                    row["mean_proportion_b"],
                ]
            )
        effects.append(effect_row)
        significance.append(significance_row)
        customdata.append(custom_row)

    max_effect = max(abs(float(value)) for row in effects for value in row)
    color_limit = max(max_effect, 1e-9)
    fig = go.Figure(
        go.Heatmap(
            z=effects,
            x=comparisons,
            y=populations,
            customdata=customdata,
            text=significance,
            texttemplate="%{text}",
            textfont={"size": 18, "color": "black"},
            colorscale="RdBu_r",
            zmid=0,
            zmin=-color_limit,
            zmax=color_limit,
            colorbar={"title": "ALR<br>difference"},
            hovertemplate=(
                "<b>%{y}</b><br>"
                "%{customdata[1]} − %{customdata[0]}<br>"
                "ALR difference: %{z:.3f}<br>"
                "95% CI: [%{customdata[2]:.3f}, %{customdata[3]:.3f}]<br>"
                "p-value: %{customdata[4]:.3g}<br>"
                "FDR: %{customdata[5]:.3g}<br>"
                "%{customdata[0]} mean proportion: %{customdata[8]:.2%}<br>"
                "%{customdata[1]} mean proportion: %{customdata[9]:.2%}<br>"
                "Samples: %{customdata[6]} vs %{customdata[7]}"
                "<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title={
            "text": (
                "ALR + pairwise Welch t-tests"
                f"<br><sup>* FDR &lt; {alpha:g}</sup>"
            ),
            "x": 0.01,
            "xanchor": "left",
        },
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=max(340, min(900, 34 * len(populations) + 190)),
        margin={"t": 85, "r": 90, "b": 90, "l": 130},
        xaxis={"title": "Pairwise x-axis comparison", "automargin": True},
        yaxis={"title": "Population", "automargin": True},
    )
    return fig
