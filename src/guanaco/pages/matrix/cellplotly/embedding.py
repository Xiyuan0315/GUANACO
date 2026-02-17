import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px
from PIL import Image
from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation

EMBEDDING_PREFIXES = {
    "X_umap": "UMAP",
    "X_pca": "PCA",
    "X_tsne": "t-SNE",
    "X_diffmap": "DiffMap",
    "X_phate": "PHATE",
    "X_draw_graph_fa": "FA",
}


def _resolve_continuous_colorscale(color_map):
    """
    Accept Plotly colormap names or `cc:<name>` for colorcet palettes.
    Returns a Plotly-compatible colorscale object.
    """
    if not isinstance(color_map, str):
        return color_map

    if not color_map.startswith("cc:"):
        return color_map

    cc_name = color_map.split(":", 1)[1].strip()
    if not cc_name:
        return "Viridis"

    try:
        import colorcet as cc
    except Exception:
        return "Viridis"

    palette = cc.palette.get(cc_name)
    if palette is None:
        return "Viridis"

    colors = list(palette)
    if len(colors) < 2:
        return "Viridis"

    denom = len(colors) - 1
    return [[i / denom, c] for i, c in enumerate(colors)]


def _resolve_spatial_context(
    adata, embedding_key, img_key=None, library_id=None, auto_select=False
):
    if embedding_key != "spatial" or "spatial" not in adata.uns:
        return None, 1.0

    spatial = adata.uns.get("spatial", {})
    if not spatial:
        return None, 1.0

    if library_id is None:
        if len(spatial) == 1:
            library_id = next(iter(spatial.keys()))
        elif "library_id" in adata.obs.columns:
            obs_library_ids = [str(x) for x in pd.Series(adata.obs["library_id"]).dropna().unique()]
            matching = [lib for lib in obs_library_ids if lib in spatial]
            if len(matching) == 1:
                library_id = matching[0]
        if library_id is None:
            if auto_select:
                library_id = sorted(spatial.keys())[0]
            else:
                raise ValueError("Multiple spatial libraries found; please provide library_id.")

    if library_id not in spatial:
        raise KeyError(f"library_id '{library_id}' not found in adata.uns['spatial'].")

    spatial_lib = spatial[library_id]
    images = spatial_lib.get("images", {})
    if img_key is None and auto_select and images:
        preferred_keys = [
            "hires",
            "highres",
            "lowres",
            "fullres",
            "full_res",
            "image",
            "img",
        ]
        for preferred in preferred_keys:
            if preferred in images:
                img_key = preferred
                break
        if img_key is None:
            img_key = next(iter(images.keys()))

    if img_key is None:
        return None, 1.0

    spatial_image = images.get(img_key)
    if spatial_image is None:
        raise KeyError(
            f"img_key '{img_key}' not found in adata.uns['spatial'][{library_id!r}]['images']."
        )

    scalefactors = spatial_lib.get("scalefactors", {})
    candidate_scale_keys = [
        f"tissue_{img_key}_scalef",
        f"{img_key}_scalef",
        "tissue_hires_scalef",
        "tissue_lowres_scalef",
    ]
    spatial_scale = 1.0
    for key in candidate_scale_keys:
        if key in scalefactors:
            spatial_scale = scalefactors[key]
            break
    return spatial_image, spatial_scale


def _prepare_embedding_dataframe(
    adata,
    embedding_key,
    x_axis=None,
    y_axis=None,
    *,
    img_key=None,
    library_id=None,
    auto_select_spatial=False,
):
    embedding_data = adata.obsm[embedding_key]
    embedding_prefix = EMBEDDING_PREFIXES.get(embedding_key, embedding_key.upper())
    embedding_columns = [f"{embedding_prefix}{i + 1}" for i in range(embedding_data.shape[1])]
    col_to_idx = {name: i for i, name in enumerate(embedding_columns)}

    spatial_image, spatial_scale = _resolve_spatial_context(
        adata,
        embedding_key,
        img_key=img_key,
        library_id=library_id,
        auto_select=auto_select_spatial,
    )
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])

    if x_axis not in col_to_idx or y_axis not in col_to_idx:
        raise ValueError(f"Invalid x/y axis for {embedding_key}: x={x_axis}, y={y_axis}")

    x_values = np.asarray(embedding_data[:, col_to_idx[x_axis]], dtype=np.float32)
    if y_axis == x_axis:
        y_values = x_values.copy()
    else:
        y_values = np.asarray(embedding_data[:, col_to_idx[y_axis]], dtype=np.float32)

    if spatial_image is not None and embedding_key == "spatial":
        if x_axis == embedding_columns[0] or x_axis == embedding_columns[1]:
            x_values = x_values * spatial_scale
        if y_axis == embedding_columns[0] or y_axis == embedding_columns[1]:
            y_values = y_values * spatial_scale

    embedding_df = pd.DataFrame({x_axis: x_values, y_axis: y_values})
    return embedding_df, x_axis, y_axis, embedding_columns, spatial_image


def _resolve_continuous_values(adata, key):
    # Follow Scanpy-like precedence: obs -> var_names
    if key in adata.obs.columns:
        values = pd.to_numeric(adata.obs[key], errors="coerce").to_numpy()
    elif key in adata.var_names:
        values = extract_gene_expression(adata, key)
    else:
        raise ValueError(
            f"Invalid key '{key}'. Expected an obs column or a gene in var_names."
        )

    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        raise ValueError(f"Key '{key}' resolved to an empty vector.")
    if np.isnan(values).all():
        raise ValueError(f"Key '{key}' does not contain numeric values for continuous plotting.")
    return values


def _transform_continuous_values(values, transformation):
    if not transformation:
        return values

    if transformation in {"log", "log1p"}:
        return np.log1p(values)

    return apply_transformation(values, transformation, copy=False)


def _apply_spatial_background(fig, spatial_image, img_alpha=1.0):
    if spatial_image is None:
        return

    if isinstance(spatial_image, np.ndarray):
        img_array = spatial_image
        if img_array.dtype != np.uint8:
            img_array = np.clip(img_array, 0, 1)
            img_array = (img_array * 255).astype(np.uint8)
        spatial_image = Image.fromarray(img_array)
        img_h, img_w = img_array.shape[0], img_array.shape[1]
    else:
        img_w, img_h = spatial_image.size

    fig.update_layout(
        images=[dict(
            source=spatial_image,
            xref="x", yref="y",
            x=0, y=0,
            sizex=img_w, sizey=img_h,
            xanchor="left", yanchor="top",
            opacity=img_alpha,
            layer="below",
        )]
    )
    fig.update_xaxes(range=[0, img_w])
    fig.update_yaxes(range=[img_h, 0], scaleanchor="x", scaleratio=1)


def _apply_embedding_layout(fig, *, title_text, x_axis, y_axis, axis_show, margin, legend=None):
    tick_color = "black" if axis_show else "rgba(0,0,0,0)"
    layout_kwargs = dict(
        autosize=True,
        plot_bgcolor="white",
        paper_bgcolor="white",
        title=dict(text=f"<b>{title_text}</b>", x=0.5, y=0.98, xanchor="center", yanchor="top"),
        xaxis=dict(title=x_axis, showgrid=False, zeroline=False, tickfont=dict(color=tick_color)),
        yaxis=dict(title=y_axis, showgrid=False, zeroline=False, tickfont=dict(color=tick_color)),
        margin=margin,
    )
    if legend is not None:
        layout_kwargs["legend"] = legend
    fig.update_layout(**layout_kwargs)
    fig.update_xaxes(showline=True, linewidth=2, linecolor="black")
    fig.update_yaxes(showline=True, linewidth=2, linecolor="black")


def _build_datashader_continuous_figure(
    x_values,
    y_values,
    color_values,
    *,
    value_key,
    x_axis,
    y_axis,
    color_map,
    axis_show,
    colorbar_title,
    spatial_image=None,
    img_alpha=1.0,
):
    try:
        import datashader as ds
    except Exception as exc:
        return None, f"datashader unavailable: {exc}"

    valid_mask = np.isfinite(x_values) & np.isfinite(y_values) & np.isfinite(color_values)
    if not np.any(valid_mask):
        return None, "all values are NaN/inf after filtering"

    x_vals = np.asarray(x_values[valid_mask], dtype=np.float32)
    y_vals = np.asarray(y_values[valid_mask], dtype=np.float32)
    c_vals = np.asarray(color_values[valid_mask], dtype=np.float32)

    x_min, x_max = float(np.nanmin(x_vals)), float(np.nanmax(x_vals))
    y_min, y_max = float(np.nanmin(y_vals)), float(np.nanmax(y_vals))
    if not np.isfinite(x_min) or not np.isfinite(x_max) or not np.isfinite(y_min) or not np.isfinite(y_max):
        return None, "invalid axis bounds"
    if x_min == x_max:
        x_max = x_min + 1e-6
    if y_min == y_max:
        y_max = y_min + 1e-6

    x_span = max(x_max - x_min, 1e-6)
    y_span = max(y_max - y_min, 1e-6)
    base = 900
    if x_span >= y_span:
        plot_width = 1000
        plot_height = int(max(350, min(1400, base * (y_span / x_span))))
    else:
        plot_height = 1000
        plot_width = int(max(350, min(1400, base * (x_span / y_span))))

    ds_df = pd.DataFrame({"x": x_vals, "y": y_vals, "v": c_vals})
    canvas = ds.Canvas(
        plot_width=int(plot_width),
        plot_height=int(plot_height),
        x_range=(x_min, x_max),
        y_range=(y_min, y_max),
    )
    agg = canvas.points(ds_df, "x", "y", agg=ds.mean("v"))
    z = np.asarray(agg.values, dtype=np.float32)
    if z.size == 0 or np.isnan(z).all():
        return None, "datashader aggregation produced empty grid"

    x_grid = np.asarray(agg.coords["x"].values, dtype=np.float32)
    y_grid = np.asarray(agg.coords["y"].values, dtype=np.float32)
    cmin = float(np.nanmin(c_vals))
    cmax = float(np.nanmax(c_vals))
    colorscale = _resolve_continuous_colorscale(color_map)

    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=x_grid,
            y=y_grid,
            colorscale=colorscale,
            zmin=cmin,
            zmax=cmax,
            hoverinfo="skip",
            colorbar=dict(title=colorbar_title, len=0.8),
        )
    )

    _apply_embedding_layout(
        fig,
        title_text=value_key,
        x_axis=x_axis,
        y_axis=y_axis,
        axis_show=axis_show,
        margin=dict(t=40, r=20, l=50, b=50),
    )
    _apply_spatial_background(fig, spatial_image, img_alpha=img_alpha)
    return fig, None


def _build_continuous_embedding_figure(
    x_values,
    y_values,
    color_values,
    cell_idx,
    *,
    value_key,
    x_axis,
    y_axis,
    color_map,
    marker_size,
    opacity,
    axis_show,
    colorbar_title,
    spatial_image=None,
    img_alpha=1.0,
    render_backend="scattergl",
):
    datashader_err = None
    if render_backend == "datashader":
        ds_fig, datashader_err = _build_datashader_continuous_figure(
            x_values,
            y_values,
            color_values,
            value_key=value_key,
            x_axis=x_axis,
            y_axis=y_axis,
            color_map=color_map,
            axis_show=axis_show,
            colorbar_title=colorbar_title,
            spatial_image=spatial_image,
            img_alpha=img_alpha,
        )
        if ds_fig is not None:
            return ds_fig

    cmin = float(np.nanmin(color_values))
    cmax = float(np.nanmax(color_values))
    colorscale = _resolve_continuous_colorscale(color_map)

    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=x_values,
            y=y_values,
            mode="markers",
            marker=dict(
                color=color_values,
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax,
                size=marker_size,
                opacity=opacity,
                colorbar=dict(title=colorbar_title, len=0.8),
            ),
            customdata=cell_idx,
            hoverinfo="skip",
            selectedpoints=None,
            selected=dict(marker=dict(opacity=1)),
            unselected=dict(marker=dict(opacity=0.2)),
        )
    )

    _apply_embedding_layout(
        fig,
        title_text=value_key,
        x_axis=x_axis,
        y_axis=y_axis,
        axis_show=axis_show,
        margin=dict(t=40, r=20, l=50, b=50),
    )

    _apply_spatial_background(fig, spatial_image, img_alpha=img_alpha)
    if datashader_err:
        fig.add_annotation(
            x=0,
            y=1.08,
            xref="paper",
            yref="paper",
            text=f"Datashader fallback to scattergl ({datashader_err})",
            showarrow=False,
            font=dict(size=11, color="darkred"),
            xanchor="left",
            yanchor="bottom",
        )

    return fig


def plot_embedding(
    adata,
    embedding_key,
    color,
    x_axis=None,
    y_axis=None,
    *,
    mode="auto",
    adata_full=None,
    transformation=None,
    order=None,
    continuous_color_map="Viridis",
    discrete_color_map=None,
    marker_size=5,
    opacity=1,
    render_backend="scattergl",
    legend_show="on legend",
    annotation=None,
    axis_show=True,
    img_key=None,
    library_id=None,
    img_alpha=1.0,
):
    """
    Unified embedding plotter for both continuous and categorical coloring.
    """
    if mode == "auto":
        mode = "continuous" if color in adata.var_names else "categorical"
    if mode not in {"continuous", "categorical"}:
        raise ValueError(f"Unsupported mode '{mode}'. Expected 'continuous' or 'categorical'.")

    auto_select_spatial = True
    embedding_df, x_axis, y_axis, _, spatial_image = _prepare_embedding_dataframe(
        adata,
        embedding_key,
        x_axis=x_axis,
        y_axis=y_axis,
        img_key=img_key,
        library_id=library_id,
        auto_select_spatial=auto_select_spatial,
    )

    if mode == "continuous":
        values = _resolve_continuous_values(adata, color)
        values = _transform_continuous_values(values, transformation)
        x_values = embedding_df[x_axis].to_numpy(dtype=np.float32, copy=False)
        y_values = embedding_df[y_axis].to_numpy(dtype=np.float32, copy=False)
        color_values = np.asarray(values, dtype=np.float32)
        cell_idx = np.arange(color_values.size, dtype=np.int32)
        if order == "max":
            order_idx = np.argsort(color_values, kind="mergesort")
            x_values, y_values, color_values, cell_idx = (
                x_values[order_idx], y_values[order_idx], color_values[order_idx], cell_idx[order_idx]
            )
        elif order == "min":
            order_idx = np.argsort(color_values, kind="mergesort")[::-1]
            x_values, y_values, color_values, cell_idx = (
                x_values[order_idx], y_values[order_idx], color_values[order_idx], cell_idx[order_idx]
            )
        elif order == "random":
            rng = np.random.default_rng(315)
            order_idx = rng.permutation(color_values.size)
            x_values, y_values, color_values, cell_idx = (
                x_values[order_idx], y_values[order_idx], color_values[order_idx], cell_idx[order_idx]
            )
        colorbar_title = transformation if transformation else color
        return _build_continuous_embedding_figure(
            x_values,
            y_values,
            color_values,
            cell_idx,
            x_axis=x_axis,
            y_axis=y_axis,
            value_key=color,
            color_map=continuous_color_map,
            marker_size=marker_size,
            opacity=opacity,
            axis_show=axis_show,
            colorbar_title=colorbar_title,
            spatial_image=spatial_image,
            img_alpha=img_alpha,
            render_backend=render_backend,
        )

    adata_full = adata if adata_full is None else adata_full
    all_unique_labels = sorted(adata_full.obs[color].unique())
    palette = discrete_color_map or px.colors.qualitative.Plotly
    label_to_color_dict = {
        label: palette[i % len(palette)]
        for i, label in enumerate(all_unique_labels)
    }
    on_data = legend_show == "on data"

    x_values = embedding_df[x_axis].to_numpy(dtype=np.float32, copy=False)
    y_values = embedding_df[y_axis].to_numpy(dtype=np.float32, copy=False)
    color_values = adata.obs[color].to_numpy()
    unique_labels_filtered = sorted(pd.unique(color_values))
    label_to_indices = {}
    for i, label in enumerate(color_values):
        label_to_indices.setdefault(label, []).append(i)

    fig = go.Figure()

    fig.add_trace(go.Scattergl(
        x=x_values,
        y=y_values,
        mode="markers",
        marker=dict(size=marker_size, color="lightgrey", opacity=opacity * 0.3),
        name="Background",
        hoverinfo="skip",
        showlegend=False,
        visible=True,
    ))

    for label in unique_labels_filtered:
        indices = np.asarray(label_to_indices.get(label, []), dtype=np.int32)
        if indices.size == 0:
            continue
        fig.add_trace(go.Scattergl(
            x=x_values[indices],
            y=y_values[indices],
            mode="markers",
            marker=dict(
                size=marker_size,
                color=label_to_color_dict[label],
                opacity=opacity,
            ),
            name=str(label),
            hoverinfo="skip",
            showlegend=not on_data,
            legendgroup=str(label),
        ))

    if on_data:
        for label in unique_labels_filtered:
            indices = np.asarray(label_to_indices.get(label, []), dtype=np.int32)
            if indices.size == 0:
                continue
            fig.add_annotation(
                x=float(np.median(x_values[indices])),
                y=float(np.median(y_values[indices])),
                text=f"<b>{label}</b>",
                showarrow=False,
                font=dict(size=12, color="black"),
                xanchor="center",
                yanchor="middle",
                opacity=0.9,
            )

    _apply_embedding_layout(
        fig,
        title_text=color,
        x_axis=x_axis,
        y_axis=y_axis,
        axis_show=axis_show,
        margin=dict(t=40, r=100, l=50, b=50),
        legend=dict(
            orientation="v",
            itemsizing="constant",
            x=1.02,
            y=0.5,
            bgcolor="rgba(0,0,0,0)",
            itemclick="toggle",
            itemdoubleclick="toggleothers",
            font=dict(size=12),
        ) if not on_data else None,
    )

    if embedding_key == "spatial":
        _apply_spatial_background(fig, spatial_image, img_alpha=img_alpha)
    if render_backend == "datashader":
        fig.add_annotation(
            x=0,
            y=1.08,
            xref="paper",
            yref="paper",
            text="Datashader mode currently applies to continuous embedding only.",
            showarrow=False,
            font=dict(size=11, color="darkred"),
            xanchor="left",
            yanchor="bottom",
        )
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
    embedding_df, x_axis, y_axis, _, _ = _prepare_embedding_dataframe(
        adata, embedding_key, x_axis=x_axis, y_axis=y_axis
    )

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
    categories = np.zeros(len(adata), dtype=object)
    categories[:] = 'Neither'
    categories[gene1_expressed & ~gene2_expressed] = f'{gene1} only'
    categories[~gene1_expressed & gene2_expressed] = f'{gene2} only'
    categories[gene1_expressed & gene2_expressed] = 'Co-expressed'
    
    x_values = embedding_df[x_axis].to_numpy(dtype=np.float32, copy=False)
    y_values = embedding_df[y_axis].to_numpy(dtype=np.float32, copy=False)
    all_cell_idx = np.arange(len(adata), dtype=np.int32)

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
        mask = categories == category
        if mask.any():
            fig.add_trace(go.Scattergl(
                x=x_values[mask],
                y=y_values[mask],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=color_map.get(category, '#808080'),
                    opacity=opacity,
                ),
                name=category,
                customdata=all_cell_idx[mask],  # Add cell indices for selection
                showlegend=(legend_show == 'right'),
                hoverinfo='skip',
                selectedpoints=None,  # Enable selection
                selected=dict(marker=dict(opacity=1)),  # Keep selected points fully visible
                unselected=dict(marker=dict(opacity=0.2))  # Dim unselected points
            ))

    if legend_show == 'on data':
        for category in category_order:
            mask = categories == category
            if mask.any():
                median_x = float(np.median(x_values[mask]))
                median_y = float(np.median(y_values[mask]))
                fig.add_annotation(
                    x=median_x, y=median_y,
                    text=f'<b>{category}</b>',
                    showarrow=False,
                    font=dict(size=10, color='black'),
                    xanchor='center', yanchor='middle',
                    opacity=0.9,
                )

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
