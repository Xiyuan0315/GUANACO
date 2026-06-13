import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px
from PIL import Image
from guanaco.utils.colors import resolve_continuous_colorscale
from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation

EMBEDDING_PREFIXES = {
    "X_umap": "UMAP",
    "X_pca": "PCA",
    "X_tsne": "t-SNE",
    "X_diffmap": "DiffMap",
    "X_phate": "PHATE",
    "X_draw_graph_fa": "FA",
}


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


def _resolve_embedding_coords(
    adata,
    embedding_key,
    x_axis=None,
    y_axis=None,
    *,
    img_key=None,
    library_id=None,
    auto_select_spatial=False,
):
    """Return (x_values, y_values, x_axis, y_axis, embedding_columns, spatial_image)."""
    embedding_data = adata.obsm[embedding_key]
    embedding_prefix = EMBEDDING_PREFIXES.get(embedding_key, embedding_key.removeprefix("X_"))
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

    return x_values, y_values, x_axis, y_axis, embedding_columns, spatial_image


def _resolve_continuous_values(adata, key, *, source_adata=None, cell_indices=None, layer=None):
    if source_adata is None:
        source_adata = adata
    row_idx = None if cell_indices is None else np.asarray(cell_indices, dtype=np.int64)

    # Follow Scanpy-like precedence: obs -> var_names
    if key in adata.obs.columns:
        values = pd.to_numeric(adata.obs[key], errors="coerce").to_numpy()
    elif key in source_adata.var_names:
        values = extract_gene_expression(source_adata, key, layer=layer)
        if row_idx is not None:
            values = values[row_idx]
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

    return apply_transformation(values, transformation, copy=False)


# Cap the displayed tissue-image resolution. The browser re-rasterizes the
# background image on every zoom/pan frame and re-encodes it on each figure
# rebuild; a multi-thousand-pixel hires/fullres image is far larger than the plot
# can ever show, so it makes spatial views noticeably slower than a plain WebGL
# embedding for no visible benefit. ~1600 px on the longest side keeps it crisp
# (including a few zoom-in steps) while cutting the per-frame raster + payload cost.
SPATIAL_IMAGE_MAX_DIM = 1600


def _embedding_scatter_trace():
    """Trace class used for point embeddings.

    Spatial backgrounds are Plotly layout images with ``layer="below"``. In the
    supported Plotly 5.x range, Scattergl composites above that image while
    preserving transparent gaps between points.
    """
    return go.Scattergl


def _categorical_customdata(obs_names, color_values, indices):
    """Return [cell id, label, row position] for selection and cross-highlight."""
    return np.column_stack([obs_names[indices], color_values[indices], indices])


def _order_continuous_points(x_values, y_values, color_values, cell_idx, highlighted_mask, order):
    if order == "max":
        order_idx = np.argsort(color_values, kind="mergesort")
    elif order == "min":
        order_idx = np.argsort(color_values, kind="mergesort")[::-1]
    elif order == "random":
        rng = np.random.default_rng(315)
        order_idx = rng.permutation(color_values.size)
    else:
        return x_values, y_values, color_values, cell_idx, highlighted_mask

    x_values = x_values[order_idx]
    y_values = y_values[order_idx]
    color_values = color_values[order_idx]
    cell_idx = cell_idx[order_idx]
    if highlighted_mask is not None:
        highlighted_mask = highlighted_mask[order_idx]
    return x_values, y_values, color_values, cell_idx, highlighted_mask


def _apply_spatial_background(fig, spatial_image, img_alpha=1.0):
    if spatial_image is None:
        return

    if isinstance(spatial_image, np.ndarray):
        img_array = spatial_image
        if img_array.dtype != np.uint8:
            img_array = np.clip(img_array, 0, 1)
            img_array = (img_array * 255).astype(np.uint8)
        spatial_image = Image.fromarray(img_array)

    # Original pixel dims define the coordinate space (coords were scaled into it),
    # so sizex/sizey and the axis ranges must stay at the original size.
    img_w, img_h = spatial_image.size

    # Downscale only the *bitmap* shown as the background; keep sizex/sizey at the
    # original dims so the points stay aligned. This shrinks the per-zoom raster and
    # the figure payload without moving anything.
    longest = max(img_w, img_h)
    display_image = spatial_image
    if longest > SPATIAL_IMAGE_MAX_DIM:
        scale = SPATIAL_IMAGE_MAX_DIM / longest
        display_image = spatial_image.resize(
            (max(1, int(round(img_w * scale))), max(1, int(round(img_h * scale)))),
            Image.BILINEAR,
        )

    fig.update_layout(
        images=[dict(
            source=display_image,
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


def _datashader_canvas(ds, x_vals, y_vals):
    """Finite data bounds + a proportionally-sized datashader Canvas.

    ``x_vals``/``y_vals`` must already be finite float arrays. Returns
    ``(canvas, x_min, x_max, y_min, y_max, x_span, y_span)``, or ``None`` if the
    bounds are not finite.
    """
    x_min, x_max = float(np.nanmin(x_vals)), float(np.nanmax(x_vals))
    y_min, y_max = float(np.nanmin(y_vals)), float(np.nanmax(y_vals))
    if not all(np.isfinite(v) for v in (x_min, x_max, y_min, y_max)):
        return None
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

    canvas = ds.Canvas(
        plot_width=int(plot_width),
        plot_height=int(plot_height),
        x_range=(x_min, x_max),
        y_range=(y_min, y_max),
    )
    return canvas, x_min, x_max, y_min, y_max, x_span, y_span


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

    result = _datashader_canvas(ds, x_vals, y_vals)
    if result is None:
        return None, "invalid axis bounds"
    canvas = result[0]

    ds_df = pd.DataFrame({"x": x_vals, "y": y_vals, "v": c_vals})
    agg = canvas.points(ds_df, "x", "y", agg=ds.mean("v"))
    z = np.asarray(agg.values, dtype=np.float32)
    if z.size == 0 or np.isnan(z).all():
        return None, "datashader aggregation produced empty grid"

    x_grid = np.asarray(agg.coords["x"].values, dtype=np.float32)
    y_grid = np.asarray(agg.coords["y"].values, dtype=np.float32)
    cmin = float(np.nanmin(c_vals))
    cmax = float(np.nanmax(c_vals))
    colorscale = resolve_continuous_colorscale(color_map)

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


def _build_datashader_categorical_figure(
    x_values,
    y_values,
    cat_values,
    *,
    label_to_color_dict,
    value_key,
    x_axis,
    y_axis,
    axis_show,
    legend_show="on legend",
    marker_size=5,
):
    """Server-side categorical rasterization via datashader.

    Aggregates per category with ``ds.count_cat`` and blends them into a single
    RGBA image with ``tf.shade`` (using the same ``label_to_color_dict`` as the
    scattergl path so colors match). The image is embedded at data coordinates,
    so the figure payload is constant regardless of cell count. Zero-point dummy
    traces provide the legend. Returns ``(fig, None)`` on success or
    ``(None, reason)`` so the caller can fall back to scattergl.

    Note: a rasterized layer has no per-point hover or lasso, so the caller only
    uses this when nothing is highlighted/selected.
    """
    try:
        import datashader as ds
        import datashader.transfer_functions as tf
    except Exception as exc:
        return None, f"datashader unavailable: {exc}"

    valid_mask = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(valid_mask):
        return None, "all coordinates are NaN/inf after filtering"

    x_vals = np.asarray(x_values[valid_mask], dtype=np.float32)
    y_vals = np.asarray(y_values[valid_mask], dtype=np.float32)
    cat_strs = np.asarray([str(c) for c in cat_values[valid_mask]], dtype=object)

    result = _datashader_canvas(ds, x_vals, y_vals)
    if result is None:
        return None, "invalid axis bounds"
    canvas, x_min, x_max, y_min, y_max, x_span, y_span = result

    df = pd.DataFrame({"x": x_vals, "y": y_vals, "cat": pd.Categorical(cat_strs)})
    present_categories = list(df["cat"].cat.categories)
    if not present_categories:
        return None, "no categories to render"

    # Align the color key with the scattergl mapping; grey for anything missing.
    color_key = {
        cat: label_to_color_dict.get(cat, "#808080") for cat in present_categories
    }

    agg = canvas.points(df, "x", "y", agg=ds.count_cat("cat"))
    if int(agg.sum()) == 0:
        return None, "datashader aggregation produced empty grid"

    shaded = tf.shade(agg, color_key=color_key, how="eq_hist")
    pil_img = shaded.to_pil()  # top row = max y, ready for layout-image placement

    fig = go.Figure()
    # The rasterized embedding as a single background image at data coordinates.
    fig.update_layout(
        images=[dict(
            source=pil_img,
            xref="x", yref="y",
            x=x_min, y=y_max,
            sizex=x_span, sizey=y_span,
            xanchor="left", yanchor="top",
            sizing="stretch",
            layer="below",
        )]
    )

    on_data = legend_show == "on data"
    # Zero-point dummy traces give a categorical legend without any per-cell data.
    if not on_data:
        for cat in present_categories:
            fig.add_trace(go.Scattergl(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=marker_size, color=color_key[cat]),
                name=str(cat),
                showlegend=True,
                hoverinfo="skip",
            ))
    else:
        # Label each category at its centroid directly on the plot.
        for cat in present_categories:
            sel = cat_strs == cat
            if not np.any(sel):
                continue
            fig.add_annotation(
                x=float(np.median(x_vals[sel])),
                y=float(np.median(y_vals[sel])),
                text=f"<b>{cat}</b>",
                showarrow=False,
                font=dict(size=12, color="black"),
                xanchor="center",
                yanchor="middle",
                opacity=0.9,
            )

    _apply_embedding_layout(
        fig,
        title_text=value_key,
        x_axis=x_axis,
        y_axis=y_axis,
        axis_show=axis_show,
        margin=dict(t=40, r=20, l=50, b=50),
        # The raster bakes the cells into an image, so the legend can't hide/show
        # categories. Make it a non-interactive color reference: a legend click
        # would otherwise fire restyleData -> the highlight cascade, which forces
        # the linked plot into an all-grey scattergl. Disabling itemclick stops it.
        legend=dict(itemclick=False, itemdoubleclick=False),
    )
    # Pin axes to the data bounds so the image lines up with the raster.
    fig.update_xaxes(range=[x_min, x_max])
    fig.update_yaxes(range=[y_min, y_max])
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
    highlighted_mask=None,
):
    datashader_err = None
    has_highlight = highlighted_mask is not None
    if render_backend == "datashader" and not has_highlight:
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
    colorscale = resolve_continuous_colorscale(color_map)

    ScatterTrace = _embedding_scatter_trace()

    fig = go.Figure()
    if has_highlight:
        highlighted_mask = np.asarray(highlighted_mask, dtype=bool)
        fig.add_trace(
            ScatterTrace(
                x=x_values,
                y=y_values,
                mode="markers",
                marker=dict(size=marker_size, color="lightgrey", opacity=max(0.15, opacity * 0.3)),
                hoverinfo="skip",
                showlegend=False,
                selectedpoints=None,
            )
        )
        plot_mask = highlighted_mask
    else:
        plot_mask = np.ones(color_values.size, dtype=bool)

    fig.add_trace(
        ScatterTrace(
            x=x_values[plot_mask],
            y=y_values[plot_mask],
            mode="markers",
            marker=dict(
                color=color_values[plot_mask],
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax,
                size=marker_size,
                opacity=opacity,
                colorbar=dict(title=colorbar_title, len=0.8),
            ),
            customdata=cell_idx[plot_mask],
            hoverinfo="skip",
            selectedpoints=None,
            selected=dict(marker=dict(opacity=1)),
            unselected=dict(marker=dict(color="lightgrey", opacity=0.5)),
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


def _color_is_continuous(adata, color, *, threshold=50):
    """Decide whether ``color`` should be drawn with a continuous colormap.

    A gene (in ``var_names``) is always continuous. An ``obs`` column is
    continuous only when it is numeric with many distinct values (e.g.
    ``n_counts`` / ``n_genes``); a low-cardinality numeric (cluster ids stored as
    ints) or any non-numeric column stays categorical. Mirrors the web app's
    ``is_continuous_annotation`` so the notebook API and app agree -- without it,
    a continuous obs like ``n_count_all`` gets thousands of discrete colors,
    which is wrong and very slow.
    """
    if color in adata.var_names:
        return True
    if color in adata.obs.columns:
        dtype = adata.obs[color].dtype
        if getattr(dtype, "kind", None) in ("i", "u", "f"):
            return adata.obs[color].nunique() >= threshold
    return False


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
    layer=None,
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
    source_adata=None,
    cell_indices=None,
    highlighted_cell_ids=None,
    show_background_layer=False,
):
    """
    Unified embedding plotter for both continuous and categorical coloring.

    ``show_background_layer`` (categorical only) always draws a grey "Background"
    layer of every cell beneath the colored groups. When the user hides a group via
    the legend, its cells show through as grey instead of vanishing, and because the
    background spans all cells the axis range stays put -- all handled client-side by
    the native legend toggle, with no server round-trip.
    """
    if mode == "auto":
        mode = "continuous" if _color_is_continuous(adata, color) else "categorical"
    if mode not in {"continuous", "categorical"}:
        raise ValueError(f"Unsupported mode '{mode}'. Expected 'continuous' or 'categorical'.")

    auto_select_spatial = True
    x_values, y_values, x_axis, y_axis, _, spatial_image = _resolve_embedding_coords(
        adata,
        embedding_key,
        x_axis=x_axis,
        y_axis=y_axis,
        img_key=img_key,
        library_id=library_id,
        auto_select_spatial=auto_select_spatial,
    )
    obs_names = np.asarray(adata.obs_names, dtype=object)
    highlighted_mask = None
    if highlighted_cell_ids is not None:
        highlighted_mask = np.isin(obs_names, highlighted_cell_ids)

    if mode == "continuous":
        values = _resolve_continuous_values(
            adata,
            color,
            source_adata=source_adata if source_adata is not None else adata,
            cell_indices=cell_indices,
            layer=layer,
        )
        values = _transform_continuous_values(values, transformation)
        color_values = np.asarray(values, dtype=np.float32)
        cell_idx = np.arange(color_values.size, dtype=np.int32)
        x_values, y_values, color_values, cell_idx, highlighted_mask = _order_continuous_points(
            x_values,
            y_values,
            color_values,
            cell_idx,
            highlighted_mask,
            order,
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
            highlighted_mask=highlighted_mask,
        )

    adata_full = adata if adata_full is None else adata_full
    all_unique_labels = sorted(adata_full.obs[color].unique())
    palette = discrete_color_map or px.colors.qualitative.Plotly
    label_to_color_dict = {
        label: palette[i % len(palette)]
        for i, label in enumerate(all_unique_labels)
    }
    on_data = legend_show == "on data"

    color_values = adata.obs[color].to_numpy()

    # Server-side rasterization for large datasets: only when nothing is
    # highlighted (the raster has no per-point lasso/hover) and there's no
    # spatial background to compose with. Falls back to scattergl on any issue.
    if render_backend == "datashader" and highlighted_mask is None and spatial_image is None:
        ds_fig, _ds_err = _build_datashader_categorical_figure(
            x_values,
            y_values,
            color_values,
            label_to_color_dict={str(k): v for k, v in label_to_color_dict.items()},
            value_key=color,
            x_axis=x_axis,
            y_axis=y_axis,
            axis_show=axis_show,
            legend_show=legend_show,
            marker_size=marker_size,
        )
        if ds_fig is not None:
            return ds_fig

    unique_labels_filtered = sorted(pd.unique(color_values))
    if highlighted_mask is None:
        active_indices = np.arange(color_values.size, dtype=np.int64)
        active_labels = color_values
    else:
        active_indices = np.flatnonzero(highlighted_mask)
        active_labels = color_values[active_indices]
    label_to_indices = {}
    for label in unique_labels_filtered:
        idx = active_indices[active_labels == label]
        if idx.size:
            label_to_indices[label] = idx.tolist()

    ScatterTrace = _embedding_scatter_trace()

    fig = go.Figure()

    # The grey "Background" trace only serves to show cells that aren't in a
    # colored trace. Without a highlight, the colored per-label traces already
    # cover every cell, so the background is fully overdrawn and invisible --
    # skip it. Drawing it would serialize a second full set of per-cell arrays
    # (x/y + all obs_names), which is a major contributor to the first-page-load
    # memory peak in backed mode (all cells served). When highlighting, the
    # background is needed to grey out non-highlighted cells; it carries no
    # customdata since hover is skipped and selection targets the colored traces.
    if highlighted_mask is not None or show_background_layer:
        fig.add_trace(ScatterTrace(
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
        fig.add_trace(ScatterTrace(
            x=x_values[indices],
            y=y_values[indices],
            mode="markers",
            marker=dict(
                size=marker_size,
                color=label_to_color_dict[label],
                opacity=opacity,
            ),
            name=str(label),
            customdata=_categorical_customdata(obs_names, color_values, indices),
            hovertemplate=f"{color}: %{{customdata[1]}}<extra></extra>",
            showlegend=not on_data,
            legendgroup=str(label),
            selectedpoints=None,
            selected=dict(marker=dict(opacity=opacity)),
            unselected=dict(marker=dict(color="lightgrey", opacity=0.5)),
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

    return fig


def plot_coexpression_embedding(
    adata, embedding_key, gene1, gene2,
    x_axis=None, y_axis=None,
    threshold1=0.5, threshold2=0.5,
    transformation=None,
    layer=None,
    color_map=None,
    marker_size=5, opacity=1,
    legend_show='right', axis_show=True,
    img_key=None, library_id=None, img_alpha=1.0,
    source_adata=None, cell_indices=None, highlighted_cell_ids=None,
):
    """
    Plot co-expression of two genes on a 2D embedding.
    Cells are categorized into 4 groups based on expression thresholds.
    """
    # auto_select_spatial=True so a "spatial" basis resolves its tissue image (and
    # scales the coords into the image's pixel space); without it the image is
    # dropped and the points float on a white background.
    x_values, y_values, x_axis, y_axis, _, spatial_image = _resolve_embedding_coords(
        adata, embedding_key, x_axis=x_axis, y_axis=y_axis,
        img_key=img_key, library_id=library_id, auto_select_spatial=True,
    )

    if source_adata is None:
        source_adata = adata
    row_idx = None if cell_indices is None else np.asarray(cell_indices, dtype=np.int64)

    # Extract gene expression for both genes from source adata, then subset rows.
    gene1_expr = extract_gene_expression(source_adata, gene1, layer=layer)
    gene2_expr = extract_gene_expression(source_adata, gene2, layer=layer)
    if row_idx is not None:
        gene1_expr = gene1_expr[row_idx]
        gene2_expr = gene2_expr[row_idx]
    
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
    
    all_cell_idx = np.arange(len(adata), dtype=np.int32)
    obs_names = np.asarray(adata.obs_names, dtype=object)
    highlighted_mask = None
    if highlighted_cell_ids is not None:
        highlighted_mask = np.isin(obs_names, highlighted_cell_ids)

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
    
    ScatterTrace = _embedding_scatter_trace()

    fig = go.Figure()
    if highlighted_mask is not None:
        fig.add_trace(ScatterTrace(
            x=x_values,
            y=y_values,
            mode='markers',
            marker=dict(
                size=marker_size,
                color='lightgrey',
                opacity=max(0.15, opacity * 0.3),
            ),
            name='Background',
            hoverinfo='skip',
            showlegend=False,
        ))

    # Add traces for each category
    for category in category_order:
        mask = categories == category
        if highlighted_mask is not None:
            mask &= highlighted_mask
        if mask.any():
            fig.add_trace(ScatterTrace(
                x=x_values[mask],
                y=y_values[mask],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=color_map.get(category, '#808080'),
                    opacity=opacity,
                ),
                name=category,
                customdata=all_cell_idx[mask],
                showlegend=(legend_show == 'right'),
                hoverinfo='skip',
                selectedpoints=None,
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(color="lightgrey", opacity=0.5))
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

    _apply_embedding_layout(
        fig,
        title_text=f"Co-expression: {gene1} & {gene2}",
        x_axis=x_axis,
        y_axis=y_axis,
        axis_show=axis_show,
        margin=dict(t=40, r=100, l=50, b=50),
        legend=dict(
            orientation='v',
            itemsizing='constant',
            x=1.02, y=0.5,
            bgcolor='rgba(0,0,0,0)',
            itemclick='toggle',
            itemdoubleclick='toggleothers',
            font=dict(size=10),
        ) if legend_show == 'right' else None,
    )

    # Draw the tissue image beneath the points for a spatial basis (matches
    # plot_embedding). Done last so it sets the final axis ranges/aspect.
    if embedding_key == "spatial":
        _apply_spatial_background(fig, spatial_image, img_alpha=img_alpha)

    return fig
