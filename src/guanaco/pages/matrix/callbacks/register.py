import json
import hashlib
import time
from collections import OrderedDict
import warnings
import numpy as np
from dash import dcc, html, Input, Output, State, callback_context, ALL, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from guanaco.pages.matrix.plots.embedding import (
    plot_embedding,
    plot_coexpression_embedding,
)
from guanaco.pages.matrix.plots.heatmap import (
    plot_unified_heatmap,
)
from guanaco.pages.matrix.plots.violin1 import plot_violin1
from guanaco.pages.matrix.plots.violin2 import plot_violin2_new
from guanaco.pages.matrix.plots.stacked_bar import plot_stacked_bar
from guanaco.pages.matrix.plots.dotmatrix import plot_dot_matrix
from guanaco.pages.matrix.plots.pseudotime import plot_genes_in_pseudotime
from guanaco.pages.matrix.plots.paga import build_paga_cytoscape
from guanaco.pages.matrix.plots.volcano import has_volcano_data
from guanaco.pages.matrix.layouts.heatmap_layout import generate_heatmap_layout
from guanaco.pages.matrix.layouts.violin_layout import generate_violin_layout, generate_split_violin_layout
from guanaco.pages.matrix.layouts.dotplot_layout import generate_dotplot_layout
from guanaco.pages.matrix.layouts.stacked_bar_layout import generate_stacked_bar_layout
from guanaco.pages.matrix.layouts.pseudotime_layout import generate_pseudotime_layout
from guanaco.pages.matrix.layouts.paga_layout import generate_paga_layout
from guanaco.pages.matrix.layouts.volcano_layout import generate_volcano_layout
from guanaco.pages.matrix.layouts.grn_demo_layout import generate_grn_demo_layout
from guanaco.pages.matrix.layouts.atac_browser_layout import generate_atac_browser_layout
from guanaco.pages.matrix.layouts.embedding_layout import (
    generate_embedding_plots as build_embedding_plots,
    initialize_scatter_components as build_initialize_scatter_components,
)
from guanaco.pages.matrix.callbacks.scatter_callbacks import register_scatter_callbacks
from guanaco.pages.matrix.callbacks.heatmap_callbacks import register_heatmap_callbacks
from guanaco.pages.matrix.callbacks.dotplot_callbacks import register_dotplot_callbacks
from guanaco.pages.matrix.callbacks.pseudotime_callbacks import register_pseudotime_callbacks
from guanaco.pages.matrix.callbacks.violin_callbacks import register_violin_callbacks
from guanaco.pages.matrix.callbacks.stacked_bar_callbacks import register_stacked_bar_callbacks
from guanaco.pages.matrix.callbacks.paga_callbacks import register_paga_callbacks
from guanaco.pages.matrix.callbacks.volcano_callbacks import register_volcano_callbacks
from guanaco.pages.matrix.callbacks.grn_demo_callbacks import register_grn_demo_callbacks
from guanaco.pages.matrix.callbacks.atac_browser_callbacks import register_atac_browser_callbacks
from guanaco.pages.matrix.plots.atac_browser import has_genomic_peak_features
from guanaco.utils.colors import discrete_palette_config
from guanaco.utils.obs_utils import sorted_categories
from guanaco.data.registry import color_config as _default_color_config
warnings.filterwarnings('ignore', message='.*observed=False.*')

palette_json = discrete_palette_config()


class FigureMemoCache:
    """Small in-process TTL+LRU cache for expensive callback figure payloads."""

    def __init__(self, max_items=24, ttl_seconds=300):
        self.max_items = max_items
        self.ttl_seconds = ttl_seconds
        self._store = OrderedDict()

    def _prune_expired(self):
        now = time.time()
        expired = [k for k, (_, ts) in self._store.items() if now - ts > self.ttl_seconds]
        for k in expired:
            self._store.pop(k, None)

    def get(self, key):
        self._prune_expired()
        item = self._store.get(key)
        if item is None:
            return None
        value, ts = item
        # Refresh LRU position
        self._store.move_to_end(key)
        return value

    def set(self, key, value):
        self._prune_expired()
        self._store[key] = (value, time.time())
        self._store.move_to_end(key)
        while len(self._store) > self.max_items:
            self._store.popitem(last=False)


# Each cached figure dict can embed per-cell arrays (scatter x/y/color/obs_names),
# so a handful of large figures dominates RAM.  24 slots covers all plot types
# across typical multi-tab usage without evicting on every tab switch.
_figure_cache = FigureMemoCache(max_items=24, ttl_seconds=300)


class _FilteredDataCache:
    """Cache AnnData views by their filter parameters.

    AnnData slice creation (adata[mask]) is cheap, but it involves rebuilding obs/var
    DataFrame subsets and rerunning pandas isin() on every call.  For label-based
    filtering of 200 k-cell datasets that adds ~10 ms per trigger; for backed/lazy
    stores the same view object also benefits from OS-level read caching.

    Note: for in-memory sparse .X, each downstream view.X access still materialises
    a row-slice copy -- this cache eliminates only the view-creation and isin()
    overhead, not per-gene extraction cost.  The largest win is for backed/lazy data
    and for label-based filtering paths.
    """

    def __init__(self, max_items=16):
        self.max_items = max_items
        self._store = OrderedDict()

    @staticmethod
    def _make_key(adata, annotation, selected_labels, selected_cells):
        adata_id = id(adata)
        if selected_cells is not None and len(selected_cells) > 0:
            n = len(selected_cells)
            payload = json.dumps(list(selected_cells), separators=(",", ":"), default=str)
            digest = hashlib.md5(payload.encode()).hexdigest()
            return (adata_id, "cells", n, digest)
        if selected_labels and annotation:
            return (adata_id, "labels", annotation, tuple(sorted(str(label) for label in selected_labels)))
        return None

    def get_or_create(self, adata, annotation, selected_labels, selected_cells):
        key = self._make_key(adata, annotation, selected_labels, selected_cells)
        if key is None:
            return adata
        item = self._store.get(key)
        if item is not None:
            self._store.move_to_end(key)
            return item
        if selected_cells is not None and len(selected_cells) > 0:
            filtered = adata[selected_cells]
        else:
            filtered = adata[adata.obs[annotation].isin(selected_labels)]
        self._store[key] = filtered
        self._store.move_to_end(key)
        while len(self._store) > self.max_items:
            self._store.popitem(last=False)
        return filtered


_filter_data_cache = _FilteredDataCache(max_items=16)


def _hash_list_signature(values):
    """Compact signature for potentially large lists."""
    if values is None:
        return None
    if not isinstance(values, (list, tuple)):
        return values
    n = len(values)
    if n == 0:
        return {"len": 0, "hash": None}
    # Hash full content for correctness while keeping key compact.
    payload = json.dumps(list(values), sort_keys=False, default=str, separators=(",", ":"))
    digest = hashlib.md5(payload.encode("utf-8")).hexdigest()
    return {"len": n, "hash": digest}


def _make_cache_key(kind, adata, **kwargs):
    payload = {"kind": kind, "adata_id": id(adata), **kwargs}
    payload_json = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.md5(payload_json.encode("utf-8")).hexdigest()


def _cached_figure_get(cache_key):
    cached = _figure_cache.get(cache_key)
    if cached is None:
        return None
    return go.Figure(cached)


def _cached_figure_set(cache_key, fig):
    _figure_cache.set(cache_key, fig.to_dict())


def apply_relayout(fig, relayout):
    if not relayout:
        return fig

    if all(k in relayout for k in ["xaxis.range[0]", "xaxis.range[1]", "yaxis.range[0]", "yaxis.range[1]"]):
        fig.update_layout(
            xaxis=dict(range=[relayout["xaxis.range[0]"], relayout["xaxis.range[1]"]]),
            yaxis=dict(range=[relayout["yaxis.range[0]"], relayout["yaxis.range[1]"]]),
        )
        return fig

    if "xaxis.range" in relayout and "yaxis.range" in relayout:
        fig.update_layout(
            xaxis=dict(range=relayout["xaxis.range"]),
            yaxis=dict(range=relayout["yaxis.range"]),
        )
        return fig

    # 3) Reset axes (autorange)
    if "xaxis.autorange" in relayout or "autosize" in relayout:

        fig.update_layout(
            xaxis=dict(autorange=True),
            yaxis=dict(autorange=True),
        )
        return fig

    return fig


def generate_embedding_plots(adata, prefix):
    """Build the embedding plot layout for a matrix dataset."""
    return build_embedding_plots(adata, prefix)


def filter_data(adata, annotation, selected_labels, selected_cells=None):
    return _filter_data_cache.get_or_create(adata, annotation, selected_labels, selected_cells)


def generate_left_control(default_gene_markers, label_list, prefix):
    # Radio buttons for input mode selection
    input_mode_radio = dbc.RadioItems(
        id=f'{prefix}-gene-input-mode',
        options=[
            {'label': 'Dropdown', 'value': 'dropdown'},
            {'label': 'Text', 'value': 'text'}
        ],
        value='dropdown',
        inline=True,
        style={'fontSize': '14px', 'marginBottom': '10px'}
    )
    
    # Dropdown for gene selection
    genes_dropdown = dcc.Dropdown(
        id=f'{prefix}-single-cell-genes-selection',
        options=[{'label': gene, 'value': gene} for gene in default_gene_markers],
        value=default_gene_markers,
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )
    
    # Text area for gene input (initially hidden)
    genes_textarea = dcc.Textarea(
        id=f'{prefix}-single-cell-genes-textarea',
        placeholder='Enter genes separated by commas (e.g., Gene1, Gene2, Gene3)',
        value=', '.join(default_gene_markers),
        style={'width': '100%', 'height': '80px', 'marginBottom': '10px', 'display': 'none'},
        className='custom-textarea'
    )
    
    error_message = html.Div(
        id=f'{prefix}-gene-input-error',
        style={'color': 'red', 'fontSize': '12px', 'marginBottom': '10px'}
    )
    
    genes_selection = html.Div([
        input_mode_radio,
        genes_dropdown,
        genes_textarea,
        error_message
    ])
    
    annotation_filter = dcc.Dropdown(
        id=f'{prefix}-single-cell-annotation-dropdown',
        options=[{'label': label, 'value': label} for label in label_list],
        value=label_list[len(label_list) // 2],
        style={'marginBottom': '15px'},
        clearable=False,
        className='custom-dropdown'
    )
    
    label_list_selection = dcc.Dropdown(
        id=f'{prefix}-single-cell-label-selection',
        multi=True,
        style={'marginBottom': '15px', 'font-size': '12px'},
        className='custom-dropdown'
    )

    
    return html.Div([
        html.Label('Select Variables:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        genes_selection,
        html.Label('Select Annotation:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        annotation_filter,
        html.Label('Select Labels:', style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        label_list_selection,
    ])


def generate_other_plots(
    adata,
    default_gene_markers,
    discrete_label_list,
    prefix,
    optional_plot_components=None,
    gene_annotation_path=None,
):
    component_aliases = {
        'vocano': 'volcano',
        'grn-demo': 'grn',
        'expression-trend': 'pseudotime',
        'expression_trend': 'pseudotime',
        'stacked-violin': 'violin',
        'violin2': 'split-violin',
        'split_violin': 'split-violin',
        'splitviolin': 'split-violin',
        'grouped-violin': 'split-violin',
        'group-violin': 'split-violin',
        'trackplot': 'peak-browser',
        'track-plot': 'peak-browser',
        'genome-browser': 'peak-browser',
        'atac_browser': 'peak-browser',
        'atac-browser': 'peak-browser',  # legacy key
    }
    valid_optional_components = {
        'heatmap', 'violin', 'split-violin', 'dotplot', 'stacked-bar',
        'pseudotime', 'paga', 'volcano', 'grn', 'peak-browser',
    }
    if optional_plot_components is None:
        selected_optional_components = ['heatmap', 'violin', 'split-violin', 'dotplot', 'stacked-bar']
        if has_genomic_peak_features(adata):
            selected_optional_components.append('peak-browser')
    else:
        selected_optional_components = [
            component_aliases.get(comp, comp)
            for comp in optional_plot_components
            if component_aliases.get(comp, comp) in valid_optional_components
        ]
    if 'peak-browser' in selected_optional_components and not has_genomic_peak_features(adata):
        selected_optional_components = [comp for comp in selected_optional_components if comp != 'peak-browser']
    if 'paga' in selected_optional_components and not (
        'paga' in adata.uns and 'connectivities' in adata.uns['paga']
    ):
        selected_optional_components = [comp for comp in selected_optional_components if comp != 'paga']
    if 'volcano' in selected_optional_components and not has_volcano_data(adata):
        selected_optional_components = [comp for comp in selected_optional_components if comp != 'volcano']
    tab_order = [
        ('dotplot', 'dotplot-tab'),
        ('heatmap', 'heatmap-tab'),
        ('violin', 'violin-tab'),
        ('split-violin', 'split-violin-tab'),
        ('stacked-bar', 'stacked-bar-tab'),
        ('pseudotime', 'pseudotime-tab'),
        ('paga', 'paga-tab'),
        ('volcano', 'volcano-tab'),
        ('grn', 'grn-tab'),
        ('peak-browser', 'peak-browser-tab'),
    ]
    initial_tab_value = next(
        (tab for key, tab in tab_order if key in selected_optional_components),
        None,
    )

    tabs_children = []
    if 'heatmap' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='Heatmap', value='heatmap-tab', children=[html.Div(generate_heatmap_layout(adata, prefix))])
        )
    if 'violin' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(
                label='Violin Plot',
                value='violin-tab',
                children=[html.Div(generate_violin_layout(default_gene_markers, discrete_label_list, prefix))]
            )
        )
    if 'split-violin' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(
                label='Comparative Violin',
                value='split-violin-tab',
                children=[html.Div(generate_split_violin_layout(default_gene_markers, discrete_label_list, prefix))]
            )
        )
    if 'dotplot' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='Dotplot', value='dotplot-tab', children=[html.Div(generate_dotplot_layout(prefix))])
        )
    if 'stacked-bar' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(
                label='Stacked Bar',
                value='stacked-bar-tab',
                children=[html.Div(generate_stacked_bar_layout(discrete_label_list, prefix))]
            )
        )
    if 'pseudotime' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='Expression Trend', value='pseudotime-tab', children=[html.Div(generate_pseudotime_layout(prefix))])
        )
    if 'paga' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='PAGA', value='paga-tab', children=[html.Div(generate_paga_layout(adata, prefix))])
        )
    if 'volcano' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='Volcano Plot', value='volcano-tab', children=[html.Div(generate_volcano_layout(adata, prefix))])
        )
    if 'grn' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(label='GRN', value='grn-tab', children=[html.Div(generate_grn_demo_layout(adata, prefix))])
        )
    if 'peak-browser' in selected_optional_components:
        tabs_children.append(
            dcc.Tab(
                label='Peak Browser',
                value='peak-browser-tab',
                children=[html.Div(generate_atac_browser_layout(adata, prefix, gene_annotation_path=gene_annotation_path))]
            )
        )

    tabs = dcc.Tabs(
        tabs_children,
        id=f'{prefix}-single-cell-tabs',
        value=initial_tab_value,
        className='custom-tabs'
    )

    no_optional_plots_msg = None
    if not tabs_children:
        no_optional_plots_msg = dbc.Alert(
            "No optional plots enabled for this dataset. Set `optional_plot_components` in the config to enable them.",
            color="secondary",
            style={'marginTop': '10px'}
        )

    return dbc.Row([
        dbc.Col(
            generate_left_control(default_gene_markers, discrete_label_list, prefix),
            xs=12, sm=12, md=4, lg=4, xl=2,
            style={"borderRight": "1px solid #ddd", "padding": "10px"},
        ),
        dbc.Col([tabs, no_optional_plots_msg] if no_optional_plots_msg else [tabs], xs=12, sm=12, md=8, lg=8, xl=10)
    ], style={'marginBottom': '50px'})

# ============= Helper Functions =============

def is_continuous_annotation(adata, annotation, threshold=50):
    """Check if an annotation is continuous based on unique value count and data type."""
    if annotation not in adata.obs.columns:
        return False

    # Any numeric dtype: kind 'i'/'u' (signed/unsigned int) or 'f' (float). This
    # covers uint16/uint32/etc. (e.g. n_genes, n_umis) that an exact dtype-name
    # list misses -- otherwise they're treated as categorical and a column with
    # thousands of unique values renders thousands of traces and crashes.
    dtype = adata.obs[annotation].dtype
    if getattr(dtype, "kind", None) in ("i", "u", "f"):
        n_unique = adata.obs[annotation].nunique()
        return n_unique >= threshold
    return False

# ============= Main Callback Functions =============

def matrix_callbacks(
    app,
    adata,
    prefix,
    embedding_render_backend="scattergl",
    color_config=None,
    gene_annotation_path=None,
):
    """Combined callback registration for both scatter and other plots"""
    embedding_render_backend = str(embedding_render_backend).lower()
    if embedding_render_backend not in {"scattergl", "datashader"}:
        embedding_render_backend = "scattergl"

    # Use the dataset's palette when provided; otherwise fall back to the global default.
    color_config = color_config or _default_color_config

    obs_columns = adata.obs.columns.to_list()
    obs_columns_lower = [c.lower() for c in obs_columns]
    var_names = adata.var_names.to_list()
    var_names_lower = [g.lower() for g in var_names]

    def _resolve_plot_adata_from_filter(filtered_data):
        if (
            filtered_data
            and filtered_data.get('cell_indices') is not None
            and filtered_data.get('n_cells', adata.n_obs) < adata.n_obs
        ):
            return adata[filtered_data['cell_indices']]
        return adata

    def _search_combined(primary, primary_lower, secondary, secondary_lower, query, limit=10):
        q = query.lower()
        primary_hits = [item for item, item_l in zip(primary, primary_lower) if q in item_l]
        secondary_hits = [item for item, item_l in zip(secondary, secondary_lower) if q in item_l]
        return (primary_hits + secondary_hits)[:limit]

    # ===== Cell Selection Hash =====
    # Pre-compute a compact hash of selected_cells once, so downstream plot
    # callbacks can use it from State instead of serialising 50k+ IDs themselves.
    @app.callback(
        Output(f'{prefix}-selected-cells-hash', 'data'),
        Input(f'{prefix}-selected-cells-store', 'data'),
    )
    def update_cells_hash(selected_cells):
        if not selected_cells:
            return None
        n = len(selected_cells)
        payload = json.dumps(selected_cells, separators=(",", ":"), default=str)
        digest = hashlib.md5(payload.encode()).hexdigest()
        return {"len": n, "hash": digest}

    # ===== Global Filter Callbacks =====

    @app.callback(
        [Output(f'{prefix}-global-filter-collapse', 'is_open'),
         Output(f'{prefix}-toggle-global-filter', 'children')],
        Input(f'{prefix}-toggle-global-filter', 'n_clicks'),
        State(f'{prefix}-global-filter-collapse', 'is_open'),
        prevent_initial_call=True
    )
    def toggle_global_filter(n_clicks, is_open):
        if is_open:
            return False, "▼ Show Filters"
        else:
            return True, "▲ Hide Filters"
    
    @app.callback(
        [Output({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value'),
         Output(f'{prefix}-filter-preview', 'children')],
        [Input(f'{prefix}-select-all-filters', 'n_clicks'),
         Input(f'{prefix}-clear-all-filters', 'n_clicks'),
         Input({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value')],
        [State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'options'),
         State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'id')],
        prevent_initial_call=True
    )
    def update_all_filters_and_preview(select_clicks, clear_clicks, current_values, all_options, all_ids):
        ctx = callback_context
        if not ctx.triggered:
            raise PreventUpdate
        
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        values = []
        
        if f'{prefix}-select-all-filters' in trigger_id:
            # Select all values for each filter
            for options in all_options:
                if options:
                    values.append([opt['value'] for opt in options])
                else:
                    values.append([])
        elif f'{prefix}-clear-all-filters' in trigger_id:
            # Clear all filters
            values = [[] for _ in all_options]
        else:
            # Use current values for real-time preview
            values = current_values or []
        
        # Calculate preview cell count in real-time
        if values and all_ids:
            mask = np.ones(adata.n_obs, dtype=bool)
            
            # Apply each filter
            for i, (filter_values, filter_id) in enumerate(zip(values, all_ids)):
                if filter_values:  # Only apply if values are selected
                    column = filter_id['column']
                    col_values = adata.obs[column].astype(str).to_numpy()
                    mask &= np.isin(col_values, filter_values)
            
            # Get preview count
            preview_count = int(mask.sum())
            preview_text = f"Preview: {preview_count:,} cells will be selected"
        else:
            preview_text = ""
        
        return values, preview_text
    
    @app.callback(
        [Output(f'{prefix}-global-filtered-data', 'data'),
         Output(f'{prefix}-global-cell-count', 'children'),
         Output(f'{prefix}-filter-preview', 'children', allow_duplicate=True)],
        Input(f'{prefix}-apply-global-filter', 'n_clicks'),
        [State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'value'),
         State({'type': f'{prefix}-global-metadata-filter', 'column': ALL}, 'id')],
        prevent_initial_call=True
    )
    def apply_global_filter(n_clicks, filter_values, filter_ids):
        if not n_clicks:
            raise PreventUpdate
        
        # Start with all cells
        mask = np.ones(adata.n_obs, dtype=bool)
        
        # Apply each filter
        for i, (values, filter_id) in enumerate(zip(filter_values, filter_ids)):
            if values:  # Only apply if values are selected
                column = filter_id['column']
                col_values = adata.obs[column].astype(str).to_numpy()
                mask &= np.isin(col_values, values)
        
        # Store positional indices (smaller than index labels in payload and faster to apply)
        filtered_indices = np.flatnonzero(mask).astype(np.int32).tolist()
        n_filtered = len(filtered_indices)
        
        # Simple status message without "Filtered by:" details
        preview_text = f"Applied: {n_filtered:,} cells selected"
        
        # Update cell count display
        cell_count_text = f"{n_filtered:,}"
        
        return {
            'cell_indices': filtered_indices,
            'n_cells': n_filtered
        }, cell_count_text, preview_text
    
    # ===== Gene Selection Callbacks =====
    
    @app.callback(
        [Output(f'{prefix}-single-cell-genes-selection', 'style'),
         Output(f'{prefix}-single-cell-genes-textarea', 'style'),
         Output(f'{prefix}-single-cell-genes-selection', 'options', allow_duplicate=True)],
        [Input(f'{prefix}-gene-input-mode', 'value'),
         Input(f'{prefix}-single-cell-genes-selection', 'value')],
        prevent_initial_call=True
    )
    def toggle_gene_input_display(input_mode, selected_genes):
        if input_mode == 'dropdown':
            # When switching to dropdown, make sure the currently selected
            # genes appear as options so the dropdown displays them.
            options = [{'label': g, 'value': g} for g in (selected_genes or [])]
            return {'marginBottom': '15px', 'font-size': '12px'}, {'display': 'none'}, options
        else:
            return {'display': 'none'}, {'width': '100%', 'height': '80px', 'marginBottom': '10px'}, no_update
    
    @app.callback(
        [Output(f'{prefix}-single-cell-genes-selection', 'value', allow_duplicate=True),
         Output(f'{prefix}-gene-input-error', 'children')],
        [Input(f'{prefix}-single-cell-genes-textarea', 'value'),
         Input(f'{prefix}-gene-input-mode', 'value')],
        prevent_initial_call=True
    )
    def validate_gene_input(textarea_value, input_mode):
        if input_mode != 'text' or not textarea_value:
            return no_update, ''
        
        # Parse the textarea input - handle various formats
        # First, check if it looks like a Python list (has brackets)
        text = textarea_value.strip()
        if text.startswith('[') and text.endswith(']'):
            text = text[1:-1]  # Remove brackets
        elif text.startswith('(') and text.endswith(')'):
            text = text[1:-1]  # Remove parentheses
        
        # Split by comma and clean each gene name
        input_genes = []
        for item in text.split(','):
            # Remove various quote types and whitespace
            gene = item.strip()
            # Remove quotes (single, double, backticks, smart quotes)
            gene = gene.strip('"\'`''""')
            if gene:  # Only add non-empty strings
                input_genes.append(gene)
        
        # Get all available genes (case-insensitive mapping)
        available_genes = list(adata.var_names)
        gene_map = {gene.upper(): gene for gene in available_genes}
        
        valid_genes = []
        invalid_genes = []
        seen_genes = set()
        duplicate_genes = []
        
        for gene in input_genes:
            gene_upper = gene.upper()
            if gene_upper in gene_map:
                actual_gene = gene_map[gene_upper]
                if actual_gene not in seen_genes:
                    valid_genes.append(actual_gene)
                    seen_genes.add(actual_gene)
                else:
                    duplicate_genes.append(gene)
            else:
                invalid_genes.append(gene)
        
        error_messages = []
        if invalid_genes:
            error_messages.append(f"Invalid genes not found: {', '.join(invalid_genes)}")
        if duplicate_genes:
            error_messages.append(f"Duplicate genes removed: {', '.join(duplicate_genes)}")
        
        error_message = ' | '.join(error_messages)
        
        return valid_genes, error_message
    
    @app.callback(
        Output(f'{prefix}-single-cell-genes-textarea', 'value'),
        Input(f'{prefix}-single-cell-genes-selection', 'value'),
        State(f'{prefix}-gene-input-mode', 'value'),
        prevent_initial_call=True
    )
    def sync_dropdown_to_textarea(dropdown_value, input_mode):
        if input_mode == 'dropdown' and dropdown_value:
            return ', '.join(dropdown_value)
        return no_update
    
    # ===== Scatter Plot Callbacks =====
    register_scatter_callbacks(
        app,
        adata,
        prefix,
        embedding_render_backend=embedding_render_backend,
        initialize_scatter_components=build_initialize_scatter_components,
        apply_relayout=apply_relayout,
        is_continuous_annotation=is_continuous_annotation,
        resolve_plot_adata_from_filter=_resolve_plot_adata_from_filter,
        search_combined=_search_combined,
        obs_columns=obs_columns,
        obs_columns_lower=obs_columns_lower,
        var_names=var_names,
        var_names_lower=var_names_lower,
        palette_json=palette_json,
        color_config=color_config,
        plot_embedding=plot_embedding,
        plot_coexpression_embedding=plot_coexpression_embedding,
    )

    # ===== Other Plots Callbacks =====

    @app.callback(
        Output(f'{prefix}-single-cell-genes-selection', 'options', allow_duplicate=True),
        Input(f'{prefix}-single-cell-genes-selection', 'search_value'),
        State(f'{prefix}-single-cell-genes-selection', 'value'),
        prevent_initial_call=True
    )
    def update_genes_dropdown(search_value, value):
        if not search_value:
            raise PreventUpdate
        label_list = adata.var_names.to_list()
        matching_labels = [label for label in label_list if search_value.lower() in label.lower()]
        selected_labels = value if value else []
        all_labels = list(set(selected_labels + matching_labels[:10]))
        return [{'label': label, 'value': label} for label in all_labels]
    
    @app.callback(
        [Output(f'{prefix}-single-cell-label-selection', 'options'),
         Output(f'{prefix}-single-cell-label-selection', 'value')],
        [Input(f'{prefix}-single-cell-annotation-dropdown', 'value'),
         Input(f'{prefix}-selected-cells-store', 'data')]
    )
    def update_labels_based_on_annotation(selected_annotation, selected_cells):
        # Filter adata if cells are selected
        src = adata[selected_cells] if selected_cells else adata
        unique_labels = sorted_categories(src, selected_annotation)
        label_options = [{'label': label, 'value': label} for label in unique_labels]
        return label_options, list(unique_labels)

    register_heatmap_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_unified_heatmap=plot_unified_heatmap,
        palette_json=palette_json,
        color_config=color_config,
        make_cache_key=_make_cache_key,
        hash_list_signature=_hash_list_signature,
        cached_figure_get=_cached_figure_get,
        cached_figure_set=_cached_figure_set,
    )

    register_dotplot_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_dot_matrix=plot_dot_matrix,
        make_cache_key=_make_cache_key,
        hash_list_signature=_hash_list_signature,
        cached_figure_get=_cached_figure_get,
        cached_figure_set=_cached_figure_set,
    )

    register_pseudotime_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_genes_in_pseudotime=plot_genes_in_pseudotime,
        palette_json=palette_json,
        color_config=color_config,
        make_cache_key=_make_cache_key,
        hash_list_signature=_hash_list_signature,
        cached_figure_get=_cached_figure_get,
        cached_figure_set=_cached_figure_set,
    )

    register_violin_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_violin1=plot_violin1,
        plot_violin2_new=plot_violin2_new,
        palette_json=palette_json,
        var_names=var_names,
        var_names_lower=var_names_lower,
        color_config=color_config,
    )

    register_stacked_bar_callbacks(
        app,
        adata,
        prefix,
        filter_data=filter_data,
        plot_stacked_bar=plot_stacked_bar,
        palette_json=palette_json,
        color_config=color_config,
    )

    register_paga_callbacks(
        app,
        adata,
        prefix,
        build_paga_cytoscape=build_paga_cytoscape,
        color_config=color_config,
    )

    register_volcano_callbacks(
        app,
        adata,
        prefix,
    )

    register_grn_demo_callbacks(
        app,
        adata,
        prefix,
    )

    if has_genomic_peak_features(adata):
        register_atac_browser_callbacks(
            app,
            adata,
            prefix,
            gene_annotation_path=gene_annotation_path,
            color_config=color_config,
        )

    # Add callback to update filter status
    @app.callback(
        Output(f'{prefix}-filter-status', 'children'),
        Input(f'{prefix}-selected-cells-store', 'data')
    )
    def update_filter_status(selected_cells):
        if selected_cells:
            n_selected = len(selected_cells)
            n_total = adata.n_obs
            return dbc.Alert(
                f"Showing {n_selected} of {n_total} cells ({n_selected/n_total*100:.1f}%) based on scatter plot selection",
                color="info",
                dismissable=False,
                style={'margin': '0'}
            )
        return None
