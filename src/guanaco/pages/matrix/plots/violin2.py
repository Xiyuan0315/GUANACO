import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from statsmodels.formula.api import ols, mixedlm
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
import warnings

from guanaco.pages.matrix.plots.violin1 import (
    VIOLIN_FILL_ALPHA,
    OVERLAY_FILL_COLOR,
    OVERLAY_LINE_COLOR,
    _to_rgba,
)
from guanaco.utils.gene_extraction_utils import extract_gene_expression, apply_transformation
from guanaco.data.loader import obs_col

DEFAULT_COLORS = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"]


def _violin_style(fillcolor, show_box, show_points=False):
    """Per-trace styling shared with the stacked violin (violin1).

    Gives both violin tabs the same look: a translucent category fill, a
    DarkSlateGrey outline, a mean line, and -- when requested -- a neutral, narrow
    box overlay (median + IQR). Individual points are never drawn; a violin is a KDE
    of the distribution.
    """
    return dict(
        fillcolor=_to_rgba(fillcolor, VIOLIN_FILL_ALPHA),
        line_color='DarkSlateGrey',
        meanline_visible=True,
        points='all' if show_points else False,
        box=dict(
            visible=show_box,
            width=0.12,
            fillcolor=OVERLAY_FILL_COLOR,
            line=dict(color=OVERLAY_LINE_COLOR, width=1),
        ),
    )


# Human-readable names for a manually chosen test (auto-selection already returns a
# description from determine_test_method). Used to label the plot with the test run.
_TEST_LABELS = {
    'mwu-test': 'Mann-Whitney U test',
    'ttest': "Welch's t-test",
    'kw-test': 'Kruskal-Wallis test',
    'anova': 'One-way ANOVA',
    'linear-model': 'Linear model (expression ~ obs1 + obs2)',
    'linear-model-interaction': 'Linear model with interaction',
    'mixed-model': 'Mixed model (expression ~ obs1 + (1|obs2))',
}


def _p_to_stars(p):
    """Significance stars for a p/q-value (publication convention)."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ''
    if p <= 0.001:
        return '***'
    if p <= 0.01:
        return '**'
    if p <= 0.05:
        return '*'
    return 'ns'


def _sig_inline(p, label='p'):
    """One-line significance text, e.g. ``p = 0.003 **`` (single top annotation)."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ''
    return f"{label} = {p:.3g} {_p_to_stars(p)}".strip()


def _sig_stacked(p, label='p'):
    """Stars over the numeric value, for per-facet annotations above each group."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ''
    return f"<b>{_p_to_stars(p)}</b><br><span style='font-size:9px;color:gray'>{label}={p:.2g}</span>"


def _fdr_adjust(p_values):
    """Benjamini-Hochberg FDR correction across mode2's per-facet comparisons.

    Returns ``(p_values, adjusted)``. Only numeric facet entries are corrected;
    bookkeeping keys (``error`` / ``model_summary`` / ...) are left untouched. With
    fewer than two valid comparisons there is nothing to correct.
    """
    skip = {'error', 'model_summary', 'method', 'overall'}
    keys = [
        k for k, v in p_values.items()
        if k not in skip and isinstance(v, (int, float)) and not (isinstance(v, float) and np.isnan(v))
    ]
    if len(keys) < 2:
        return p_values, False
    adj = multipletests([float(p_values[k]) for k in keys], method='fdr_bh')[1]
    out = dict(p_values)
    for k, q in zip(keys, adj):
        out[k] = float(q)
    return out, True


def determine_test_method(meta1_levels, meta2_levels, mode, test_override=None):

    if test_override and test_override != 'auto':
        return test_override, ""
    
    if mode == 'mode1':  # One metadata only
        if meta1_levels == 2:
            return 'mwu-test', "Mann-Whitney U test (2 groups)"
        else:
            return 'kw-test', "Kruskal-Wallis test (>2 groups)"
    
    elif mode == 'mode2':  # Facet by meta1, compare meta2
        if meta2_levels == 2:
            return 'mwu-test', "Mann-Whitney U test within each facet"
        else:
            return 'kw-test', "Kruskal-Wallis test within each facet"
    
    elif mode == 'mode3':  # Linear model
        return 'linear-model', "Linear model: expression ~ meta1 + meta2"
    
    elif mode == 'mode4':  # Mixed model
        return 'mixed-model', "Mixed model: expression ~ meta1 + (1|meta2)"
    
    return 'none', ""


def _expression_groups(df, group_col, categories=None):
    groups = []
    group_names = []
    for group in (categories if categories else sorted(df[group_col].unique())):
        group_data = df[df[group_col] == group]['Expression'].values
        if len(group_data) > 0:
            groups.append(group_data)
            group_names.append(group)
    return groups, group_names


def _safe_stat_test(test_func, *groups, **kwargs):
    try:
        _, p_val = test_func(*groups, **kwargs)
        return p_val
    except Exception:
        return np.nan


def _p_value_for_groups(groups, test_method):
    if len(groups) < 2:
        return None
    if test_method == 'mwu-test' and len(groups) == 2:
        return _safe_stat_test(mannwhitneyu, groups[0], groups[1], alternative='two-sided')
    if test_method == 'ttest' and len(groups) == 2:
        return _safe_stat_test(ttest_ind, groups[0], groups[1], equal_var=False)
    if test_method == 'kw-test':
        return _safe_stat_test(kruskal, *groups)
    if test_method == 'anova':
        return _safe_stat_test(f_oneway, *groups)
    return None


def _calculate_mode2_p_values(df, meta1, meta2, test_method, labels=None):
    p_values = {}
    for m1_group in (labels if labels else sorted(df[meta1].unique())):
        facet_data = df[df[meta1] == m1_group]
        groups, _ = _expression_groups(facet_data, meta2)
        p_val = _p_value_for_groups(groups, test_method)
        if p_val is not None:
            p_values[m1_group] = p_val
    return p_values


def _linear_model_p_values(df, meta1, meta2, interaction=False):
    p_values = {}
    try:
        df_model = df.copy()
        df_model[meta1] = pd.Categorical(df_model[meta1])
        df_model[meta2] = pd.Categorical(df_model[meta2])

        if interaction:
            formula = f'Expression ~ C({meta1}) * C({meta2})'
            model = ols(formula, data=df_model).fit()
            anova_table = anova_lm(model, typ=2)
            p_values['model_summary'] = {
                'meta1_p': anova_table.loc[f'C({meta1})', 'PR(>F)'],
                'meta2_p': anova_table.loc[f'C({meta2})', 'PR(>F)'],
                'interaction_p': anova_table.loc[f'C({meta1}):C({meta2})', 'PR(>F)'],
                'r_squared': model.rsquared,
                'aic': model.aic
            }
            p_values['overall'] = p_values['model_summary']['interaction_p']
            return p_values

        formula = f'Expression ~ C({meta1}) + C({meta2})'
        model = ols(formula, data=df_model).fit()
        p_values['model_summary'] = {
            'meta1_p': model.pvalues[f'C({meta1})[T.{df_model[meta1].cat.categories[1]}]'],
            'meta2_p': model.pvalues[f'C({meta2})[T.{df_model[meta2].cat.categories[1]}]'] if len(df_model[meta2].unique()) > 1 else np.nan,
            'r_squared': model.rsquared,
            'aic': model.aic
        }
        p_values['overall'] = p_values['model_summary']['meta1_p']
    except Exception as e:
        p_values['error'] = str(e)
        p_values['overall'] = np.nan
    return p_values


def _pseudobulk_p_values(df, meta1, meta2):
    pseudobulk = df.groupby([meta1, meta2])['Expression'].mean().reset_index()
    groups, _ = _expression_groups(pseudobulk, meta1)
    if len(groups) == 2:
        return {
            'overall': _safe_stat_test(ttest_ind, groups[0], groups[1], equal_var=False),
            'method': 'pseudobulk_ttest',
        }
    if len(groups) > 2:
        return {
            'overall': _safe_stat_test(f_oneway, *groups),
            'method': 'pseudobulk_anova',
        }
    return {'overall': np.nan, 'method': 'pseudobulk'}


def _mixed_model_p_values(df, meta1, meta2):
    p_values = {}
    try:
        df_model = df.copy()
        df_model[meta1] = pd.Categorical(df_model[meta1])
        df_model[meta2] = pd.Categorical(df_model[meta2])

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                formula = f'Expression ~ C({meta1})'
                model = mixedlm(formula, data=df_model, groups=df_model[meta2]).fit(reml=False)

                p_values['model_summary'] = {
                    'meta1_p': model.pvalues[f'C({meta1})[T.{df_model[meta1].cat.categories[1]}]'],
                    'log_likelihood': model.llf,
                    'converged': model.converged
                }
                p_values['overall'] = p_values['model_summary']['meta1_p']
        except Exception:
            p_values.update(_pseudobulk_p_values(df, meta1, meta2))
    except Exception as e:
        p_values['error'] = str(e)
        p_values['overall'] = np.nan
    return p_values


def calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels=None, design_type='crossed'):
    if mode == 'mode1':  # One metadata only
        groups, _ = _expression_groups(df, meta1, labels)
        p_val = _p_value_for_groups(groups, test_method)
        return {'overall': p_val} if p_val is not None else {}

    if mode == 'mode2':  # Facet by meta1, compare meta2
        return _calculate_mode2_p_values(df, meta1, meta2, test_method, labels)

    if mode == 'mode3':  # Linear model
        return _linear_model_p_values(
            df,
            meta1,
            meta2,
            interaction=test_method == 'linear-model-interaction',
        )

    if mode == 'mode4':  # Mixed model
        return _mixed_model_p_values(df, meta1, meta2)

    return {}


def assign_colors(levels, color_map=None, palette=None):
    color_map = dict(color_map or {})
    
    if palette:
        default_colors = palette
    else:
        default_colors = DEFAULT_COLORS
        
    for i, level in enumerate(sorted(levels)):
        if level not in color_map:
            color_map[level] = default_colors[i % len(default_colors)]
    return color_map


def _empty_message_figure(message):
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=14, color='gray'),
    )
    return fig


def _build_expression_frame(adata, key, meta1, meta2, mode, transformation=None, layer=None, labels=None):
    expression_data = extract_gene_expression(adata, key, layer=layer)
    if transformation:
        expression_data = apply_transformation(expression_data, transformation, copy=False)

    data = {
        'Expression': expression_data,
        meta1: obs_col(adata.obs, meta1).values,
    }
    if mode != 'mode1' and meta2:
        data[meta2] = obs_col(adata.obs, meta2).values

    df = pd.DataFrame(data)
    if labels:
        df = df[df[meta1].isin(labels)]
    return df


def _resolve_color_map(design_type, meta1_cats, meta2_cats, color_map=None, palette=None):
    if design_type == 'nested' or not meta2_cats:
        return assign_colors(meta1_cats, color_map, palette=palette)
    return assign_colors(meta2_cats, color_map, palette=palette)


def _resolve_test_method(test_method, meta1_levels, meta2_levels, mode):
    if test_method == 'auto':
        return determine_test_method(meta1_levels, meta2_levels, mode)
    return test_method, _TEST_LABELS.get(test_method, "")


def _create_violin_for_mode(
    fig,
    mode,
    design_type,
    df,
    key,
    meta1,
    meta2,
    meta1_cats,
    meta2_cats,
    color_map,
    show_box,
    show_points,
):
    if mode == 'mode1':
        return _create_single_meta_violin(
            fig,
            df,
            key,
            meta1,
            meta1_cats,
            show_box,
            color_map,
            show_points,
        )
    if design_type == 'nested':
        return _create_nested_violin(
            fig,
            df,
            key,
            meta1,
            meta2,
            meta1_cats,
            meta2_cats,
            color_map,
            show_box,
            show_points,
        )
    if len(meta2_cats) == 2:
        return _create_split_violin(
            fig,
            df,
            key,
            meta1,
            meta2,
            meta1_cats,
            meta2_cats,
            color_map,
            show_box,
            show_points,
        )
    return _create_crossed_violin(
        fig,
        df,
        key,
        meta1,
        meta2,
        meta1_cats,
        meta2_cats,
        color_map,
        show_box,
        show_points,
    )


def _uses_split_violin(mode, meta2, meta2_cats, design_type):
    return (
        mode in ['mode2', 'mode3', 'mode4']
        and meta2
        and len(meta2_cats) == 2
        and design_type in ['crossed', 'sparse']
    )


def _expression_axis_range(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return [0, 1]

    y_min = float(values.min())
    y_max = float(values.max())
    if y_min == y_max:
        pad = 1e-3 if y_max == 0 else abs(y_max) * 0.05
        if y_min >= 0:
            return [0, y_max + pad]
        return [y_min - pad, y_max + pad]

    if y_min >= 0:
        return [0, y_max * 1.15]

    pad = (y_max - y_min) * 0.05
    return [y_min - pad, y_max + pad]


def _annotation_y(values):
    return _expression_axis_range(values)[1]


P_VALUE_METADATA_KEYS = {'error', 'model_summary', 'method', 'overall'}


def _is_valid_facet_p_value(facet, p_val, meta1_groups):
    if facet in P_VALUE_METADATA_KEYS or facet not in meta1_groups:
        return False
    try:
        return not np.isnan(p_val)
    except TypeError:
        return False


def _add_paper_annotation(fig, text, size=12):
    fig.add_annotation(
        x=0.5,
        y=1.02,
        xref='paper',
        yref='paper',
        text=text,
        showarrow=False,
        font=dict(size=size, color='DarkSlateGrey')
    )


def _add_data_annotation(fig, x, y, text, xref=None):
    annotation = dict(
        x=x,
        y=y,
        text=text,
        showarrow=False,
        font=dict(size=10, color='DarkSlateGrey'),
    )
    if xref is not None:
        annotation['xref'] = xref
    fig.add_annotation(**annotation)


def _sorted_nested_meta2(df, meta1, meta2, meta1_groups):
    meta2_groups = sorted(df[meta2].unique())
    meta2_to_meta1 = df.groupby(meta2)[meta1].first().to_dict()
    return sorted(meta2_groups, key=lambda x: (
        meta1_groups.index(meta2_to_meta1.get(x, meta1_groups[0]))
        if meta2_to_meta1.get(x) in meta1_groups
        else 999,
        str(x)
    )), meta2_to_meta1


def _add_nested_mode2_annotations(fig, p_values, df, meta1, meta2, meta1_groups, y_max, p_label):
    sorted_meta2, meta2_to_meta1 = _sorted_nested_meta2(df, meta1, meta2, meta1_groups)

    for m1_group in meta1_groups:
        if m1_group not in p_values:
            continue
        p_val = p_values[m1_group]
        if np.isnan(p_val):
            continue

        m2_in_group = [m2 for m2 in sorted_meta2 if meta2_to_meta1.get(m2) == m1_group]
        if not m2_in_group:
            continue

        center_m2 = m2_in_group[len(m2_in_group) // 2]
        _add_data_annotation(fig, center_m2, y_max, _sig_stacked(p_val, p_label))


def _add_split_mode2_annotations(fig, p_values, meta1_groups, y_max, p_label):
    for facet, p_val in p_values.items():
        if _is_valid_facet_p_value(facet, p_val, meta1_groups):
            _add_data_annotation(fig, facet, y_max, _sig_stacked(p_val, p_label))


def _add_crossed_mode2_annotations(fig, p_values, df, meta1_groups, meta2, y_max, p_label):
    meta2_groups = sorted(df[meta2].unique())
    center_m2 = meta2_groups[len(meta2_groups) // 2]
    for facet, p_val in p_values.items():
        if _is_valid_facet_p_value(facet, p_val, meta1_groups):
            _add_data_annotation(fig, f"{facet}_{center_m2}", y_max, _sig_stacked(p_val, p_label))


def _add_grouped_mode2_annotations(fig, p_values, df, meta1_groups, meta2, y_max, p_label):
    meta2_groups = sorted(df[meta2].unique())
    n_meta2 = len(meta2_groups)
    for facet, p_val in p_values.items():
        if _is_valid_facet_p_value(facet, p_val, meta1_groups):
            facet_idx = meta1_groups.index(facet)
            center_pos = facet_idx * n_meta2 + (n_meta2 - 1) / 2
            _add_data_annotation(fig, center_pos, y_max, _sig_stacked(p_val, p_label), xref='x')


def _add_mode2_annotations(fig, p_values, df, meta1, meta2, split_violin, design_type, y_max, p_label):
    meta1_groups = sorted(df[meta1].unique())

    if design_type == 'nested':
        _add_nested_mode2_annotations(fig, p_values, df, meta1, meta2, meta1_groups, y_max, p_label)
    elif split_violin:
        _add_split_mode2_annotations(fig, p_values, meta1_groups, y_max, p_label)
    elif design_type in ['crossed', 'sparse']:
        _add_crossed_mode2_annotations(fig, p_values, df, meta1_groups, meta2, y_max, p_label)
    else:
        _add_grouped_mode2_annotations(fig, p_values, df, meta1_groups, meta2, y_max, p_label)


def _model_summary_text(summary, mode, meta1, meta2):
    p_text_parts = []
    if 'meta1_p' in summary:
        p1 = summary['meta1_p']
        p_text_parts.append(f"{meta1}: {p1:.3g} {_p_to_stars(p1)}".strip())
    if mode == 'mode3' and 'meta2_p' in summary and not np.isnan(summary['meta2_p']):
        p2 = summary['meta2_p']
        p_text_parts.append(f"{meta2}: {p2:.3g} {_p_to_stars(p2)}".strip())
    if 'interaction_p' in summary and not np.isnan(summary['interaction_p']):
        p_int = summary['interaction_p']
        p_text_parts.append(f"interaction: {p_int:.3g} {_p_to_stars(p_int)}".strip())
    return ", ".join(p_text_parts)


def _add_model_annotations(fig, p_values, mode, meta1, meta2):
    if 'model_summary' in p_values:
        _add_paper_annotation(fig, _model_summary_text(p_values['model_summary'], mode, meta1, meta2))
    elif 'overall' in p_values:
        _add_paper_annotation(fig, _sig_inline(p_values['overall'], 'p'))


def add_p_value_annotations_new(fig, p_values, df, mode, meta1=None, meta2=None, split_violin=False, design_type='crossed', adjusted=False):
    if not p_values:
        return

    # Per-facet values shown above each group are q-values when FDR-adjusted.
    p_label = 'q' if adjusted else 'p'

    y_max = _annotation_y(df['Expression'])

    if mode == 'mode1':
        # Single p-value annotation
        if 'overall' in p_values and not np.isnan(p_values['overall']):
            _add_paper_annotation(fig, _sig_inline(p_values['overall'], 'p'))
    
    elif mode == 'mode2':
        _add_mode2_annotations(fig, p_values, df, meta1, meta2, split_violin, design_type, y_max, p_label)
    
    elif mode in ['mode3', 'mode4']:
        _add_model_annotations(fig, p_values, mode, meta1, meta2)


def plot_violin2_new(adata, key, meta1, meta2, mode, transformation=None, layer=None,
                     show_box=False, test_method='auto',
                     labels=None, color_map=None, palette=None, show_points=False):
    """
    Create violin plots for gene expression data with support for multiple metadata designs.
    
    Handles three design patterns:
    1. Crossed (Orthogonal): Every level of meta1 co-occurs with every level of meta2
       Example: Cluster x Condition (all clusters have both Control and Treated)
    2. Nested (Hierarchical): Levels of meta2 belong exclusively to specific levels of meta1
       Example: Condition x Patient (P1,P2 in Control; P3,P4 in Disease)
    3. Sparse (Mixed): Some overlaps exist, but not all combinations are present
       Example: Tissue x CellType (Neurons in Brain but not Liver)
    
    """
    
    # Early validation
    if not key or not meta1:
        return go.Figure()
    
    if key not in adata.var_names:
        return go.Figure()

    if mode not in ['mode1', 'mode2', 'mode3', 'mode4']:
        return go.Figure()

    if mode in ['mode2', 'mode3', 'mode4'] and not meta2:
        return go.Figure()

    df = _build_expression_frame(adata, key, meta1, meta2, mode, transformation, layer, labels)

    # Handle empty dataframe
    if len(df) == 0:
        return _empty_message_figure("No data available for selected filters")
    
    # Get unique categories (preserve order if categorical)
    meta1_cats = _get_ordered_categories(df, meta1)
    meta2_cats = _get_ordered_categories(df, meta2) if mode != 'mode1' and meta2 else []
    
    # Detect design pattern to determine plotting strategy
    design_type = _detect_design_pattern(df, meta1, meta2) if meta2 and mode != 'mode1' else 'single'
    
    color_map = _resolve_color_map(design_type, meta1_cats, meta2_cats, color_map, palette)

    # Determine test method
    test_method, test_description = _resolve_test_method(
        test_method,
        len(meta1_cats),
        len(meta2_cats),
        mode,
    )

    # Subsample the rows that feed the violin shapes (KDE) so the payload and
    # client-side cost stay bounded on large datasets -- the same trick violin1
    # uses. The full ``df`` is still used for the statistical tests, p-value
    # annotations and layout below; only the drawn violins use ``df_plot``.
    group_cols = [meta1] if (mode == 'mode1' or not meta2) else [meta1, meta2]
    df_plot = _downsample_for_violins(df, group_cols)

    # Create figure based on mode
    fig = go.Figure()
    fig, title, xaxis_title = _create_violin_for_mode(
        fig,
        mode,
        design_type,
        df_plot,
        key,
        meta1,
        meta2,
        meta1_cats,
        meta2_cats,
        color_map,
        show_box,
        show_points,
    )

    # Calculate and add p-values
    if test_method and test_method != 'none':
        p_values = calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels, design_type)
        adjusted = False
        if mode == 'mode2':
            # mode2 runs one comparison per meta1 facet -> correct for multiple testing.
            p_values, adjusted = _fdr_adjust(p_values)
            if adjusted and test_description:
                test_description += " · FDR-adjusted (BH)"
        split_violin = _uses_split_violin(mode, meta2, meta2_cats, design_type)
        add_p_value_annotations_new(fig, p_values, df, mode, meta1, meta2, split_violin, design_type, adjusted=adjusted)

    _apply_common_layout(fig, df, key, title, xaxis_title, meta2, mode, test_description)

    return fig


def _get_ordered_categories(df, col):
    """Get ordered categories from a column, preserving categorical order if available."""
    if col not in df.columns:
        return []
    
    if hasattr(df[col], 'cat'):
        return list(df[col].cat.categories)
    else:
        return sorted(df[col].unique())


def _detect_design_pattern(df, meta1, meta2):
    """
    Detect the relationship between two metadata columns.
    
    Returns:
    --------
    str: 'crossed', 'nested', or 'sparse'
    """
    if meta2 not in df.columns:
        return 'single'
    
    # Create contingency table
    contingency = pd.crosstab(df[meta1], df[meta2])
    
    # Count non-zero cells
    n_meta1 = len(contingency.index)
    n_meta2 = len(contingency.columns)
    n_nonzero = (contingency > 0).sum().sum()
    n_total = n_meta1 * n_meta2
    
    # Determine pattern
    if n_nonzero == n_total:
        return 'crossed'  # All combinations present
    elif n_nonzero == max(n_meta1, n_meta2):
        return 'nested'   # One-to-one or hierarchical
    else:
        return 'sparse'   # Some combinations missing


def _get_violin_spanmode(data):
    """Determine appropriate spanmode based on data variance."""
    variance = np.var(data)
    if variance < 1e-10:
        return 'soft'
    return 'hard'


# A violin is a KDE, so a uniform random subsample per drawn violin reproduces the
# same shape while bounding the figure payload and client-side KDE cost -- the same
# trick violin1 uses. The statistical tests still run on the full data; only the
# rows fed to the go.Violin traces are subsampled. scalemode='width' means fewer
# points don't change the drawn width.
MAX_POINTS_PER_VIOLIN = 2_000


def _downsample_for_violins(df, group_cols, cap=MAX_POINTS_PER_VIOLIN, seed=315):
    """Keep at most ``cap`` rows per (drawn-violin) group; preserves every group.

    ``group_cols`` are the columns that define one violin (``[meta1]`` for the
    single-meta mode, ``[meta1, meta2]`` otherwise). Groups with <= ``cap`` rows are
    kept whole, so every non-empty combination still produces a violin (and nested
    designs keep their meta2->meta1 mapping). A fixed seed keeps shapes stable
    across re-renders.
    """
    if len(df) <= cap:
        return df
    rng = np.random.default_rng(seed)
    keep = []
    for _, grp in df.groupby(group_cols, observed=True, sort=False):
        idx = grp.index.to_numpy()
        if idx.size > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        keep.append(idx)
    if not keep:
        return df
    return df.loc[np.concatenate(keep)]


def _create_single_meta_violin(fig, df, key, meta1, meta1_cats, show_box, color_map=None, show_points=False):
    """Create violin plot with single metadata (mode1)."""
    color_map = color_map or {}

    for i, group in enumerate(meta1_cats):
        group_data = df[df[meta1] == group]['Expression']

        if len(group_data) == 0:
            continue

        spanmode = _get_violin_spanmode(group_data)

        fillcolor = color_map.get(group, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        fig.add_trace(go.Violin(
            x=[group] * len(group_data),
            y=group_data,
            name=str(group),
            width=0.8,
            spanmode=spanmode,
            scalemode='width',
            hoveron='violins',
            hoverinfo='y',
            jitter=0.05,
            showlegend=False,
            **_violin_style(fillcolor, show_box, show_points),
        ))
    
    title = f"Expression of {key} by {meta1}"
    xaxis_title = meta1
    
    return fig, title, xaxis_title


def _create_split_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                         color_map, show_box, show_points=False):
    """
    Create split violin plot for binary meta2 comparison.
    
    Handles: Crossed, Nested, and Sparse designs by checking data existence.
    """
    
    # Track which meta2 groups have been added to legend
    legend_added = set()
    
    for m1_group in meta1_cats:
        facet_data = df[df[meta1] == m1_group]
        
        if len(facet_data) == 0:
            continue
        
        for j, m2_group in enumerate(meta2_cats):
            subset = facet_data[facet_data[meta2] == m2_group]
            
            # CRITICAL: Skip empty combinations (handles Nested/Sparse designs)
            if len(subset) == 0:
                continue
            
            side = 'negative' if j == 0 else 'positive'
            spanmode = _get_violin_spanmode(subset['Expression'])
            
            # Only show legend for first occurrence of each meta2 group
            show_legend = m2_group not in legend_added
            if show_legend:
                legend_added.add(m2_group)
            
            fig.add_trace(go.Violin(
                x=[m1_group] * len(subset),
                y=subset['Expression'],
                name=str(m2_group),
                legendgroup=str(m2_group),  # Group legend items
                scalegroup=str(m1_group),   # Scale within each meta1 group
                scalemode='width',
                side=side,
                width=0.8,
                spanmode=spanmode,
                hoveron='violins',
                hoverinfo='y',
                jitter=0.05,
                showlegend=show_legend,
                **_violin_style(color_map.get(m2_group, DEFAULT_COLORS[j % len(DEFAULT_COLORS)]), show_box, show_points),
            ))
    
    fig.update_layout(violinmode='overlay')
    
    title = f"Expression of {key} by {meta1} and {meta2}"
    xaxis_title = meta1
    
    return fig, title, xaxis_title


def _create_nested_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                          color_map, show_box, show_points=False):
    """
    Create violin plot for nested designs (e.g., Condition x Patient).
    
    In nested designs, meta2 values belong exclusively to specific meta1 values.
    Example: Patient1, Patient2 belong to 'Healthy'; Patient3, Patient4 belong to 'Severe'
    
    Strategy: Use meta2 (Patient) as x-axis, color by meta1 (Condition).
    This ensures each patient gets its own violin, grouped visually by condition.
    """
    
    # Track which meta1 groups have been added to legend
    legend_added = set()
    
    # Build mapping of meta2 -> meta1 for proper ordering
    meta2_to_meta1 = df.groupby(meta2)[meta1].first().to_dict()
    
    # Sort meta2 by their parent meta1 category for visual grouping
    sorted_meta2 = sorted(meta2_cats, key=lambda x: (
        meta1_cats.index(meta2_to_meta1.get(x, meta1_cats[0])) if meta2_to_meta1.get(x) in meta1_cats else 999,
        str(x)
    ))
    
    for m2_group in sorted_meta2:
        subset = df[df[meta2] == m2_group]
        
        if len(subset) == 0:
            continue
        
        # Get the parent meta1 category for this meta2
        m1_group = meta2_to_meta1.get(m2_group)
        if m1_group is None:
            continue
        
        spanmode = _get_violin_spanmode(subset['Expression'])
        
        # Get color index for meta1
        m1_idx = meta1_cats.index(m1_group) if m1_group in meta1_cats else 0
        
        # Only show legend for first occurrence of each meta1 group
        show_legend = m1_group not in legend_added
        if show_legend:
            legend_added.add(m1_group)
        
        fig.add_trace(go.Violin(
            x=[str(m2_group)] * len(subset),
            y=subset['Expression'],
            name=str(m1_group),  # Legend shows meta1 (Condition)
            legendgroup=str(m1_group),  # Group by meta1 for toggle
            width=0.8,
            scalemode='width',
            spanmode=spanmode,
            hoveron='violins',
            hoverinfo='y+name',
            jitter=0.05,
            showlegend=show_legend,
            **_violin_style(color_map.get(m1_group, DEFAULT_COLORS[m1_idx % len(DEFAULT_COLORS)]), show_box, show_points),
        ))
    
    # No grouping needed - each meta2 has its own x position
    fig.update_layout(
        violingap=0.1
    )
    
    title = f"Expression of {key} by {meta2} (colored by {meta1})"
    xaxis_title = meta2
    
    return fig, title, xaxis_title


def _create_crossed_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                           color_map, show_box, show_points=False):
    """
    Create violin plot for crossed designs (e.g., Days × Patients).
    
    In crossed designs, each meta2 value exists for all meta1 values.
    Example: Patient1, Patient2, Patient3 each have data for Day1, Day3, Day7
    
    Strategy: Create separate x-axis position for each combination (9 violins),
    grouped visually by meta1, colored by meta2.
    X-axis labels: day1_patient1, day1_patient2, day1_patient3, day3_patient1, ...
    """
    
    # Track which meta2 groups have been added to legend
    legend_added = set()
    
    # Create x-axis labels for each combination
    x_labels = []
    for m1_group in meta1_cats:
        for m2_group in meta2_cats:
            x_labels.append(f"{m1_group}_{m2_group}")
    
    for m1_group in meta1_cats:
        facet_data = df[df[meta1] == m1_group]
        
        if len(facet_data) == 0:
            continue
        
        for j, m2_group in enumerate(meta2_cats):
            subset = facet_data[facet_data[meta2] == m2_group]
            
            # Skip empty combinations
            if len(subset) == 0:
                continue
            
            spanmode = _get_violin_spanmode(subset['Expression'])
            
            # Create combined x-label
            x_label = f"{m1_group}_{m2_group}"
            
            # Only show legend for first occurrence of each meta2 group
            show_legend = m2_group not in legend_added
            if show_legend:
                legend_added.add(m2_group)
            
            fig.add_trace(go.Violin(
                x=[x_label] * len(subset),
                y=subset['Expression'],
                name=str(m2_group),
                legendgroup=str(m2_group),  # Group legend items for toggle
                width=0.8,
                scalemode='width',
                spanmode=spanmode,
                hoveron='violins',
                hoverinfo='y+name',
                jitter=0.05,
                showlegend=show_legend,
                **_violin_style(color_map.get(m2_group, DEFAULT_COLORS[j % len(DEFAULT_COLORS)]), show_box, show_points),
            ))
    
    # Set x-axis category order to maintain grouping
    fig.update_xaxes(categoryorder='array', categoryarray=x_labels)
    
    fig.update_layout(
        violingap=0.1
    )
    
    title = f"Expression of {key} by {meta1} and {meta2}"
    xaxis_title = f"{meta1} × {meta2}"
    
    return fig, title, xaxis_title


def _apply_common_layout(fig, df, key, title, xaxis_title, meta2, mode, test_description=""):
    """Apply common layout settings to the figure."""
    y_range = _expression_axis_range(df['Expression'])

    # Surface which statistical test was run (esp. for "Auto") as a title subtitle.
    title_text = title
    if test_description:
        title_text = f"{title}<br><span style='font-size:11px;color:gray'>Test: {test_description}</span>"

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12, color='DarkSlateGrey'),
        title=dict(
            text=title_text,
            x=0.5,
            xanchor='center',
            y=0.95,
            font=dict(size=16, color='DarkSlateGrey')
        ),
        yaxis_range=y_range,
        xaxis_title=xaxis_title,
        yaxis_title='Expression',
        legend=dict(
            font=dict(size=12),
            title=dict(text=meta2 if mode != 'mode1' and meta2 else "", font=dict(size=12)),
            tracegroupgap=5
        ),
        margin=dict(b=80, t=100, l=60, r=40)
    )

    fig.update_yaxes(
        showgrid=False,
        showline=True, linewidth=2, linecolor='black',
        title=dict(text=key, font=dict(size=14, color='DarkSlateGrey'))
    )
    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickangle=-45 if len(fig.data) > 6 else 0
    )
