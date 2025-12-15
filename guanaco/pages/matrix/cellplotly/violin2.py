import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from statsmodels.formula.api import ols, mixedlm
import warnings

default_color = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"]

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


def calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels=None, design_type='crossed'):

    p_values = {}
    
    if mode == 'mode1':  # One metadata only
        # Simple comparison across meta1 groups
        groups = []
        group_names = []
        
        for group in (labels if labels else sorted(df[meta1].unique())):
            group_data = df[df[meta1] == group]['Expression'].values
            if len(group_data) > 0:
                groups.append(group_data)
                group_names.append(group)
        
        if len(groups) >= 2:
            if test_method == 'mwu-test' and len(groups) == 2:
                _, p_val = mannwhitneyu(groups[0], groups[1], alternative='two-sided')
                p_values['overall'] = p_val
            elif test_method == 'ttest' and len(groups) == 2:
                _, p_val = ttest_ind(groups[0], groups[1], equal_var=False)
                p_values['overall'] = p_val
            elif test_method in ['kw-test', 'anova']:
                test_func = kruskal if test_method == 'kw-test' else f_oneway
                _, p_val = test_func(*groups)
                p_values['overall'] = p_val
    
    elif mode == 'mode2':  # Facet by meta1, compare meta2
        meta1_groups = labels if labels else sorted(df[meta1].unique())
        
        if design_type == 'nested':

            for m1_group in meta1_groups:
                facet_data = df[df[meta1] == m1_group]
                groups = []
                group_names = []
                
                for m2_group in sorted(facet_data[meta2].unique()):
                    group_data = facet_data[facet_data[meta2] == m2_group]['Expression'].values
                    if len(group_data) > 0:
                        groups.append(group_data)
                        group_names.append(m2_group)
                
                if len(groups) >= 2:
                    if test_method in ['mwu-test', 'ttest'] and len(groups) == 2:
                        test_func = mannwhitneyu if test_method == 'mwu-test' else ttest_ind
                        kwargs = {'alternative': 'two-sided'} if test_method == 'mwu-test' else {'equal_var': False}
                        try:
                            _, p_val = test_func(groups[0], groups[1], **kwargs)
                            p_values[m1_group] = p_val
                        except:
                            p_values[m1_group] = np.nan
                            
                    elif test_method in ['kw-test', 'anova'] and len(groups) > 1:
                        test_func = kruskal if test_method == 'kw-test' else f_oneway
                        try:
                            _, p_val = test_func(*groups)
                            p_values[m1_group] = p_val
                        except:
                            p_values[m1_group] = np.nan
        else:
            # Compare meta2 within each meta1 facet (crossed/sparse designs)
            for m1_group in meta1_groups:
                facet_data = df[df[meta1] == m1_group]
                groups = []
                group_names = []
                
                for m2_group in sorted(facet_data[meta2].unique()):
                    group_data = facet_data[facet_data[meta2] == m2_group]['Expression'].values
                    if len(group_data) > 0:
                        groups.append(group_data)
                        group_names.append(m2_group)
                
                if len(groups) >= 2:
                    if test_method in ['mwu-test', 'ttest'] and len(groups) == 2:
                        test_func = mannwhitneyu if test_method == 'mwu-test' else ttest_ind
                        kwargs = {'alternative': 'two-sided'} if test_method == 'mwu-test' else {'equal_var': False}
                        try:
                            _, p_val = test_func(groups[0], groups[1], **kwargs)
                            p_values[m1_group] = p_val
                        except:
                            p_values[m1_group] = np.nan
                            
                    elif test_method in ['kw-test', 'anova'] and len(groups) > 1:
                        test_func = kruskal if test_method == 'kw-test' else f_oneway
                        try:
                            _, p_val = test_func(*groups)
                            p_values[m1_group] = p_val
                        except:
                            p_values[m1_group] = np.nan
    
    elif mode == 'mode3':  # Linear model
        if test_method == 'linear-model-interaction':
            # Linear model with interaction (Two-way ANOVA): expression ~ meta1 * meta2
            try:
                df_model = df.copy()
                df_model[meta1] = pd.Categorical(df_model[meta1])
                df_model[meta2] = pd.Categorical(df_model[meta2])
                
                formula = f'Expression ~ C({meta1}) * C({meta2})'
                model = ols(formula, data=df_model).fit()
                
                from statsmodels.stats.anova import anova_lm
                anova_table = anova_lm(model, typ=2)
                
                p_values['model_summary'] = {
                    'meta1_p': anova_table.loc[f'C({meta1})', 'PR(>F)'],
                    'meta2_p': anova_table.loc[f'C({meta2})', 'PR(>F)'],
                    'interaction_p': anova_table.loc[f'C({meta1}):C({meta2})', 'PR(>F)'],
                    'r_squared': model.rsquared,
                    'aic': model.aic
                }
                
                p_values['overall'] = p_values['model_summary']['interaction_p']
                
            except Exception as e:
                p_values['error'] = str(e)
                p_values['overall'] = np.nan
        else:
            # Default linear model: expression ~ meta1 + meta2
            try:
                df_model = df.copy()
                df_model[meta1] = pd.Categorical(df_model[meta1])
                df_model[meta2] = pd.Categorical(df_model[meta2])
                
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
    
    elif mode == 'mode4':  # Mixed model
        # Mixed model: expression ~ meta1 + (1|meta2)
        try:
            # Ensure categorical variables
            df_model = df.copy()
            df_model[meta1] = pd.Categorical(df_model[meta1])
            df_model[meta2] = pd.Categorical(df_model[meta2])
            
            # Try mixed effects model
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    formula = f'Expression ~ C({meta1})'
                    model = mixedlm(formula, data=df_model, groups=df_model[meta2]).fit(reml=False)
                    
                    # Extract p-values
                    p_values['model_summary'] = {
                        'meta1_p': model.pvalues[f'C({meta1})[T.{df_model[meta1].cat.categories[1]}]'],
                        'log_likelihood': model.llf,
                        'converged': model.converged
                    }
                    
                    p_values['overall'] = p_values['model_summary']['meta1_p']
                    
            except:
                # Fallback to pseudobulk approach
                pseudobulk = df.groupby([meta1, meta2])['Expression'].mean().reset_index()
                groups = []
                for m1 in sorted(pseudobulk[meta1].unique()):
                    group_data = pseudobulk[pseudobulk[meta1] == m1]['Expression'].values
                    if len(group_data) > 0:
                        groups.append(group_data)
                
                if len(groups) == 2:
                    _, p_val = ttest_ind(groups[0], groups[1], equal_var=False)
                    p_values['overall'] = p_val
                    p_values['method'] = 'pseudobulk_ttest'
                else:
                    _, p_val = f_oneway(*groups)
                    p_values['overall'] = p_val
                    p_values['method'] = 'pseudobulk_anova'
                    
        except Exception as e:
            p_values['error'] = str(e)
            p_values['overall'] = np.nan
    
    return p_values


def assign_colors(levels, color_map=None):
    if color_map is None:
        color_map = {}
    default_colors = px.colors.qualitative.Plotly
    for i, level in enumerate(sorted(levels)):
        if level not in color_map:
            color_map[level] = default_colors[i % len(default_colors)]
    return color_map


def add_p_value_annotations_new(fig, p_values, df, mode, meta1=None, meta2=None, split_violin=False, design_type='crossed'):
    if not p_values:
        return
    
    # Safe y_max calculation
    if len(df['Expression']) > 0:
        y_max = df['Expression'].max() * 1.1
    else:
        y_max = 1
    
    if mode == 'mode1':
        # Single p-value annotation
        if 'overall' in p_values and not np.isnan(p_values['overall']):
            p_val = p_values['overall']
            fig.add_annotation(
                x=0.5,
                y=1.02,
                xref='paper',
                yref='paper',
                text=f"<b>p = {p_val:.3g}</b>" if p_val < 0.05 else f"p = {p_val:.3g}",
                showarrow=False,
                font=dict(size=12, color='DarkSlateGrey')
            )
    
    elif mode == 'mode2':
        # P-values above each facet group
        meta1_groups = sorted(df[meta1].unique())
        
        if design_type == 'nested':

            meta2_groups = sorted(df[meta2].unique())
            
            meta2_to_meta1 = df.groupby(meta2)[meta1].first().to_dict()
            
            sorted_meta2 = sorted(meta2_groups, key=lambda x: (
                meta1_groups.index(meta2_to_meta1.get(x, meta1_groups[0])) if meta2_to_meta1.get(x) in meta1_groups else 999,
                str(x)
            ))
            
            # Find center position for each meta1 group
            for m1_group in meta1_groups:
                if m1_group not in p_values:
                    continue
                p_val = p_values[m1_group]
                if np.isnan(p_val):
                    continue
                
                # Find all meta2 that belong to this meta1
                m2_in_group = [m2 for m2 in sorted_meta2 if meta2_to_meta1.get(m2) == m1_group]
                
                if len(m2_in_group) == 0:
                    continue
                
                # Get the center meta2 name as x position
                center_idx = len(m2_in_group) // 2
                center_m2 = m2_in_group[center_idx]
                
                fig.add_annotation(
                    x=center_m2,
                    y=y_max,
                    text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                    showarrow=False,
                    font=dict(size=10, color='DarkSlateGrey')
                )
        
        elif split_violin:
            # For split violins, x position is the facet name directly
            for facet, p_val in p_values.items():
                if facet not in ['error', 'model_summary'] and not np.isnan(p_val) and facet in meta1_groups:
                    fig.add_annotation(
                        x=facet,
                        y=y_max,
                        text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                        showarrow=False,
                        font=dict(size=10, color='DarkSlateGrey')
                    )
        
        elif design_type in ['crossed', 'sparse']:
            # For crossed/sparse designs with combined x-labels (day1_patient1, day1_patient2, ...)
            # Position p-value above the center of each meta1 group
            meta2_groups = sorted(df[meta2].unique())
            n_meta2 = len(meta2_groups)
            
            for facet, p_val in p_values.items():
                if facet not in ['error', 'model_summary'] and not np.isnan(p_val) and facet in meta1_groups:
                    # Find the center x-label for this meta1 group
                    center_idx = n_meta2 // 2
                    center_m2 = meta2_groups[center_idx]
                    center_x_label = f"{facet}_{center_m2}"
                    
                    fig.add_annotation(
                        x=center_x_label,
                        y=y_max,
                        text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                        showarrow=False,
                        font=dict(size=10, color='DarkSlateGrey')
                    )
        else:
            # For grouped violins (sparse design), need to find center position of each group
            meta2_groups = sorted(df[meta2].unique())
            n_meta2 = len(meta2_groups)
            
            for i, (facet, p_val) in enumerate(p_values.items()):
                if facet not in ['error', 'model_summary'] and not np.isnan(p_val) and facet in meta1_groups:
                    # Calculate center position for this facet group
                    facet_idx = meta1_groups.index(facet)
                    start_pos = facet_idx * n_meta2
                    center_pos = start_pos + (n_meta2 - 1) / 2
                    
                    fig.add_annotation(
                        x=center_pos,
                        y=y_max,
                        text=f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                        showarrow=False,
                        font=dict(size=10, color='DarkSlateGrey'),
                        xref='x'
                    )
    
    elif mode in ['mode3', 'mode4']:
        # Single line p-value display for models
        if 'model_summary' in p_values:
            summary = p_values['model_summary']
            p_text_parts = []
            
            # Add meta1 p-value
            if 'meta1_p' in summary:
                p1 = summary['meta1_p']
                p1_text = f"{meta1}: <b>{p1:.3g}</b>" if p1 < 0.05 else f"{meta1}: {p1:.3g}"
                p_text_parts.append(p1_text)
            
            # Add meta2 p-value for linear model
            if mode == 'mode3' and 'meta2_p' in summary and not np.isnan(summary['meta2_p']):
                p2 = summary['meta2_p']
                p2_text = f"{meta2}: <b>{p2:.3g}</b>" if p2 < 0.05 else f"{meta2}: {p2:.3g}"
                p_text_parts.append(p2_text)
            
            # Add interaction p-value if present
            if 'interaction_p' in summary and not np.isnan(summary['interaction_p']):
                p_int = summary['interaction_p']
                p_int_text = f"interaction: <b>{p_int:.3g}</b>" if p_int < 0.05 else f"interaction: {p_int:.3g}"
                p_text_parts.append(p_int_text)
            
            # Join with comma
            p_text = ", ".join(p_text_parts)
            
            fig.add_annotation(
                x=0.5,
                y=1.02,
                xref='paper',
                yref='paper',
                text=p_text,
                showarrow=False,
                font=dict(size=12, color='DarkSlateGrey')
            )
        elif 'overall' in p_values:
            # Fallback for simple p-value (e.g., pseudobulk)
            p_val = p_values['overall']
            text = f"<b>{p_val:.3g}</b>" if p_val < 0.05 else f"{p_val:.3g}"
            
            fig.add_annotation(
                x=0.5,
                y=1.02,
                xref='paper',
                yref='paper',
                text=text,
                showarrow=False,
                font=dict(size=12, color='DarkSlateGrey')
            )


def plot_violin2_new(adata, key, meta1, meta2, mode, transformation='log',
                     show_box=False, show_points=False, test_method='auto', 
                     labels=None, color_map=None):
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
    
    # Extract and transform expression data
    from .gene_extraction_utils import extract_gene_expression, apply_transformation
    expression_data = extract_gene_expression(adata, key)
    
    if transformation:
        expression_data = apply_transformation(expression_data, transformation, copy=False)
    
    # Build dataframe
    df = pd.DataFrame({
        'Expression': expression_data,
        meta1: adata.obs[meta1].values
    })
    
    if mode != 'mode1' and meta2:
        df[meta2] = adata.obs[meta2].values
    
    # Filter by selected labels if provided
    if labels:
        df = df[df[meta1].isin(labels)]
    
    # Handle empty dataframe
    if len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available for selected filters",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color='gray')
        )
        return fig
    
    # Get unique categories (preserve order if categorical)
    meta1_cats = _get_ordered_categories(df, meta1)
    meta2_cats = _get_ordered_categories(df, meta2) if mode != 'mode1' and meta2 else []
    
    # Detect design pattern to determine plotting strategy
    design_type = _detect_design_pattern(df, meta1, meta2) if meta2 and mode != 'mode1' else 'single'
    
    # Assign colors based on design type
    # For nested designs, color by meta1 (the parent category)
    # For crossed/sparse designs, color by meta2
    if design_type == 'nested':
        color_map = assign_colors(meta1_cats, color_map)
    elif meta2_cats:
        color_map = assign_colors(meta2_cats, color_map)
    
    # Determine test method
    meta1_levels = len(meta1_cats)
    meta2_levels = len(meta2_cats)
    
    if test_method == 'auto':
        test_method, test_description = determine_test_method(meta1_levels, meta2_levels, mode)
    else:
        test_description = ""
    
    # Create figure based on mode
    fig = go.Figure()
    points_mode = 'all' if show_points else False
    
    if mode == 'mode1':
        # Single metadata violin plot
        fig, title, xaxis_title = _create_single_meta_violin(
            fig, df, key, meta1, meta1_cats, 
            show_box, points_mode
        )
    
    elif mode == 'mode2':
        # Two metadata comparison
        if design_type == 'nested':
            # Nested design: use meta2 as x-axis, color by meta1
            fig, title, xaxis_title = _create_nested_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
        elif design_type in ['crossed', 'sparse']:
            # Crossed/Sparse design: each combination gets its own x position
            fig, title, xaxis_title = _create_crossed_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
        elif len(meta2_cats) == 2:
            # Split violin for binary comparison
            fig, title, xaxis_title = _create_split_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
        else:
            # Fallback to crossed violin
            fig, title, xaxis_title = _create_crossed_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
    
    elif mode in ['mode3', 'mode4']:
        if design_type == 'nested':
            fig, title, xaxis_title = _create_nested_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
        else:
            # Use crossed violin for crossed/sparse designs
            fig, title, xaxis_title = _create_crossed_violin(
                fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                color_map, show_box, points_mode
            )
    
    # Calculate and add p-values
    if test_method and test_method != 'none':
        p_values = calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels, design_type)
        split_violin = (mode == 'mode2' and meta2 and len(meta2_cats) == 2 and design_type != 'nested')
        add_p_value_annotations_new(fig, p_values, df, mode, meta1, meta2, split_violin, design_type)
    
    _apply_common_layout(fig, df, key, title, xaxis_title, meta2, mode, show_points)
    
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


def _create_single_meta_violin(fig, df, key, meta1, meta1_cats, show_box, points_mode):
    """Create violin plot with single metadata (mode1)."""
    
    for i, group in enumerate(meta1_cats):
        group_data = df[df[meta1] == group]['Expression']
        
        if len(group_data) == 0:
            continue
        
        spanmode = _get_violin_spanmode(group_data)
        
        fig.add_trace(go.Violin(
            x=[group] * len(group_data),
            y=group_data,
            name=str(group),
            box_visible=show_box,
            points=points_mode,
            meanline_visible=True,
            width=0.8,
            spanmode=spanmode,
            fillcolor=default_color[i % len(default_color)],
            line_color='DarkSlateGrey',
            hoveron='violins',
            jitter=0.05,
            showlegend=False
        ))
    
    title = f"Expression of {key} by {meta1}"
    xaxis_title = meta1
    
    return fig, title, xaxis_title


def _create_split_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                         color_map, show_box, points_mode):
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
                box_visible=show_box,
                points=points_mode,
                meanline_visible=True,
                spanmode=spanmode,
                fillcolor=color_map.get(m2_group, default_color[j % len(default_color)]),
                line_color='DarkSlateGrey',
                hoveron='violins',
                jitter=0.05,
                showlegend=show_legend
            ))
    
    fig.update_layout(violinmode='overlay')
    
    title = f"Expression of {key} by {meta1} and {meta2}"
    xaxis_title = meta1
    
    return fig, title, xaxis_title


def _create_nested_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                          color_map, show_box, points_mode):
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
            box_visible=show_box,
            points=points_mode,
            meanline_visible=True,
            spanmode=spanmode,
            fillcolor=color_map.get(m1_group, default_color[m1_idx % len(default_color)]),
            line_color='DarkSlateGrey',
            hoveron='violins',
            hoverinfo='y+name',
            jitter=0.05,
            showlegend=show_legend
        ))
    
    # No grouping needed - each meta2 has its own x position
    fig.update_layout(
        violingap=0.1
    )
    
    title = f"Expression of {key} by {meta2} (colored by {meta1})"
    xaxis_title = meta2
    
    return fig, title, xaxis_title


def _create_crossed_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                           color_map, show_box, points_mode):
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
                box_visible=show_box,
                points=points_mode,
                meanline_visible=True,
                spanmode=spanmode,
                fillcolor=color_map.get(m2_group, default_color[j % len(default_color)]),
                line_color='DarkSlateGrey',
                hoveron='violins',
                hoverinfo='y+name',
                jitter=0.05,
                showlegend=show_legend
            ))
    
    # Set x-axis category order to maintain grouping
    fig.update_xaxes(categoryorder='array', categoryarray=x_labels)
    
    fig.update_layout(
        violingap=0.1
    )
    
    title = f"Expression of {key} by {meta1} and {meta2}"
    xaxis_title = f"{meta1} × {meta2}"
    
    return fig, title, xaxis_title


def _create_grouped_violin(fig, df, key, meta1, meta2, meta1_cats, meta2_cats,
                           color_map, show_box, points_mode):
    """
    Create grouped violin plot for multiple meta2 categories.
    
    Handles: Crossed, Nested, and Sparse designs by checking data existence.
    Uses violinmode='group' for side-by-side placement.
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
            
            spanmode = _get_violin_spanmode(subset['Expression'])
            
            # Only show legend for first occurrence of each meta2 group
            show_legend = m2_group not in legend_added
            if show_legend:
                legend_added.add(m2_group)
            
            fig.add_trace(go.Violin(
                x=[m1_group] * len(subset),
                y=subset['Expression'],
                name=str(m2_group),
                legendgroup=str(m2_group),  # Group legend items for toggle
                width=0.8,
                box_visible=show_box,
                points=points_mode,
                meanline_visible=True,
                spanmode=spanmode,
                fillcolor=color_map.get(m2_group, default_color[j % len(default_color)]),
                line_color='DarkSlateGrey',
                hoveron='violins',
                hoverinfo='y+name',
                jitter=0.05,
                showlegend=show_legend
            ))
    
    # Use 'group' mode for side-by-side violins
    fig.update_layout(
        violinmode='group',
        violingap=0.1,
        violingroupgap=0.05
    )
    
    title = f"Expression of {key} by {meta1} and {meta2}"
    xaxis_title = meta1
    
    return fig, title, xaxis_title


def _apply_common_layout(fig, df, key, title, xaxis_title, meta2, mode, show_points):
    """Apply common layout settings to the figure."""
    
    # Calculate y-axis range safely
    if len(df['Expression']) > 0:
        y_max = df['Expression'].max() * 1.15
        y_range = [0, y_max]
    else:
        y_range = [0, 1]
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12, color='DarkSlateGrey'),
        title=dict(
            text=title,
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
    
    if show_points:
        fig.update_traces(marker=dict(size=1.5), selector=dict(type='violin'))
    
    fig.update_yaxes(
        showgrid=False,
        showline=True, linewidth=2, linecolor='black',
        title=dict(text=key, font=dict(size=14, color='DarkSlateGrey'))
    )
    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black',
        tickangle=-45 if len(fig.data) > 6 else 0
    )