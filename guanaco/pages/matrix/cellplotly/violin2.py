import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
import statsmodels.api as sm
from statsmodels.formula.api import ols, mixedlm
from guanaco.data_loader import color_config
from itertools import combinations
import warnings

def determine_test_method(meta1_levels, meta2_levels, mode, test_override=None):
    """
    Automatically determine the appropriate statistical test based on metadata and mode.
    
    Returns: (test_method, test_description)
    """
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


def calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels=None):
    """
    Calculate p-values based on the selected mode and test method.
    """
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
        # Compare meta2 within each meta1 facet
        meta1_groups = labels if labels else sorted(df[meta1].unique())
        
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
        # Linear model: expression ~ meta1 + meta2
        try:
            # Ensure categorical variables
            df_model = df.copy()
            df_model[meta1] = pd.Categorical(df_model[meta1])
            df_model[meta2] = pd.Categorical(df_model[meta2])
            
            formula = f'Expression ~ C({meta1}) + C({meta2})'
            model = ols(formula, data=df_model).fit()
            
            # Extract p-values
            p_values['model_summary'] = {
                'meta1_p': model.pvalues[f'C({meta1})[T.{df_model[meta1].cat.categories[1]}]'],
                'meta2_p': model.pvalues[f'C({meta2})[T.{df_model[meta2].cat.categories[1]}]'] if len(df_model[meta2].unique()) > 1 else np.nan,
                'r_squared': model.rsquared,
                'aic': model.aic
            }
            
            # Store main effect p-value
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
    """Assign colors to categorical levels."""
    if color_map is None:
        color_map = {}
    default_colors = px.colors.qualitative.Plotly
    for i, level in enumerate(sorted(levels)):
        if level not in color_map:
            color_map[level] = color_config[i] if i < len(color_config) else default_colors[i % len(default_colors)]
    return color_map


def add_p_value_annotations_new(fig, p_values, df, mode, meta1=None, meta2=None, split_violin=False):
    """Add p-value annotations based on mode."""
    if not p_values:
        return
    
    y_max = max(df['Expression']) * 1.1  # Increased from 1.05 to 1.1 for more space
    
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
        
        if split_violin:
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
        else:
            # For grouped violins, need to find center position of each group
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
                     show_box=True, show_points=True, test_method='auto', 
                     labels=None, color_map=None):
    """
    Plot violin plots with new mode-based analysis.
    
    Parameters:
    -----------
    mode : str
        'mode1': One metadata only (meta2 is None)
        'mode2': Facet by meta1, compare meta2
        'mode3': Linear model with meta1 + meta2
        'mode4': Mixed model with meta1 + (1|meta2)
    """
    if not key or not meta1:
        return go.Figure()
    
    # Extract expression data
    if key not in adata.var_names:
        return go.Figure()
    
    # Use gene extraction utility to handle views properly
    from .gene_extraction_utils import extract_gene_expression, apply_transformation
    expression_data = extract_gene_expression(adata, key)
    
    # Apply transformation if specified
    if transformation:
        expression_data = apply_transformation(expression_data, transformation, copy=False)
    
    # Prepare dataframe
    df = pd.DataFrame({
        'Expression': expression_data,
        meta1: adata.obs[meta1]
    })
    
    if mode != 'mode1' and meta2:
        df[meta2] = adata.obs[meta2]
    
    # Filter by selected labels if provided
    if labels:
        df = df[df[meta1].isin(labels)]
    
    # Count levels
    meta1_levels = df[meta1].nunique()
    meta2_levels = df[meta2].nunique() if mode != 'mode1' and meta2 else 0
    
    # Determine test method if auto
    if test_method == 'auto':
        test_method, test_description = determine_test_method(meta1_levels, meta2_levels, mode)
    else:
        test_description = ""
    
    # Create figure
    fig = go.Figure()
    points_mode = 'all' if show_points else False
    
    if mode == 'mode1':
        # Simple violin plot with meta1 only
        meta1_groups = sorted(df[meta1].unique())
        
        for i, group in enumerate(meta1_groups):
            group_data = df[df[meta1] == group]['Expression']
            
            # Calculate variance to determine bandwidth (seaborn-like behavior)
            variance = np.var(group_data)
            if variance < 1e-10:  # Near-zero variance
                bandwidth = 0.01
                spanmode = 'soft'
            else:
                # Normal bandwidth for data with variance
                bandwidth = 0.2
                spanmode = 'hard'
            
            fig.add_trace(go.Violin(
                x=[group] * len(group_data),
                y=group_data,
                name=group,
                box_visible=show_box,
                points=points_mode,
                meanline_visible=True,
                width=0.8,
                bandwidth=bandwidth,
                spanmode=spanmode,
                fillcolor=color_config[i % len(color_config)],
                line_color='DarkSlateGrey',
                hoveron='violins',
                jitter=0.05,
                showlegend=False
            ))
        
        title = f"Expression of {key} by {meta1}"
        xaxis_title = meta1
    
    elif mode == 'mode2':
        # Faceted plot: group by meta1, split/group by meta2
        meta1_groups = sorted(df[meta1].unique())
        meta2_groups = sorted(df[meta2].unique())
        color_map = assign_colors(meta2_groups, color_map)
        
        if len(meta2_groups) == 2:
            # Use split violins for binary comparison
            for i, m1_group in enumerate(meta1_groups):
                facet_data = df[df[meta1] == m1_group]
                
                for j, m2_group in enumerate(meta2_groups):
                    side = 'negative' if j == 0 else 'positive'
                    subset = facet_data[facet_data[meta2] == m2_group]
                    
                    if len(subset) > 0:
                        # Calculate variance to determine bandwidth (seaborn-like behavior)
                        variance = np.var(subset['Expression'])
                        if variance < 1e-10:  # Near-zero variance
                            bandwidth = 0.01
                            spanmode = 'soft'
                        else:
                            # Normal bandwidth for data with variance
                            bandwidth = 0.2
                            spanmode = 'hard'
                        
                        fig.add_trace(go.Violin(
                            x=[m1_group] * len(subset),
                            y=subset['Expression'],
                            legendgroup=m2_group,
                            width=0.8,
                            scalemode='width',
                            scalegroup=str(m1_group),
                            name=m2_group,
                            side=side,
                            box_visible=show_box,
                            points=points_mode,
                            meanline_visible=True,
                            bandwidth=bandwidth,
                            spanmode=spanmode,
                            fillcolor=color_map[m2_group],
                            line_color='DarkSlateGrey',
                            hoveron='violins',
                            showlegend=(i == 0),
                            jitter=0.05
                        ))
            
            fig.update_layout(violinmode='overlay')
        else:
            # Use grouped violins for multiple groups
            x_labels = []
            for m1_group in meta1_groups:
                for m2_group in meta2_groups:
                    x_labels.append(f"{m1_group}<br>{m2_group}")
            
            x_pos = 0
            for m1_group in meta1_groups:
                facet_data = df[df[meta1] == m1_group]
                
                for m2_group in meta2_groups:
                    subset = facet_data[facet_data[meta2] == m2_group]
                    
                    if len(subset) > 0:
                        show_legend = m2_group not in [t.legendgroup for t in fig.data]
                        
                        # Calculate variance to determine bandwidth (seaborn-like behavior)
                        variance = np.var(subset['Expression'])
                        if variance < 1e-10:  # Near-zero variance
                            bandwidth = 0.01
                            spanmode = 'soft'
                        else:
                            # Normal bandwidth for data with variance
                            bandwidth = 0.2
                            spanmode = 'hard'
                        
                        fig.add_trace(go.Violin(
                            x=[x_labels[x_pos]] * len(subset),
                            y=subset['Expression'],
                            name=m2_group,
                            legendgroup=m2_group,
                            width=0.8,
                            box_visible=show_box,
                            points=points_mode,
                            meanline_visible=True,
                            bandwidth=bandwidth,
                            fillcolor=color_map[m2_group],
                            line_color='DarkSlateGrey',
                            hoveron='violins',
                            hoverinfo='y',
                            spanmode=spanmode,
                            jitter=0.05,
                            showlegend=show_legend
                        ))
                    
                    x_pos += 1
            
            fig.update_layout(violingap=0.3, violingroupgap=0.1)
        
        title = f"Expression of {key} by {meta1} and {meta2}"
        xaxis_title = meta1 if len(meta2_groups) == 2 else ""
    
    elif mode in ['mode3', 'mode4']:
        # Same visualization as mode2 but with different statistical analysis
        # The visualization logic is identical to mode2
        meta1_groups = sorted(df[meta1].unique())
        meta2_groups = sorted(df[meta2].unique())
        color_map = assign_colors(meta2_groups, color_map)
        
        # Create faceted visualization
        x_labels = []
        for m1_group in meta1_groups:
            for m2_group in meta2_groups:
                x_labels.append(f"{m1_group}<br>{m2_group}")
        
        x_pos = 0
        for m1_group in meta1_groups:
            facet_data = df[df[meta1] == m1_group]
            
            for m2_group in meta2_groups:
                subset = facet_data[facet_data[meta2] == m2_group]
                
                if len(subset) > 0:
                    show_legend = m2_group not in [t.legendgroup for t in fig.data]
                    
                    # Calculate variance to determine bandwidth (seaborn-like behavior)
                    variance = np.var(subset['Expression'])
                    if variance < 1e-10:  # Near-zero variance
                        bandwidth = 0.01
                        spanmode = 'soft'
                    else:
                        # Normal bandwidth for data with variance
                        bandwidth = 0.2
                        spanmode = 'hard'
                    
                    fig.add_trace(go.Violin(
                        x=[x_labels[x_pos]] * len(subset),
                        y=subset['Expression'],
                        name=m2_group,
                        legendgroup=m2_group,
                        width=0.8,
                        box_visible=show_box,
                        points=points_mode,
                        meanline_visible=True,
                        bandwidth=bandwidth,
                        fillcolor=color_map[m2_group],
                        line_color='DarkSlateGrey',
                        hoveron='violins',
                        hoverinfo='y',
                        spanmode=spanmode,
                        jitter=0.05,
                        showlegend=show_legend
                    ))
                
                x_pos += 1
        
        fig.update_layout(violingap=0.3, violingroupgap=0.1)
        
        title = f"Expression of {key} by {meta1} and {meta2}"
        xaxis_title = ""
    
    # Calculate p-values
    if test_method and test_method != 'none':
        p_values = calculate_p_values_by_mode(df, meta1, meta2, mode, test_method, labels)
        split_violin = (mode == 'mode2' and meta2 and df[meta2].nunique() == 2)
        add_p_value_annotations_new(fig, p_values, df, mode, meta1, meta2, split_violin)
    
    # Shared layout updates
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12, color='DarkSlateGrey'),
        title=dict(
            text=title,
            x=0.5,
            xanchor='center',
            y=0.9,
            font=dict(size=16, color='DarkSlateGrey')
        ),
        yaxis_range=[0, max(df['Expression']) * 1.1],
        xaxis_title=xaxis_title,
        yaxis_title='Expression',
        legend=dict(
            font=dict(size=12),
            title=dict(text=meta2 if mode != 'mode1' else "", font=dict(size=12))
        ),
        margin=dict(b=80, t=80)
    )
    
    if show_points:
        fig.update_traces(marker=dict(size=1.5), selector=dict(type='violin'))
    
    fig.update_yaxes(
        showgrid=False,
        showline=True, linewidth=2, linecolor='black',
        title=dict(text=key, font=dict(size=14, color='DarkSlateGrey'))
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    
    return fig