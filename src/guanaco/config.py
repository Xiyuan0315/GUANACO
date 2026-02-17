common_config = {
    'toImageButtonOptions': {
        'format': 'svg',
        'scale': 1,
        'filename': 'plot'
    },
    'displaylogo': False,
    'modeBarButtonsToRemove': ['select', 'lasso2d', 'pan2d', 'zoom2d', 'autoScale2d', 'zoomIn2d', 'zoomOut2d', 'pan', 'zoom', 'resetScale2d']
}

# Annotation scatter config (left side) - with selection tools
scatter_config = {
    'responsive': True,
    'toImageButtonOptions': {
        'format': 'svg',
        'scale': 1,
        'filename': 'annotation_scatter'
    },
    'displaylogo': False,
    'scrollZoom': False,  # Disable scroll to zoom
    'doubleClick': 'reset',  # Double click to reset zoom
    'showTips': False,  # Reduce hover updates
    'frameMargins': 0,  # Reduce margin calculations
    'modeBarButtonsToRemove': ['select',  'autoScale2d','zoomIn','zoomOut'],
}

# Gene scatter config (right side) - pan only, no selection tools
gene_scatter_config = {
    'responsive': True,
    'toImageButtonOptions': {
        'format': 'svg',
        'scale': 1,
        'filename': 'gene_scatter'
    },
    'displaylogo': False,
    'scrollZoom': False,  # Disable scroll to zoom
    'doubleClick': 'reset',  # Double click to reset zoom
    'showTips': False,  # Reduce hover updates
    'frameMargins': 0,  # Reduce margin calculations
    'modeBarButtonsToRemove': ['lasso2d','select',  'autoScale2d','zoomIn','zoomOut'],
}

# Debounced scatter config - use this for callbacks that need debouncing
scatter_config_debounced = {
    **scatter_config,
    'edits': {
        'shapePosition': True,
        'annotationPosition': True
    },
    'queueLength': 1  # Limit update queue to reduce lag
}

# Performance optimizations for scatter plots
def optimize_scatter_performance(fig):
    """Apply performance optimizations to scatter plot figures"""
    fig.update_layout(
        # Reduce update frequency
        updatemenus=None,
        # Set pan as default drag mode
        # Smooth zoom transitions
        transition={
            'duration': 150,  # Faster transitions
            'easing': 'cubic-in-out'
        },
        # Reduce reflow calculations
        autosize=True,
        # Optimize for large datasets
        hovermode='closest'  # Reduces hover calculations
    )
    
    # Optimize trace rendering
    for trace in fig.data:
        if hasattr(trace, 'mode') and 'markers' in str(trace.mode):
            trace.update(
                marker=dict(
                    line=dict(width=0),  # Remove marker borders to speed up rendering
                    sizemode='diameter'  # More efficient size mode
                )
            )
    
    return fig
