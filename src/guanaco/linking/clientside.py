"""Browser-only highlight for an already-rendered pair of cell embeddings."""

CELL_HIGHLIGHT_JS = """
function(click, selection, figure, current) {
    const ctx = window.dash_clientside.callback_context;
    const events = (ctx.triggered || []).filter(item =>
        item.prop_id === '__SOURCE_ID__.clickData' ||
        item.prop_id === '__SOURCE_ID__.selectedData');
    let ids = current == null ? null : current;
    if (events.length) {
        const picked = new Set();
        for (const event of events) {
            for (const point of (event.value && event.value.points) || []) {
                // Embedding traces carry stable obs_names as Plotly ids, not
                // row positions: the two modalities may have different orders.
                if (point.id != null) picked.add(String(point.id));
            }
        }
        ids = picked.size ? Array.from(picked) : null;
    }
    const wrap = document.getElementById('__TARGET_ID__');
    const gd = wrap && (wrap.classList.contains('js-plotly-plot') ? wrap :
        wrap.querySelector('.js-plotly-plot'));
    if (gd && gd.data && window.Plotly) {
        const wanted = ids === null ? null : new Set(ids);
        const indices = [], selected = [];
        gd.data.forEach((trace, t) => {
            if (trace.type !== 'scattergl' || !trace.ids) return;
            indices.push(t);
            const matches = [];
            if (wanted) trace.ids.forEach((id, i) => {
                if (wanted.has(String(id))) matches.push(i);
            });
            selected.push(wanted === null ? null : matches);
        });
        if (indices.length) window.Plotly.restyle(gd, {
            selectedpoints: selected,
            'selected.marker.opacity': 1,
            'unselected.marker.opacity': 0.14,
            'unselected.marker.color': '#D1D5DB'
        }, indices);
    }
    return ids;
}
"""
