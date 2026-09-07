/* Preset display only. No cross-panel events, data aggregation, or network I/O. */
(() => {
    "use strict";
    const panels = JSON.parse(document.getElementById("guanaco-static-data").textContent);
    const dashboard = document.getElementById("dashboard");
    const ready = [];
    function element(tag, text, className) {
        const node = document.createElement(tag);
        if (text != null) node.textContent = text;
        if (className) node.className = className;
        return node;
    }
    for (const panel of panels) {
        const card = element("section", null, "panel" + (panel.wide ? " wide" : ""));
        card.id = "panel-" + panel.id;
        const heading = element("div", null, "panel-header");
        heading.append(element("h2", panel.title), element("p", panel.description));
        const controls = element("div", null, "controls");
        const label = element("label", "Prepared grouping");
        const select = element("select");
        select.id = "preset-" + panel.id;
        label.htmlFor = select.id;
        panel.presets.forEach((preset, index) => {
            const option = element("option", preset.label);
            option.value = String(index);
            select.append(option);
        });
        const reset = element("button", "Reset view");
        reset.type = "button";
        reset.setAttribute("aria-label", "Reset " + panel.title);
        controls.append(label, select, reset);
        heading.append(controls);
        const plot = element("div", null, "plot");
        plot.id = "plot-" + panel.id;
        const plotFrame = element("div", null, "plot-frame");
        plotFrame.append(plot);
        const status = element("span", "Opening saved figure…", "panel-status");
        status.setAttribute("role", "status");
        card.append(heading, plotFrame, status);
        dashboard.append(card);
        async function render() {
            select.disabled = reset.disabled = true;
            const preset = panel.presets[Number(select.value)];
            // Plotly mutates inputs; keep the saved preset immutable for resets.
            const figure = JSON.parse(JSON.stringify(preset.figure));
            const layout = {...figure.layout, autosize: true};
            delete layout.width;
            delete layout.uirevision;
            layout.title = {text: ""};
            try {
                await Plotly.react(plot, figure.data, layout, {
                    responsive: true, displaylogo: false, scrollZoom: false,
                    modeBarButtonsToRemove: ["select2d", "lasso2d"],
                    toImageButtonOptions: {format: "png", filename: "guanaco-" + panel.id}
                });
                status.textContent = preset.label + " · saved result · no recalculation";
                status.classList.remove("error");
            } catch (error) {
                status.textContent = "Unable to display this preset: " + error.message;
                status.classList.add("error");
                throw error;
            } finally {
                select.disabled = reset.disabled = false;
            }
        }
        select.addEventListener("change", () => render().catch(console.error));
        reset.addEventListener("click", () => render().catch(console.error));
        ready.push(render());
    }
    Promise.all(ready).then(() => {
        document.documentElement.dataset.viewerReady = "true";
    }).catch(() => { document.documentElement.dataset.viewerReady = "error"; });
})();
