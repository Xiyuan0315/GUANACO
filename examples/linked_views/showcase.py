"""Small public gallery using the existing LinkedView components and callbacks."""

import inspect
from pathlib import Path
from urllib.parse import parse_qs

from dash import Dash, Input, Output, dcc, html

from cases import CASES


def example_page(case, demo):
    # Capture each pristine layout once, before any visitors make selections.
    # Never rebuild from demo._state: that is the notebook's last-selection mirror,
    # not a browser session. Each mounted dcc.Store gets its own initial state.
    return html.Article([
        html.Header([
            html.Span("INTERACTIVE EXAMPLE", className="eyebrow"),
            html.H1(case.title),
            html.P(case.question, className="question"),
        ]),
        html.Div([html.Strong("Try it"), html.Span(case.instruction)], className="try-it"),
        html.Div(demo.component(), className="linked-example"),
        html.Section([
            html.H2("What is linked?"), html.P(case.explanation),
            html.Details([
                html.Summary("Show Python code"),
                html.P("This is the same builder used by the live example. Import it from cases.py to use it in a notebook."),
                dcc.Markdown(f"```python\nfrom cases import {case.build.__name__}\n\n"
                             f"demo = {case.build.__name__}()\ndemo.show_jupyter()\n```"),
                dcc.Markdown(f"```python\n{inspect.getsource(case.build)}\n```"),
            ]),
            html.Details([
                html.Summary("About the demonstration data"),
                html.P(case.data_note),
                html.P([
                    "Spatial data: ",
                    html.A("Squidpy Visium H&E tutorial", href="https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html",
                           target="_blank", rel="noopener noreferrer"),
                    ". The other four examples use deterministic synthetic data. These are interaction demonstrations, not performance benchmarks.",
                ]),
            ]),
        ], className="example-notes"),
    ])


def create_app():
    app = Dash(
        __name__, suppress_callback_exceptions=True,
        assets_folder=str(Path(__file__).with_name("showcase_assets")),
        title="GUANACO · Interactive examples",
        update_title=None,
        meta_tags=[
            {"name": "viewport", "content": "width=device-width, initial-scale=1"},
            {"name": "description", "content": "Try five GUANACO examples connecting cell populations, RNA, protein, spatial relationships and chromatin accessibility."},
        ],
    )
    # Deliberately bounded: five small datasets, one mounted example per visitor.
    # Register every callback once at startup, never inside a navigation callback.
    pages = {}
    for case in CASES:
        demo = case.build()
        pages[case.slug] = example_page(case, demo)
        demo.register(app)

    app.layout = html.Div([
        dcc.Location(id="gallery-location", refresh=False),
        html.A("Skip to example", href="#example-content", className="skip-link"),
        html.Aside([
            html.Div([html.Strong("GUANACO"), html.Span("Linked-view gallery")], className="brand"),
            html.P("Follow a biological question across connected plots.", className="sidebar-intro"),
            html.Nav([
                dcc.Link([
                    html.Span(f"{index:02d}", className="nav-number"),
                    html.Span(case.title),
                ], href=f"?demo={case.slug}", id=f"nav-{case.slug}", className="demo-link")
                for index, case in enumerate(CASES, start=1)
            ], **{"aria-label": "Examples"}),
            html.Div([
                html.Span("Real spatial data · 4 synthetic demos", className="data-badge"),
                html.P("No uploads or installation needed. Click or select directly in the plots."),
                html.A("GUANACO on GitHub ↗", href="https://github.com/Systems-Immunometabolism-Lab/guanaco-viz",
                       target="_blank", rel="noopener noreferrer"),
            ], className="sidebar-footer"),
        ], className="sidebar"),
        html.Main([
            html.Div([
                html.Span("LIVE DEMO · 5 LINKED WORKFLOWS"),
                html.Button("Reset example", id="reset-example", n_clicks=0),
            ], className="toolbar"),
            dcc.Loading(
                html.Div(id="example-content", tabIndex=-1),
                id="example-loading", type="circle", color="#176B67",
                # Only navigation/reset replaces the page. Nested figure and
                # selection-store updates must not hide the entire example.
                target_components={"example-content": "children"},
                delay_show=250, show_initially=False,
            ),
        ], className="main-content"),
    ], className="gallery")

    @app.callback(
        Output("example-content", "children"),
        *[Output(f"nav-{case.slug}", "className") for case in CASES],
        Input("gallery-location", "search"), Input("reset-example", "n_clicks"),
    )
    def navigate(search, resets):
        slug = parse_qs((search or "").lstrip("?")).get("demo", [CASES[0].slug])[0]
        if slug not in pages:
            slug = CASES[0].slug
        # A React key forces a fresh mount on reset, clearing graphs AND stores.
        page = html.Div(pages[slug], key=f"{slug}-{resets or 0}")
        return page, *[
            "demo-link active" if case.slug == slug else "demo-link" for case in CASES
        ]

    return app
