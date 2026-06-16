"""Desktop entry point for GUANACO.

Flow for the packaged app:

    1. Show the tkinter config wizard so the user can create / pick a config,
       while the heavy (config-independent) libraries import in the background.
    2. Open a native pywebview window with a loading splash.
    3. Load the data and start the Dash server on a background thread, then
       swap the splash for the live app once the server answers.

The data is loaded lazily inside the webview worker thread (not before the
window opens) so a double-clicked ``.app`` shows the splash within a second,
instead of appearing to hang while AnnData/MuData loads. The library prewarm in
step 1 means that by launch time only the data load remains.
"""

from __future__ import annotations

import argparse
import threading
from pathlib import Path

from guanaco.cli import (
    SPLASH_HTML,
    _load_settings,
    _wait_for_server,
    apply_settings_env,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GUANACO desktop app: create a config, then explore your data in a native window.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=None,
        help="Path to an existing config JSON. Pre-fills the config builder; with --no-builder it launches directly.",
    )
    parser.add_argument(
        "--no-builder",
        action="store_true",
        help="Skip the config builder and launch straight from --config (which must exist).",
    )
    return parser.parse_args()


def _prewarm_libraries() -> None:
    """Import the heavy, config-independent libraries in the background.

    The Dash + anndata/muon + matplotlib import cost (~30s, plus matplotlib's
    one-time font-cache build) does not depend on which config the user picks,
    so we pay it *while* they fill out the wizard. By the time they launch,
    these modules are already in ``sys.modules`` and ``guanaco.main`` imports
    them for free — only the data load remains. Best-effort: any failure here
    is swallowed, and the real import later surfaces the error properly.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # non-interactive; safe off the main thread
        import matplotlib.pyplot  # noqa: F401  (triggers the font-cache build)
        import anndata  # noqa: F401
        import muon  # noqa: F401
        import plotly  # noqa: F401
        import dash  # noqa: F401
        import dash_bootstrap_components  # noqa: F401
        import dash_bio  # noqa: F401  (largest single chunk)
        import dash_ag_grid  # noqa: F401
    except Exception:
        pass


def _start_prewarm() -> threading.Thread:
    thread = threading.Thread(target=_prewarm_libraries, daemon=True)
    thread.start()
    return thread


def _resolve_config_via_wizard(initial: Path | None) -> Path | None:
    """Run the tkinter wizard. Returns the chosen config path if the user asked
    to launch, or None if they only saved / closed the window."""
    from guanaco.config_wizard import launch_config_wizard

    default = initial or (Path.home() / "guanaco.json")
    result = launch_config_wizard(default)
    if result is None:
        return None
    config_path, should_launch = result
    return config_path if should_launch else None


def _run_window(settings: dict) -> None:
    """Open the pywebview window and bring the Dash server up behind a splash."""
    # Imported dynamically (not `import webview`) so static import scanners -- e.g.
    # the Dash cloud builder's dependency detection -- don't pick it up and try to
    # install the unrelated PyPI `webview` package, which needs native GTK/WebKit
    # and fails to build on a server. Desktop mode is the only path that needs it.
    import importlib
    try:
        webview = importlib.import_module("webview")
    except ImportError:
        raise SystemExit(
            "Desktop mode requires pywebview. Install it with:\n"
            "    pip install 'guanaco-viz[desktop]'"
        )

    # Off by default in pywebview, so the webview silently declines every
    # download (CSV exports, plot PNGs). Enabling it lets the native webview
    # save files to ~/Downloads.
    webview.settings["ALLOW_DOWNLOADS"] = True

    host = settings["host"]
    port = settings["port"]
    # pywebview must hit a loopback address; 0.0.0.0 is only a bind address.
    window_host = "127.0.0.1" if host in ("0.0.0.0", "", None) else host
    url = f"http://{window_host}:{port}"

    window = webview.create_window(
        "GUANACO",
        html=SPLASH_HTML,
        width=1440,
        height=900,
        min_size=(1024, 700),
        text_select=True,
    )

    def boot() -> None:
        # Imported here (not at module scope) so the splash is already on screen
        # while the data loads — registry.initialize_data() runs at import time.
        try:
            from guanaco.main import app
        except Exception as exc:  # surface load errors inside the window
            import html as _html
            window.load_html(
                "<div style='font-family:sans-serif;padding:2rem;color:#b00'>"
                f"<h2>GUANACO failed to load</h2><pre>{_html.escape(str(exc))}</pre></div>"
            )
            return

        # Optional public sharing. Only tunnelled requests are challenged, so the
        # local webview below loads without a password prompt.
        if settings.get("share"):
            from guanaco import share

            public_url, username, password = share.enable_sharing(
                app.server, port, settings["share_username"], settings["share_password"]
            )
            share.print_share_banner(public_url, username, password)

        threading.Thread(
            target=lambda: app.run_server(
                host=host, port=port, debug=False, use_reloader=False
            ),
            daemon=True,
        ).start()

        # Data is already loaded by this point, so the server comes up quickly.
        if _wait_for_server(url, timeout=60):
            try:
                window.set_title(getattr(app, "title", None) or "GUANACO")
            except Exception:
                pass
            window.load_url(url)
        else:
            window.load_html(
                "<h2 style='font-family:sans-serif;padding:2rem'>"
                "GUANACO failed to start the server within 60 seconds.</h2>"
            )

    # gui=None lets pywebview pick the native backend (Cocoa / EdgeChromium / GTK).
    webview.start(boot)


def desktop_main() -> None:
    """Entry point: wizard → desktop window."""
    import gc

    args = _parse_args()

    config_path = None
    if args.config is not None:
        config_path = args.config.expanduser()
        if not config_path.is_absolute():
            config_path = config_path.resolve()

    launch_directly = args.no_builder and config_path is not None and config_path.exists()

    if not launch_directly:
        # The wizard creates a Tk interpreter on this (main) thread, and Tk
        # requires it to be torn down on the same thread. The prewarm/boot
        # worker threads run cyclic GC during their heavy imports, which would
        # otherwise finalize the wizard's Tk objects (kept alive only by
        # reference cycles) off-thread and abort with a Tcl_AsyncDelete panic.
        # So: suspend automatic GC across the wizard, then force a collection on
        # this thread once it has closed — finalizing Tk where it's safe.
        gc.disable()
        try:
            # Warm the heavy libraries while the user is busy in the wizard, so
            # the post-launch wait only covers loading their data.
            _start_prewarm()
            config_path = _resolve_config_via_wizard(config_path)
        finally:
            gc.collect()  # finalize the wizard's Tk interpreter on the main thread
            gc.enable()
        if config_path is None:
            return  # user closed the wizard without choosing to launch

    if config_path is None or not config_path.exists():
        raise SystemExit(f"Config file not found: {config_path}")

    settings = _load_settings(config_path)
    apply_settings_env(settings, config_path)
    _run_window(settings)


if __name__ == "__main__":
    desktop_main()
