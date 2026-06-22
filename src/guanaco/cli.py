#!/usr/bin/env python3
"""Command-line interface for GUANACO visualization tool."""

import argparse
import contextlib
import io
import json
import logging
import os
import socket
import sys
import warnings
from pathlib import Path


DEFAULT_SETTINGS = {
    "host": "0.0.0.0",
    "port": 4399,
    "max_cells": 10000,
    "lazy_load": True,
    "backed_mode": False,
    "embedding_render_backend": "scattergl",
    # Public sharing via a cloudflared quick tunnel, gated by HTTP Basic Auth.
    "share": False,
    "share_username": "guanaco",
    "share_password": "",  # blank -> a random password is generated at launch
}

warnings.filterwarnings("ignore", category=FutureWarning, message=".*__version__.*deprecated.*")


def _as_bool(value, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "on"}:
            return True
        if normalized in {"false", "0", "no", "off"}:
            return False
    raise ValueError(f"settings.{field} must be a boolean")


def _load_settings(config_path: Path) -> dict:
    config = json.loads(config_path.read_text())
    settings = config.get("settings", {})
    if settings is None:
        settings = {}
    if not isinstance(settings, dict):
        raise ValueError("settings must be a JSON object")

    merged = DEFAULT_SETTINGS | settings
    merged["lazy_load"] = _as_bool(merged["lazy_load"], field="lazy_load")

    merged["backed_mode"] = _as_bool(merged["backed_mode"], field="backed_mode")

    backend = str(merged["embedding_render_backend"]).lower()
    if backend not in {"scattergl", "datashader"}:
        raise ValueError("settings.embedding_render_backend must be 'scattergl' or 'datashader'")
    merged["embedding_render_backend"] = backend

    merged["port"] = int(merged["port"])
    merged["max_cells"] = None if merged["max_cells"] is None else int(merged["max_cells"])

    merged["share"] = _as_bool(merged["share"], field="share")
    merged["share_username"] = str(merged["share_username"]).strip() or "guanaco"
    merged["share_password"] = str(merged["share_password"])
    return merged


def _dataset_names(config_path: Path) -> list[str]:
    config = json.loads(config_path.read_text())
    return [
        key for key, value in config.items()
        if key not in {"title", "color", "genome", "settings"} and isinstance(value, dict)
    ]


def _loading_message(config_path: Path) -> str:
    names = _dataset_names(config_path)
    if not names:
        return "Loading GUANACO data..."
    preview = ", ".join(names[:2])
    if len(names) > 2:
        preview += f" +{len(names) - 2} more"
    return f"Loading GUANACO data: {preview}"


def _display_urls(host: str, port: int) -> list[str]:
    if host in {"0.0.0.0", "", None}:
        urls = [f"http://127.0.0.1:{port}"]
        try:
            local_ip = socket.gethostbyname(socket.gethostname())
        except Exception:
            local_ip = None
        if local_ip and not local_ip.startswith("127."):
            urls.append(f"http://{local_ip}:{port}")
        return urls
    return [f"http://{host}:{port}"]


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="GUANACO: Interactive visualization tool for single-cell and genome browser data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=Path("guanaco.json"),
        help="Path to the configuration JSON file"
    )

    parser.add_argument(
        "--config-builder",
        action="store_true",
        help="Open a GUI to build a GUANACO config file"
    )

    args = parser.parse_args()
    
    # Validate config file exists. Relative config paths are resolved from the
    # current working directory.
    config_path = args.config.expanduser()
    if not config_path.is_absolute():
        config_path = config_path.resolve()

    if args.config_builder:
        args.config = config_path
        args.config_dir = config_path.parent
        args.settings = None
        return args

    if not config_path.exists():
        parser.error(f"Config file '{config_path}' not found.")
    args.config = config_path
    args.config_dir = config_path.parent
    try:
        args.settings = _load_settings(config_path)
    except ValueError as exc:
        parser.error(str(exc))

    return args


SPLASH_HTML = """
<!DOCTYPE html><html><head><meta charset="utf-8"><style>
  html,body{height:100%;margin:0}
  body{display:flex;flex-direction:column;align-items:center;
       justify-content:center;font-family:-apple-system,Segoe UI,Roboto,
       sans-serif;background:#0f1115;color:#e6e6e6}
  .spinner{width:46px;height:46px;border:4px solid #2a2d34;
           border-top-color:#4ea1ff;border-radius:50%;
           animation:spin 0.9s linear infinite;margin-bottom:22px}
  @keyframes spin{to{transform:rotate(360deg)}}
  .t{font-size:20px;font-weight:600;letter-spacing:.5px}
  .s{font-size:13px;color:#8a8f98;margin-top:8px}
</style></head><body>
  <div class="spinner"></div>
  <div class="t">🧬 GUANACO</div>
  <div class="s">Starting up, loading your data…</div>
</body></html>
"""


def apply_settings_env(settings: dict, config_path: Path) -> None:
    """Publish the active config + settings as environment variables.

    Data loading happens at import time inside ``guanaco.data.registry``, which
    reads these variables, so they must be set *before* ``guanaco.main`` (or
    ``guanaco.app``) is imported.
    """
    os.environ['GUANACO_CONFIG'] = str(config_path)
    if settings["max_cells"] is not None:
        os.environ['GUANACO_MAX_CELLS'] = str(settings["max_cells"])
    os.environ['GUANACO_LAZY_LOAD'] = "true" if settings["lazy_load"] else "false"
    os.environ['GUANACO_BACKED_MODE'] = str(settings["backed_mode"]).lower()
    os.environ['GUANACO_EMBEDDING_RENDER_BACKEND'] = settings["embedding_render_backend"]


def _wait_for_server(url: str, timeout: float = 30.0) -> bool:
    """Poll the local server until it responds, so the window never shows a
    blank/error page while Dash is still booting. Returns True once ready."""
    import time
    import urllib.request

    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=1):
                return True
        except Exception:
            time.sleep(0.15)
    return False


def main():
    """Main entry point for GUANACO CLI."""
    import time
    from flask import cli as flask_cli
    
    args = parse_args()

    if args.config_builder:
        from guanaco.config_wizard import launch_config_wizard

        result = launch_config_wizard(args.config)
        if result is None:
            return
        config_path, should_launch = result
        if not should_launch:
            return
        args.config = config_path
        args.config_dir = config_path.parent
        args.settings = _load_settings(config_path)
    
    # Set environment variables from CLI args (must precede importing the app)
    apply_settings_env(args.settings, args.config)

    start_time = time.time()
    
    try:
        from guanaco.utils.progress_utils import Spinner

        spinner = Spinner(_loading_message(args.config), show_result=False, stream=sys.stderr)
        spinner.start()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                import anndata
                import muon
                import matplotlib
                matplotlib.use('Agg')
                import plotly
                from guanaco.main import app
        finally:
            spinner.stop()
        
        elapsed = time.time() - start_time
        print(f"Startup completed in {elapsed:.1f} seconds")
        for url in _display_urls(args.settings["host"], args.settings["port"]):
            print(f"Open: {url}")
        
    except KeyboardInterrupt:
        print("\n\n❌ Startup cancelled by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error during startup: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Optional public sharing: gate the app behind Basic Auth, then open a tunnel.
    if args.settings["share"]:
        from guanaco import share

        public_url, username, password = share.enable_sharing(
            app.server,
            args.settings["port"],
            args.settings["share_username"],
            args.settings["share_password"],
        )
        share.print_share_banner(public_url, username, password)

    # Silence Flask/Werkzeug's development banner and per-request access logs.
    flask_cli.show_server_banner = lambda *args, **kwargs: None
    logging.getLogger("werkzeug").setLevel(logging.ERROR)

    # Now actually start the server. (dash 3 removed run_server in favour of run.)
    app.run(
        host=args.settings["host"],
        debug=False,
        port=args.settings["port"],
        dev_tools_silence_routes_logging=True,
    )

if __name__ == "__main__":
    main()
