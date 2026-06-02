#!/usr/bin/env python3
"""Command-line interface for GUANACO visualization tool."""

import argparse
import json
import os
from pathlib import Path


DEFAULT_SETTINGS = {
    "host": "0.0.0.0",
    "port": 4399,
    "max_cells": 10000,
    "lazy_load": True,
    "backed_mode": False,
    "embedding_render_backend": "scattergl",
}


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

    backed_mode = merged["backed_mode"]
    if isinstance(backed_mode, str):
        normalized = backed_mode.strip().lower()
        if normalized == "r+":
            backed_mode = "r+"
        else:
            backed_mode = _as_bool(normalized, field="backed_mode")
    elif not isinstance(backed_mode, bool):
        raise ValueError("settings.backed_mode must be a boolean or 'r+'")
    merged["backed_mode"] = backed_mode

    backend = str(merged["embedding_render_backend"]).lower()
    if backend not in {"scattergl", "datashader"}:
        raise ValueError("settings.embedding_render_backend must be 'scattergl' or 'datashader'")
    merged["embedding_render_backend"] = backend

    merged["port"] = int(merged["port"])
    merged["max_cells"] = None if merged["max_cells"] is None else int(merged["max_cells"])
    return merged


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
        "--config-wizard",
        "--generate-config",
        action="store_true",
        help="Open a GUI wizard to create a GUANACO config file"
    )
    
    args = parser.parse_args()
    
    # Validate config file exists. Relative config paths are resolved from the
    # current working directory.
    config_path = args.config.expanduser()
    if not config_path.is_absolute():
        config_path = config_path.resolve()

    if args.config_wizard:
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


def main():
    """Main entry point for GUANACO CLI."""
    import time
    import sys
    
    args = parse_args()

    if args.config_wizard:
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
    
    # Set environment variables from CLI args
    os.environ['GUANACO_DATA_DIR'] = str(args.config_dir)
    os.environ['GUANACO_CONFIG'] = str(args.config)
    if args.settings["max_cells"] is not None:
        os.environ['GUANACO_MAX_CELLS'] = str(args.settings["max_cells"])
    os.environ['GUANACO_LAZY_LOAD'] = "true" if args.settings["lazy_load"] else "false"
    os.environ['GUANACO_BACKED_MODE'] = str(args.settings["backed_mode"]).lower()
    os.environ['GUANACO_EMBEDDING_RENDER_BACKEND'] = args.settings["embedding_render_backend"]
    
    print("🧬 Starting GUANACO server...")
    print(f"📄 Config file: {args.config}")
    print(f"📁 Relative data path base: {args.config_dir}")
    print(f"🌐 Server will run on {args.settings['host']}:{args.settings['port']}")
    print(f"💾 Backed mode: {args.settings['backed_mode']}")
    print(f"📊 Embedding render backend: {args.settings['embedding_render_backend']}")
    print("─" * 60)
    
    start_time = time.time()
    
    try:
        from guanaco.utils.progress_utils import Spinner

        # Load core libraries
        with Spinner("Loading core libraries (scanpy, muon)...Please wait 5-10 seconds..."):
            import anndata 
            import muon 
        
        # Load visualization
        with Spinner("Loading visualization libraries (matplotlib, plotly)..."):
            import matplotlib
            matplotlib.use('Agg')
            import plotly
        
        
        with Spinner("Importing GUANACO modules and Loading data..."):
            from guanaco.main import app
        
        elapsed = time.time() - start_time
        print(f"\n✅ Startup completed in {elapsed:.1f} seconds")
        print("─" * 60)
        
    except KeyboardInterrupt:
        print("\n\n❌ Startup cancelled by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error during startup: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Now actually start the server
    print()  # Empty line before Flask messages
    app.run_server(
        host=args.settings["host"],
        debug=False,
        port=args.settings["port"],
    )

if __name__ == "__main__":
    main()
