# guanaco-viz - A Unified Web-Based Platform for Single-Cell Multi-Omics Data Visualization
# Version 1.0
#
# Copyright (C) 2025  Xiyuan Zhang
#
# This file is part of guanaco-viz.
#
# guanaco-viz is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# guanaco-viz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Project repository: https://github.com/Systems-Immunometabolism-Lab/guanaco-viz
# Contact: xiyuan315@outlook.com

#!/usr/bin/env python3

"""Command-line interface for GUANACO visualization tool."""

import argparse
import os
from pathlib import Path

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="GUANACO: Interactive visualization tool for single-cell and genome browser data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-c", "--config",
        type=str,
        default=Path("guanaco.json"),
        help="Name of configuration JSON file(relative to --data-dir)"
    )
    
    parser.add_argument(
        "-d", "--data-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory containing AnnData files referenced in config"
    )
    
    parser.add_argument(
        "-p", "--port",
        type=int,
        default=4399,
        help="Port to run the Dash server on"
    )
    
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host to run the Dash server on"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run server in debug mode"
    )
    
    parser.add_argument(
        "--max-cells",
        type=int,
        default=10000,
        help="Maximum number of cells to load per dataset"
    )
    
    parser.add_argument(
        "--backed-mode",
        action="store_true",
        help="Enable backed mode for memory-efficient loading of large datasets"
    )
    
    args = parser.parse_args()
    
    # Validate config file exists
    config_path = args.data_dir / args.config
    if not config_path.exists():
        parser.error(f"Config file '{config_path}' not found.")

    
    # Validate data directory exists
    if not args.data_dir.exists():
        parser.error(f"Data directory '{args.data_dir}' not found.")   
    return args


def main():
    """Main entry point for GUANACO CLI."""
    import time
    import sys
    from guanaco.progress_utils import Spinner
    
    args = parse_args()
    
    # Set environment variables from CLI args
    os.environ['GUANACO_DATA_DIR'] = str(args.data_dir)
    config_path = (args.data_dir / args.config).resolve()
    os.environ['GUANACO_CONFIG'] = str(config_path)
    os.environ['GUANACO_MAX_CELLS'] = str(args.max_cells)
    
    # Set backed mode if enabled
    if args.backed_mode:
        os.environ['GUANACO_BACKED_MODE'] = 'True'
    
    print("🧬 Starting GUANACO server...")
    print(f"📄 Config file: {args.config}")
    print(f"📁 Data directory: {args.data_dir}")
    print(f"🌐 Server will run on {args.host}:{args.port}")
    print(f"🐛 Debug mode: {args.debug}")
    print(f"💾 Backed mode: {args.backed_mode}")
    print("─" * 60)
    
    start_time = time.time()
    
    try:
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
        if args.debug:
            traceback.print_exc()
        sys.exit(1)
    
    # Now actually start the server
    print()  # Empty line before Flask messages
    app.run_server(host=args.host, debug=args.debug, port=args.port)

if __name__ == "__main__":
    main()
