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

import dash
import dash_bootstrap_components as dbc
from pathlib import Path
import os
import json

def load_config(json_path):
    if not json_path.exists():
        raise FileNotFoundError(f"Config file not found: {json_path}")
    return json.loads(json_path.read_text())

# Get config path from environment variable or use default
JSON_PATH = Path(os.environ.get("GUANACO_CONFIG", "guanaco.json"))

# Initialize the app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.LUX, "/assets/scientific_style.css"],
    suppress_callback_exceptions=True
)

config = load_config(JSON_PATH)
app.title = config.get("title", "GUANACO")
