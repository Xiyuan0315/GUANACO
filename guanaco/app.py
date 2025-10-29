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
