# **GUANACO: A Unified Web-Based Platform for Single-Cell Multi-Omics Data Visualization** 
<table border="0" cellspacing="0" cellpadding="0">
  <tr>
    <td width="120" align="center" valign="top">
      <img src="guanaco/assets/logo.png" width="100" alt="GUANACO logo" />
    </td>
    <td style="padding-left: 20px;">
      <strong>GUANACO</strong> (Graphical Unified Analysis and Navigation of Cellular Omics) is a Python-based platform that empowers biologists to explore multi-omics single-cell data directly in the browser with clicks<br><br>
      <strong>GUANACO</strong> leverages interactive visualizations to make data exploration and figure customization effortless.
    </td>
  </tr>
</table>

## Features


- **Various visualization plot types** – Support for:
  - Dimensionality reduction (UMAP / t-SNE)
  - Heatmaps
  - Violin plots
  - Dot plots / Matrix plots
  - Stacked plots
  - Pseudotime plots
  - Genome browser
- **Free cell selection** – Select cells with a click or lasso, define custom subpopulations as easily as drawing on paper.
- **Perception-aware tooltips** – Prevent misinterpretation by revealing actual values behind the visualization
- **100+ color maps** – Choose from a wide range of continuous and discrete palettes, including options optimized for color vision deficiencies.
- **Interactive layout** – Resize plots, reorder axes, and zoom in on details all directly in the browser

<img alt="figure" src="figure.png" />

Example Interface: [Launch the interactive demo](https://guanaco-demo.chen-sysimeta-lab.com/)


## Installation

### 1. Clone the repository
```bash
git clone https://github.com/Systems-Immunometabolism-Lab/guanaco-viz.git
cd guanaco-viz
```

### 2. (Recommended) Create a virtual environment

Using **venv** (Python ≥3.10):

```
python -m venv venv
source venv/bin/activate   # on Linux/Mac
venv\Scripts\activate      # on Windows
```

Using **conda**:

```
# Create a new conda environment (Python 3.10 or later)
conda create -n guanaco python=3.10
conda activate guanaco
```

Using **pixi**:

```
# Create and enter environment defined in pixi.toml
pixi shell
pixi add pip
```

### 3. Install from local directory

```bash
pip install .
```
Or for development (editable install):
```bash
pip install -e .
```
## Usage
```bash
guanaco -c config.json -d data_folder
```

### Command-line Options

- `-c, --config`: Name of configuration JSON file (relative to --data-dir) (default: guanaco.json)
- `-d, --data-dir`: Directory containing AnnData files referenced in config (default: current directory)
- `-p, --port`: Port to run the Dash server on (default: 4399)
- `--host`: Host to run the Dash server on (default: 0.0.0.0)
- `--debug`: Run server in debug mode
- `--max-cells`: Maximum number of cells to load per dataset (default: 8000)
- `--backed-mode`: Enable backed mode for memory-efficient loading of large datasets

## Configuration

Create a configuration JSON file specifying your datasets. See `example_config.json` for a complete example configuration. Simpliest case for visualizing scRNA data(.h5ad) is:
```
{
  "Demo": {"sc_data": "PBMC_int.h5ad"}
}
```

A complete user guide with detailed instructions, tutorials, and examples is available at: [**GUANACO User Guide**](https://systems-immunometabolism-lab.github.io/guanaco-viz/?utm_source=chatgpt.com)
