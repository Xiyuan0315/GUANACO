# **GUANACO: A Unified Web-Based Platform for Single-Cell Multi-Omics Data Visualization** 
<table border="0" cellspacing="0" cellpadding="0">
  <tr>
    <td width="120" align="center" valign="top">
      <img src="src/guanaco/assets/logo.png" width="100" alt="GUANACO logo" />
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


## Installation (Pixi Recommended)

### 1. Clone the repository
```bash
git clone https://github.com/Systems-Immunometabolism-Lab/guanaco-viz.git
cd guanaco-viz
```

### 2. Create the Pixi environment

```bash
pixi install
```

Python compatibility: `>=3.10,<3.13` (recommended: Python 3.11 or 3.12).

### 3. Install GUANACO package into the Pixi environment

For normal use:

```bash
pixi run install
```

For development (editable install):

```
pixi run install-dev
```

### 4. Run GUANACO

```bash
pixi run guanaco -c config.json -d data_folder
```

Or use the predefined task:

```bash
pixi run run
```

### 5. Build distributable artifacts

```bash
pixi run build
```

This produces wheel and source distributions in `dist/`.

### 6. Publish to PyPI

1. Update `version` in `pyproject.toml`.
2. Build and validate package files:

```bash
pixi run build
pixi run pypi-check
```

3. Upload to TestPyPI (recommended first):

```bash
pixi run python -m twine upload --repository testpypi dist/*
```

4. Upload to PyPI:

```bash
pixi run python -m twine upload dist/*
```

### Install from PyPI (for users)

```bash
pip install guanaco-viz
guanaco --help
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
## License

This project is licensed under the GNU General Public License v3.0 (GPLv3).  
You may freely redistribute and/or modify it under the terms of the license.  

See the [LICENSE](LICENSE) file for the full text.

## Citation

If you use **GUANACO** in your research, please cite our preprint:

> **Zhang X, Kuddus M, Xia Q, Hu Y, Chen P.**  
> *GUANACO: A Unified Web-Based Platform for Single-Cell Multi-Omics Data Visualization*  
> bioRxiv 2025.09.18.677070; doi: [https://doi.org/10.1101/2025.09.18.677070](https://doi.org/10.1101/2025.09.18.677070)

*The full citation will be updated onces the peer-reviewed publication becomes available.*
