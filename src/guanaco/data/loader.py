"""Dataloader for single‑cell data and genome browser tracks using DatasetBundle objects."""

import json
import urllib.parse
import os
from pathlib import Path
from itertools import cycle
from typing import Any, Sequence
import anndata as ad
import boto3
import numpy as np
from botocore import UNSIGNED
from botocore.config import Config

import muon as mu
mu.set_options(pull_on_update=False)


def _register_anndata_null_reader_compat() -> None:
    try:
        from anndata._io.specs.registry import IOSpec, _REGISTRY
        from anndata.compat import H5Array
    except Exception:
        return

    spec = IOSpec("null", "0.1.0")
    if (H5Array, spec, frozenset()) in _REGISTRY.read:
        return

    @_REGISTRY.register_read(H5Array, spec)
    def _read_h5_null(_elem, *, _reader):
        return None


_register_anndata_null_reader_compat()

# Get paths from environment variables set by CLI, with fallbacks
BASE_DIR = Path(os.environ.get('GUANACO_DATA_DIR', '.'))
JSON_PATH = Path(os.environ.get("GUANACO_CONFIG", "guanaco.json"))

DEFAULT_COLORS: list[str] = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]

GLOBAL_CONFIG_KEYS = {"title", "color", "genome", "settings"}

class DatasetBundle:
    def __init__(
        self,
        title: str,
        description: str,
        adata: ad.AnnData | None,
        gene_markers: list[str] | None,
        label_list: list[str] | None,
        genome_tracks: dict[str, list[dict[str, Any]]] | None,
        ref_track: dict[str, str] | None,
        color_config: list[str],
        adata_path: str | None = None,
        lazy_load: bool = False,
        backed_mode: bool | str = True,
        optional_plot_components: list[str] | None = None,
        max_cells: int | None = 10_000,
    ):

        self.title = title
        self.description = description
        self._adata = adata
        self.gene_markers = gene_markers  # Only for RNA, other modalities get first 10 dynamically
        self.label_list = label_list
        self.genome_tracks = genome_tracks
        self.ref_track = ref_track
        self.color_config = color_config
        self.adata_path = adata_path
        self.lazy_load = lazy_load
        self.backed_mode = backed_mode
        self.optional_plot_components = optional_plot_components
        self.max_cells = max_cells

    @property
    def adata(self):
        """Lazy load AnnData when first accessed."""
        if self.lazy_load and self._adata is None and self.adata_path:
            print(f"Loading dataset {self.title}...")
            backed = self.backed_mode if self.backed_mode else False
            self._adata = load_adata(
                self.adata_path,
                max_cells=self.max_cells,
                backed=backed,
            )
            # Update gene markers and labels after loading
            # if self.gene_markers is None:
            #     self.gene_markers = self._adata.var_names[:6].tolist()
            if self.label_list is None:
                self.label_list = get_discrete_labels(self._adata)
        return self._adata

    def __repr__(self):
        if self._adata is not None:
            cells_info = f"{self._adata.n_obs} cells"
        elif self.lazy_load and self.adata_path:
            cells_info = "data not loaded (lazy)"
        else:
            cells_info = "no AnnData"
        tracks_info = f"{len(self.genome_tracks)} genome tracks" if self.genome_tracks else "no genome tracks"
        return f"<DatasetBundle {self.title}: {cells_info}, {tracks_info}>"

# ----------------------------------------------------------------------------
# Config helpers
# ----------------------------------------------------------------------------

def load_config(json_path: Path) -> dict[str, Any]:
    if not json_path.exists():
        raise FileNotFoundError(f"Config file not found: {json_path}")
    return json.loads(json_path.read_text())

def load_adata(
    file: str | Path,
    *,
    max_cells: int | None = 10_000,
    seed: int | None = None,
    base_dir: Path = BASE_DIR,
    backed: bool | str = False,
) -> ad.AnnData | mu.MuData:
    """
    Load a single .h5ad or .h5mu file and optionally down-sample cells.
    
    Args:
        file: Path to the data file
        max_cells: Maximum number of cells to load (downsampling if needed)
        seed: Random seed for downsampling
        base_dir: Base directory for relative paths
        backed: If True, use backed mode (disk-based). If 'r+', use read-write backed mode.

    Returns:
        AnnData (for .h5ad) or MuData (for .h5mu)
    """

    path = Path(file)
    if not path.is_absolute():
        path = base_dir / path
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    if path.suffix == ".h5mu":
        if backed:
            mode = 'r+' if backed == 'r+' else True
            adata = mu.read_h5mu(path, backed=mode)
            print(f"Loaded {path} in backed mode (disk-based)")
        else:
            adata = mu.read_h5mu(path)
    elif path.suffix == ".h5ad":
        if backed:
            # Use backed mode for disk-based access
            mode = 'r+' if backed == 'r+' else 'r'
            adata = ad.read_h5ad(path, backed=mode)
            print(f"Loaded {path} in backed mode (disk-based)")
        else:
            adata = ad.read_h5ad(path)
    else:
        raise ValueError(f"Unsupported file extension: {path.suffix}")

    if max_cells is not None and not backed:
        # Handle both AnnData and MuData separately
        # Note: Downsampling is not supported in backed mode
        if isinstance(adata, ad.AnnData):
            if adata.n_obs > max_cells:
                rng = np.random.default_rng(seed)
                idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
                idx.sort()
                # Copy to release references to the full matrix.
                adata = adata[idx, :].copy()
                print(f"Down-sampled {path} to {max_cells} cells")
        elif isinstance(adata, mu.MuData):
            # For MuData, down-sample the primary modality
            primary_mod = list(adata.mod.keys())[0]
            n_obs = adata.mod[primary_mod].n_obs
            if n_obs > max_cells:
                rng = np.random.default_rng(seed)
                idx = rng.choice(n_obs, size=max_cells, replace=False)
                idx.sort()
                for mod in adata.mod:
                    # Copy per modality to avoid retaining the full object via views.
                    adata.mod[mod] = adata.mod[mod][idx, :].copy()
                print(f"Down-sampled MuData {path} to {max_cells} cells")
    elif max_cells is not None and backed:
        print(f"Warning: Downsampling not supported in backed mode for {path}")

    return adata

# ----------------------------------------------------------------------------
# Discrete label helpers
# ----------------------------------------------------------------------------

def get_discrete_labels(adata: ad.AnnData, *, max_unique: int = 50) -> list[str]:
    obs = adata.obs
    if obs is None or obs.empty:
        return []

    # For smaller datasets, pandas vectorized nunique is usually faster.
    if adata.n_obs <= 20_000:
        nunique = obs.nunique(dropna=True)
        return nunique[nunique < max_unique].sort_values().index.tolist()

    # For large datasets, cap counting per column and stop early once threshold is reached.
    def _count_unique_with_cap(values, cap):
        seen = set()
        for v in values:
            # Skip missing values without importing pandas.
            if v is None:
                continue
            if isinstance(v, (float, np.floating)) and np.isnan(v):
                continue
            try:
                seen.add(v)
            except TypeError:
                # Handle unhashable values defensively.
                seen.add(str(v))
            if len(seen) >= cap:
                return len(seen)
        return len(seen)

    counts: list[tuple[str, int]] = []
    for col in obs.columns:
        series = obs[col]
        if str(series.dtype) == "category":
            # Integer codes for categorical columns are compact and fast to scan.
            codes = series.cat.codes.to_numpy(copy=False)
            seen_codes = set()
            for code in codes:
                if code < 0:  # missing category
                    continue
                seen_codes.add(int(code))
                if len(seen_codes) >= max_unique:
                    break
            n_unique = len(seen_codes)
        else:
            values = series.to_numpy(copy=False)
            n_unique = _count_unique_with_cap(values, max_unique)

        if n_unique < max_unique:
            counts.append((col, n_unique))

    counts.sort(key=lambda item: item[1])
    return [col for col, _ in counts]

def get_modality_variables(adata: ad.AnnData | None, modality: str = 'RNA', n_vars: int = 10) -> list[str]:
    """
    Get the first n variables for a specific modality.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    modality : str
        The modality to get variables for ('RNA', 'Protein', etc.)
    n_vars : int
        Number of variables to return (default: 10)
        
    Returns:
    --------
    list : List of variable names for the modality
    """
    if adata is None:
        return []
    
    if modality == 'RNA' or modality is None:
        # RNA modality - use var_names
        return adata.var_names[:n_vars].tolist() if hasattr(adata, 'var_names') else []
    
    elif modality == 'Protein' and 'protein' in adata.obsm:
        # Protein modality - check for feature names in uns
        if 'protein' in adata.uns:
            if 'features' in adata.uns['protein']:
                return adata.uns['protein']['features'][:n_vars].tolist()
            elif 'var_names' in adata.uns['protein']:
                return adata.uns['protein']['var_names'][:n_vars].tolist()
        # Generate generic protein names if no names found
        n_proteins = adata.obsm['protein'].shape[1]
        return [f'Protein_{i+1}' for i in range(min(n_vars, n_proteins))]
    
    elif modality in adata.layers:
        # Layer modality - use the same var_names as RNA
        return adata.var_names[:n_vars].tolist() if hasattr(adata, 'var_names') else []
    
    else:
        # Check obsm for other modalities (like ATAC)
        for key in adata.obsm.keys():
            if key.lower() == modality.lower() or (modality == 'ATAC' and key.lower() in ['atac', 'peaks', 'chromatin']):
                # Try to find feature names in uns
                if key in adata.uns and 'features' in adata.uns[key]:
                    return adata.uns[key]['features'][:n_vars].tolist()
                # Generate generic names
                n_features = adata.obsm[key].shape[1]
                return [f'{modality}_{i+1}' for i in range(min(n_vars, n_features))]
        
        # Fallback to RNA
        return adata.var_names[:n_vars].tolist() if hasattr(adata, 'var_names') else []

# ----------------------------------------------------------------------------
# S3 genome tracks
# ----------------------------------------------------------------------------


def load_tracks_from_s3(
    bucket_urls: Sequence[str],
    max_heights: Sequence[int | None],
    atac_names: Sequence[str],
    colors: Sequence[str] = None,
) -> dict[str, list[dict[str, Any]]]:
    if colors is None:
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # fallback default

    tracks_dict: dict[str, list[dict[str, Any]]] = {}
    s3_clients: dict[str | None, Any] = {}

    for bucket_url, height, atac in zip(bucket_urls, max_heights, atac_names):
        try:
            parsed = urllib.parse.urlparse(bucket_url.rstrip("/"))

            if parsed.hostname and ".s3." in parsed.hostname:
                # AWS-style: https://bucket.s3.region.amazonaws.com
                bucket_name = parsed.hostname.split(".")[0]
                endpoint_url = None  # AWS default
                prefix = parsed.path.strip("/")
            else:
                # Path-style: https://host/bucket[/prefix...]
                endpoint_url = f"{parsed.scheme}://{parsed.hostname}"
                parts = parsed.path.strip("/").split("/", 1)
                bucket_name = parts[0]
                prefix = parts[1] if len(parts) > 1 else ""

            # Reuse S3 clients per endpoint to reduce startup/network overhead.
            endpoint_key = endpoint_url or "__aws_default__"
            if endpoint_key not in s3_clients:
                s3_clients[endpoint_key] = boto3.client(
                    "s3",
                    config=Config(signature_version=UNSIGNED),
                    endpoint_url=endpoint_url,
                )
            s3 = s3_clients[endpoint_key]
        except Exception as exc:
            print(f"Error reading S3 bucket from URL '{bucket_url}': {exc}")
            continue

        tracks: list[dict[str, Any]] = []
        colors_iter = cycle(colors)
        kwargs = {"Bucket": bucket_name}
        if prefix:
            kwargs["Prefix"] = prefix
        base_url = bucket_url.rstrip("/")
        normalized_prefix = prefix.rstrip("/") if prefix else ""

        try:
            paginator = s3.get_paginator("list_objects_v2")
            pages = paginator.paginate(**kwargs)
        except Exception as exc:
            print(f"Error listing S3 objects for URL '{bucket_url}': {exc}")
            continue

        for page in pages:
            for obj in page.get("Contents", []):
                key = obj.get("Key")
                if not key:
                    continue
                colour = next(colors_iter)

                # Avoid duplicating prefix in URLs when bucket_url already includes it.
                relative_key = key
                if normalized_prefix and key.startswith(normalized_prefix + "/"):
                    relative_key = key[len(normalized_prefix) + 1 :]
                encoded = urllib.parse.quote(relative_key, safe="/")
                lower_key = key.lower()

                if lower_key.endswith((".bigwig", ".bw")):
                    tracks.append(
                        {
                            "name": Path(key).stem,
                            "type": "wig",
                            "format": "bigwig",
                            "url": f"{base_url}/{encoded}",
                            "max": height,
                            "color": colour,
                        }
                    )
                elif lower_key.endswith(".bedpe"):
                    tracks.append(
                        {
                            "name": Path(key).stem,
                            "type": "interaction",
                            "format": "bedpe",
                            "url": f"{base_url}/{encoded}",
                            "color": colour,
                        }
                    )
                elif lower_key.endswith(".bed"):
                    tracks.append(
                        {
                            "name": Path(key).stem,
                            "type": "annotation",
                            "format": "bed",
                            "url": f"{base_url}/{encoded}",
                            "color": colour,
                        }
                    )
                elif lower_key.endswith((".bigbed", ".bb")):
                    tracks.append(
                        {
                            "name": Path(key).stem,
                            "type": "annotation",
                            "format": "bigBed",
                            "url": f"{base_url}/{encoded}",
                            "color": colour,
                        }
                    )

        tracks_dict[atac] = tracks

    return tracks_dict

# Reference genomes
_REF_URLS = {
    # Human
    "hg38": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
    "hg19": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit",
    "hg18": "https://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit",

    # Mouse
    "mm39": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit",
    "mm10": "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit",
    "mm9": "https://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit",

    # Rat
    "rn6": "https://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.2bit",
    "rn5": "https://hgdownload.cse.ucsc.edu/goldenPath/rn5/bigZips/rn5.2bit",

    # Zebrafish
    "danRer11": "https://hgdownload.cse.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.2bit",
    "danRer10": "https://hgdownload.cse.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.2bit",

    # Fruit fly
    "dm6": "https://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.2bit",
    "dm3": "https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.2bit",

    # Nematode (worm)
    "ce11": "https://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit",
    "ce10": "https://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit",

    # Yeast
    "sacCer3": "https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit",

    # Chicken
    "galGal6": "https://hgdownload.cse.ucsc.edu/goldenPath/galGal6/bigZips/galGal6.2bit",

    # Xenopus
    "xenTro9": "https://hgdownload.cse.ucsc.edu/goldenPath/xenTro9/bigZips/xenTro9.2bit",

    # Dog
    "canFam3": "https://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.2bit",

    # Cow
    "bosTau9": "https://hgdownload.cse.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.2bit",

    # Pig
    "susScr11": "https://hgdownload.cse.ucsc.edu/goldenPath/susScr11/bigZips/susScr11.2bit",

    # Macaque
    "rheMac10": "https://hgdownload.cse.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.2bit"

}

def get_ref_track(genome: str) -> dict[str, str]:
    try:
        url = _REF_URLS[genome]
    except KeyError as exc:
        raise ValueError(f"Unsupported genome: {genome}") from exc
    return {"label": genome, "url": url}

# ----------------------------------------------------------------------------
# Main data loader
# ----------------------------------------------------------------------------

def initialize_data(
    json_path: Path | None = None,
    base_dir: Path | None = None,
    *,
    max_cells: int | None = None,
    lazy_load: bool = True,
    backed_mode: bool | str = False,
    cfg: dict[str, Any] | None = None,
) -> dict[str, DatasetBundle]:
    # Use provided paths or get from environment/defaults
    if json_path is None:
        json_path = JSON_PATH
    if base_dir is None:
        base_dir = json_path.parent if json_path is not None else BASE_DIR
    
    if cfg is None:
        cfg = load_config(json_path)
    global_colors = cfg.get("color", DEFAULT_COLORS)
    genome = cfg.get("genome", "hg38")
    datasets: dict[str, DatasetBundle] = {}

    for dataset_key, dataset_cfg in cfg.items():
        if dataset_key in GLOBAL_CONFIG_KEYS:
            continue
        
        # Skip if not a dataset configuration (dict)
        if not isinstance(dataset_cfg, dict):
            continue
            
        # Handle AnnData section (optional)
        adata = None
        gene_markers = None
        label_list = None
                
        if "sc_data" in dataset_cfg and dataset_cfg["sc_data"]:
            adata_file = dataset_cfg["sc_data"]
            adata_path = Path(adata_file)
            if not adata_path.is_absolute():
                adata_path = base_dir / adata_path
            
            if lazy_load:
                # Don't load data yet, just store the path
                adata = None
                gene_markers = dataset_cfg.get("markers", None)
                label_list = None
            else:
                adata = load_adata(adata_file, max_cells=max_cells, base_dir=base_dir, backed=backed_mode)
                # Use provided markers or default to first 6 genes for RNA only
                gene_markers = dataset_cfg.get("markers", None)
                label_list = get_discrete_labels(adata) if adata else None


        # Handle genome browser section (optional)
        genome_tracks = None
        ref_track = None
        
        if "bucket_urls" in dataset_cfg and dataset_cfg["bucket_urls"]:
            # Use dataset-specific genome or global genome
            dataset_genome = dataset_cfg.get("genome", genome)
            # Set defaults for optional genome browser parameters
            max_heights = dataset_cfg.get("max_height", [None] * len(dataset_cfg["bucket_urls"]))
            atac_names = dataset_cfg.get("ATAC_name", [f"Track_{i}" for i in range(len(dataset_cfg["bucket_urls"]))])
            
            genome_tracks = load_tracks_from_s3(
                dataset_cfg["bucket_urls"],
                max_heights,
                atac_names,
                global_colors,
            )
            ref_track = get_ref_track(dataset_genome)

        # Create dataset bundle only if at least one data type is present
        if (adata is not None or genome_tracks is not None or 
            (lazy_load and "sc_data" in dataset_cfg and dataset_cfg["sc_data"])):
            dataset_bundle = DatasetBundle(
                title=dataset_key,
                description=dataset_cfg.get("description", ""),
                adata=adata,
                gene_markers=gene_markers,
                label_list=label_list,
                genome_tracks=genome_tracks,
                ref_track=ref_track,
                color_config=global_colors,
                adata_path=str(adata_path) if dataset_cfg.get("sc_data") else None,
                lazy_load=lazy_load,
                backed_mode=backed_mode,
                optional_plot_components=dataset_cfg.get("optional_plot_components"),
                max_cells=max_cells,
            )
            datasets[dataset_key] = dataset_bundle
        else:
            print(f"Warning: Dataset '{dataset_key}' has neither AnnData nor genome browser data. Skipping.")

    return datasets
