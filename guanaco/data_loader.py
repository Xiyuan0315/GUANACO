"""Dataloader for single‑cell data and genome browser tracks using DatasetBundle objects."""

import json
import urllib.parse
import os
from pathlib import Path
from typing import Any, Iterable, Sequence
import anndata as ad
import boto3
import numpy as np
from botocore import UNSIGNED
from botocore.config import Config

import muon as mu
mu.set_options(pull_on_update=False)

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

    @property
    def adata(self):
        """Lazy load AnnData when first accessed."""
        if self.lazy_load and self._adata is None and self.adata_path:
            print(f"Loading dataset {self.title}...")
            max_cells = int(os.environ.get('GUANACO_MAX_CELLS', '8000'))
            seed = int(os.environ.get('GUANACO_SEED', '42')) if 'GUANACO_SEED' in os.environ else None
            backed = self.backed_mode if self.backed_mode else False
            self._adata = load_adata(self.adata_path, max_cells=max_cells, seed=seed, backed=backed)
            # Update gene markers and labels after loading
            # if self.gene_markers is None:
            #     self.gene_markers = self._adata.var_names[:6].tolist()
            if self.label_list is None:
                self.label_list = get_discrete_labels(self._adata)
        return self._adata

    def __repr__(self):
        if self._adata:
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
                adata = adata[idx, :]
                print(f"Down-sampled {path} to {max_cells} cells")
        elif isinstance(adata, mu.MuData):
            # For MuData, down-sample the primary modality
            primary_mod = list(adata.mod.keys())[0]
            n_obs = adata.mod[primary_mod].n_obs
            if n_obs > max_cells:
                rng = np.random.default_rng(seed)
                idx = rng.choice(n_obs, size=max_cells, replace=False)
                for mod in adata.mod:
                    adata.mod[mod] = adata.mod[mod][idx, :]
                print(f"Down-sampled MuData {path} to {max_cells} cells")
    elif max_cells is not None and backed:
        print(f"Warning: Downsampling not supported in backed mode for {path}")

    return adata

# ----------------------------------------------------------------------------
# Discrete label helpers
# ----------------------------------------------------------------------------

def get_discrete_labels(adata: ad.AnnData, *, max_unique: int = 50) -> list[str]:
    nunique = adata.obs.nunique()
    return nunique[nunique < max_unique].sort_values().index.tolist()

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

def _iter_or_cycle(values: Sequence[Any], n: int) -> Iterable[Any]:
    while len(values) < n:
        values = list(values) + list(values)
    return values[:n]


def load_tracks_from_s3(
    bucket_urls: Sequence[str],
    max_heights: Sequence[int | None],
    atac_names: Sequence[str],
    colors: Sequence[str] = None,
) -> dict[str, list[dict[str, Any]]]:
    if colors is None:
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # fallback default

    tracks_dict: dict[str, list[dict[str, Any]]] = {}

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

            # Use unsigned requests (for public buckets like a3s.fi)
            s3 = boto3.client(
                "s3",
                config=Config(signature_version=UNSIGNED),
                endpoint_url=endpoint_url,
            )

            kwargs = {"Bucket": bucket_name}
            if prefix:
                kwargs["Prefix"] = prefix

            contents = s3.list_objects_v2(**kwargs).get("Contents", [])
        except Exception as exc:
            print(f"Error reading S3 bucket from URL '{bucket_url}': {exc}")
            continue

        colors_iter = _iter_or_cycle(colors, len(contents))
        tracks: list[dict[str, Any]] = []

        for obj, colour in zip(contents, colors_iter):
            key = obj["Key"]
            encoded = urllib.parse.quote(key)

            if key.endswith((".bigwig", ".bw")):
                tracks.append({
                    "name": Path(key).stem,
                    "type": "wig",
                    "format": "bigwig",
                    "url": f"{bucket_url.rstrip('/')}/{encoded}",
                    "max": height,
                    "color": colour,
                })
            elif key.endswith(".bedpe"):
                tracks.append({
                    "name": Path(key).stem,
                    "type": "interaction",
                    "format": "bedpe",
                    "url": f"{bucket_url.rstrip('/')}/{encoded}",
                    "color": colour,
                })
            elif key.endswith(".bed"):
                tracks.append({
                    "name": Path(key).stem,
                    "type": "annotation",
                    "format": "bed",
                    "url": f"{bucket_url.rstrip('/')}/{encoded}",
                    "color": colour,
                })
            elif key.endswith((".bigBed", ".bb")):
                tracks.append({
                    "name": Path(key).stem,
                    "type": "annotation",
                    "format": "bigBed",
                    "url": f"{bucket_url.rstrip('/')}/{encoded}",
                    "color": colour,
                })

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
    seed: int | None = None,
    lazy_load: bool = True,
    backed_mode: bool | str = False,
) -> dict[str, DatasetBundle]:
    # Use provided paths or get from environment/defaults
    if json_path is None:
        json_path = JSON_PATH
    if base_dir is None:
        base_dir = BASE_DIR
    if max_cells is None:
        max_cells = int(os.environ.get('GUANACO_MAX_CELLS', '8000'))
    if seed is None and 'GUANACO_SEED' in os.environ:
        seed = int(os.environ['GUANACO_SEED'])
    
    cfg = load_config(json_path)
    global_colors = cfg.get("color", DEFAULT_COLORS)
    genome = cfg.get("genome", "hg38")
    datasets: dict[str, DatasetBundle] = {}

    for dataset_key, dataset_cfg in cfg.items():
        if dataset_key == "color" or dataset_key == "genome":
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
                adata = load_adata(adata_file, max_cells=max_cells, seed=seed, base_dir=base_dir, backed=backed_mode)
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
                adata_path=str(adata_path) if "sc_data" in dataset_cfg else None,
                lazy_load=lazy_load,
                backed_mode=backed_mode,
            )
            datasets[dataset_key] = dataset_bundle
        else:
            print(f"Warning: Dataset '{dataset_key}' has neither AnnData nor genome browser data. Skipping.")

    return datasets

# Initialize datasets with lazy loading by default (can be disabled via env var)
lazy_load = os.environ.get('GUANACO_LAZY_LOAD', 'true').lower() != 'false'
backed_mode = os.environ.get('GUANACO_BACKED_MODE', 'false').lower()
if backed_mode == 'true':
    backed_mode = True
elif backed_mode == 'r+':
    backed_mode = 'r+'
else:
    backed_mode = False

datasets = initialize_data(lazy_load=lazy_load, backed_mode=backed_mode)
# Don't print dataset info at startup to avoid loading data
# Use colors from config if available, otherwise use defaults
cfg = load_config(JSON_PATH)
color_config = cfg.get("color", DEFAULT_COLORS)  