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
import pandas as pd
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
        backed_mode: bool = True,
        optional_plot_components: list[str] | None = None,
        scatter_defaults: dict[str, str | None] | None = None,
        max_cells: int | None = 10_000,
        expression_layer: str | None = None,
        gene_annotation_path: str | None = None,
        modality_configs: dict[str, dict] | None = None,
    ):

        self.title = title
        self.description = description
        self._adata = adata
        self.gene_markers = gene_markers  # Dataset-level default; per-modality overrides in modality_configs
        self.label_list = label_list
        self.genome_tracks = genome_tracks
        self.ref_track = ref_track
        self.color_config = color_config
        self.adata_path = adata_path
        self.lazy_load = lazy_load
        self.backed_mode = backed_mode
        self.optional_plot_components = optional_plot_components  # Dataset-level default
        # Optional per-dataset scatter defaults: keys 'embedding', 'annotation', 'gene'.
        self.scatter_defaults = scatter_defaults or {}  # Dataset-level default
        self.max_cells = max_cells
        self.expression_layer = expression_layer
        self.gene_annotation_path = gene_annotation_path  # Dataset-level default
        # Per-modality overrides keyed by modality name (e.g. "rna", "atac").
        # Each value is a dict with optional keys: gene_markers, scatter_defaults,
        # optional_plot_components, gene_annotation_path.
        self.modality_configs = modality_configs or {}

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
                expression_layer=self.expression_layer,
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


def _resolve_optional_local_path(value: str | None, base_dir: Path) -> str | None:
    if not value:
        return None
    if _is_remote_uri(value):
        return str(value)
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def _resolve_gene_annotation_config(value: str | None, base_dir: Path) -> str | None:
    """Resolve the ``gene_annotation`` config value.

    A genome id (``hg38``/``GRCh38``/``mm10`` ...) or an ``http(s)``/``ftp`` URL is
    passed through untouched -- ``load_gene_annotation`` downloads those itself, and
    keeping the id intact also lets the gene track badge show e.g. ``hg38`` rather
    than a joined-up filesystem path. Only a genuine local path gets ``base_dir``
    resolution.
    """
    if not value:
        return None
    text = str(value).strip()
    from guanaco.pages.matrix.plots.gene_annotation import is_known_genome_id

    if is_known_genome_id(text) or text.lower().startswith(("http://", "https://", "ftp://")):
        return text
    return _resolve_optional_local_path(text, base_dir)

def _random_row_indices(n_obs: int, max_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Sorted random subset of row indices (sorted for efficient backed reads)."""
    return np.sort(rng.choice(n_obs, size=max_cells, replace=False))


# Cloud/remote URI schemes we accept for ``sc_data``. These must never be turned
# into ``pathlib.Path`` objects -- doing so mangles the ``scheme://`` part (e.g.
# ``s3://bucket`` collapses to ``s3:/bucket``).
_REMOTE_SCHEMES: tuple[str, ...] = (
    "s3://",
    "gs://",
    "gcs://",
    "az://",
    "abfs://",
    "abfss://",
    "http://",
    "https://",
)


def _is_remote_uri(file: str | Path) -> bool:
    """True if ``file`` is a cloud/remote URI rather than a local filesystem path."""
    return isinstance(file, str) and file.lower().startswith(_REMOTE_SCHEMES)


def _source_suffix(file: str | Path) -> str:
    """Lower-cased file extension of a local path or remote URI.

    Handles trailing slashes (``.zarr`` stores are directories) and query
    strings on remote URIs (``...dataset.zarr?versionId=...``).
    """
    text = str(file).split("?", 1)[0].split("#", 1)[0].rstrip("/")
    return os.path.splitext(text)[1].lower()


def _open_zarr_group(store: str | Path):
    """Open a zarr store's root group, preferring consolidated metadata.

    Follows the anndata remote-access tutorial: a remote URI is wrapped in an
    ``fsspec`` store and opened from its *consolidated* metadata, so opening costs a
    single small metadata fetch rather than a full listing of the store -- the key
    to fast remote "backed" access. Falls back to a plain group open when the store
    isn't consolidated, and to a direct URL open if the fsspec store can't be built.
    """
    import zarr

    if _is_remote_uri(store):
        try:
            fsstore = zarr.storage.FsspecStore.from_url(str(store), read_only=True)
            try:
                return zarr.open_consolidated(fsstore, mode="r")
            except Exception:
                return zarr.open_group(fsstore, mode="r")
        except Exception:
            return zarr.open_group(str(store), mode="r")

    try:
        return zarr.open_consolidated(str(store), mode="r")
    except Exception:
        return zarr.open_group(str(store), mode="r")


def _zarr_elem(group, path_parts):
    """Navigate to a sub-element of an opened zarr group (e.g. ('layers', 'X_csc'))."""
    elem = group
    for part in path_parts:
        elem = elem[part]
    return elem


def _zarr_encoding(group, path_parts) -> str | None:
    """Encoding-type attr of a zarr element (e.g. 'csc_matrix'), or None if unreadable."""
    try:
        return str(dict(_zarr_elem(group, path_parts).attrs).get("encoding-type", "")) or None
    except Exception:
        return None


def _group_is_mudata(grp) -> bool:
    """Detect MuData vs AnnData from an already-open zarr root group.

    Looks for the MuData ``mod`` group / encoding marker. Returns False (assume
    AnnData) if the group can't be introspected, letting the reader surface any
    real error.
    """
    try:
        if "mod" in grp:
            return True
        return str(dict(grp.attrs).get("encoding-type", "")) == "MuData"
    except Exception as exc:
        print(f"[guanaco] could not introspect zarr group ({exc}); assuming AnnData.")
        return False


def _zarr_is_mudata(store: str | Path) -> bool:
    """Best-effort MuData-vs-AnnData detection for a ``.zarr`` store (opens it).

    Prefer :func:`_group_is_mudata` when the root group is already open, so the
    store isn't opened twice (a network round-trip for a remote URI).
    """
    try:
        return _group_is_mudata(_open_zarr_group(store))
    except Exception as exc:
        print(f"[guanaco] could not introspect zarr store {store} ({exc}); assuming AnnData.")
        return False


def _downsample(data, max_cells: int | None, rng: np.random.Generator):
    """Randomly down-sample an in-memory AnnData/MuData to ``max_cells`` rows."""
    if max_cells is None:
        return data

    if hasattr(data, "mod"):  # MuData
        primary = next(iter(data.mod))
        n_obs = data.mod[primary].n_obs
        if n_obs > max_cells:
            idx = _random_row_indices(n_obs, max_cells, rng)
            return data[idx].copy()
        return data

    if data.n_obs > max_cells:
        idx = _random_row_indices(data.n_obs, max_cells, rng)
        return data[idx, :].copy()
    return data


def _eager_read_elem(group, *path):
    """Eagerly read a zarr sub-element into native pandas/numpy, or None on failure.

    ``anndata.io.read_elem`` reads an ``obs``/``var``/``obsm`` group into memory in a
    single batched pass, rather than the per-column, per-chunk fetches that the
    ``read_lazy`` -> xarray -> ``to_dataframe`` path issues. Over a high-latency
    remote store those round-trips dominate startup, so this is the fast path; any
    failure returns ``None`` so the caller can fall back to the lazy conversion.
    """
    if group is None:
        return None
    try:
        from anndata.io import read_elem

        return read_elem(_zarr_elem(group, path))
    except Exception:
        return None


def _eager_load_annotations(adata: ad.AnnData, group=None) -> None:
    """Materialize a lazy AnnData's small annotations in place, keeping X/layers lazy.

    ``read_lazy`` returns ``obs``/``var`` as xarray ``Dataset2D`` and ``obsm`` as
    dask arrays, which the app's pandas/numpy code paths don't accept. These are
    small (per-cell/per-gene metadata, embeddings), so we pull them into memory
    while leaving the large ``X``/``layers`` matrices on cloud storage for on-demand
    per-gene reads.

    When the open zarr ``group`` is supplied we read each annotation with
    :func:`_eager_read_elem` (one batched native read), which is dramatically faster
    over a remote store than materializing the lazy ``Dataset2D`` column by column.
    We fall back to the lazy ``.to_dataframe()`` / ``.compute()`` conversion whenever
    the eager read is unavailable.
    """
    if not isinstance(adata.obs, pd.DataFrame):
        df = _eager_read_elem(group, "obs")
        adata.obs = df if isinstance(df, pd.DataFrame) else adata.obs.ds.to_dataframe()
    if not isinstance(adata.var, pd.DataFrame):
        df = _eager_read_elem(group, "var")
        adata.var = df if isinstance(df, pd.DataFrame) else adata.var.ds.to_dataframe()
    for key in list(adata.obsm.keys()):
        arr = adata.obsm[key]
        if isinstance(arr, np.ndarray):
            continue
        eager = _eager_read_elem(group, "obsm", key)
        if eager is not None:
            adata.obsm[key] = eager if isinstance(eager, np.ndarray) else np.asarray(eager)
        elif hasattr(arr, "compute"):
            adata.obsm[key] = arr.compute()
        else:
            adata.obsm[key] = np.asarray(arr)


def _backed_expression_layer(group, adata, preferred: str | None = None) -> str | None:
    """Pick a gene-major layer to serve expression from for fast per-gene reads.

    Per-gene reads slice a single column (``X[:, j]``). On a cell-major
    (``csr_matrix``) matrix that scans much of the data; on a gene-major
    (``csc_matrix``) matrix it reads just that column -- far cheaper, especially
    over the network. Many stores ship a CSC copy alongside a CSR ``X`` (e.g. a
    ``layers['X_csc']``). This returns the layer to use:

    - ``preferred``: an explicit ``expression_layer`` from the config (if present).
    - otherwise auto-detect: when ``X`` is CSR and a CSC layer exists, use it.

    Returns ``None`` to keep reading from ``X``.
    """
    layer_names = list(adata.layers.keys())
    if preferred:
        if preferred in layer_names:
            return preferred
        print(f"[guanaco] expression_layer '{preferred}' not found in {layer_names}; using X.")
        return None
    try:
        if str(dict(group["X"].attrs).get("encoding-type", "")) != "csr_matrix":
            return None
        for name in layer_names:
            if str(dict(group["layers"][name].attrs).get("encoding-type", "")) == "csc_matrix":
                return name
    except Exception:
        return None
    return None


def _load_zarr_backed(store: str | Path, *, group=None, expression_layer: str | None = None) -> ad.AnnData:
    """Open an AnnData ``.zarr`` store backed: ``X`` stays remote/on-disk, metadata in memory.

    This is the cloud-native form of "backed" access described in the anndata
    ``read_lazy`` tutorial. The store is opened from its consolidated metadata and
    the expression matrix is **never downloaded up front** -- ``X``/``layers`` stay
    as lazy dask arrays on the (possibly remote ``s3://``/``gs://``/``https://``)
    store. Only the small ``obs``/``var``/``obsm`` annotations are pulled into
    memory so the app's pandas/numpy code paths work; the gene columns the user
    views are fetched on demand by the extraction layer
    (``guanaco.utils.gene_extraction_utils``), which slices ``X[:, j]`` and computes
    just that column. Adding genes fetches only the new columns; previously read
    genes come from the gene cache.

    If the store provides a gene-major (CSC) layer -- either named via the config's
    ``expression_layer`` or auto-detected when ``X`` is CSR -- that layer is swapped
    into ``X`` so the existing per-gene read path transparently hits the fast,
    column-friendly encoding.

    Crucially, the expression matrix is re-opened with one **gene per chunk**
    (``chunks=(n_obs, 1)``). ``read_lazy``'s default chunking batches ~1000 genes
    per chunk, so a single-gene ``X[:, j]`` read would otherwise pull ~1000 genes'
    worth of data off the remote store. Per-gene chunking makes a single-gene read
    fetch only that one column -- the difference between seconds and sub-second on a
    remote CSC matrix (see the anndata ``read_lazy`` tutorial).
    """
    from anndata.experimental import read_elem_lazy, read_lazy

    if group is None:
        group = _open_zarr_group(store)
    adata = read_lazy(group)

    layer = _backed_expression_layer(group, adata, expression_layer)
    expr_path = ("layers", layer) if layer is not None else ("X",)

    # Serve expression one gene per chunk so a single-gene read fetches one column.
    if _zarr_encoding(group, expr_path) == "csc_matrix":
        try:
            adata.X = read_elem_lazy(_zarr_elem(group, expr_path), chunks=(adata.n_obs, 1))
            if layer is not None:
                del adata.layers[layer]
            print(
                f"[guanaco] backed expression served from gene-major '{'/'.join(expr_path)}' "
                "rechunked to one gene per column (fast per-gene reads)"
            )
        except Exception as exc:
            print(f"[guanaco] per-gene rechunk unavailable ({exc}); using default chunking.")
            if layer is not None:
                adata.X = adata.layers[layer]
                del adata.layers[layer]
    elif layer is not None:
        # Layer isn't CSC (e.g. an explicit non-CSC choice); use it as-is.
        adata.X = adata.layers[layer]
        del adata.layers[layer]
        print(f"[guanaco] backed expression served from layer '{layer}'")

    _eager_load_annotations(adata, group)
    print(
        f"Opened {store} backed: X stays remote/lazy, "
        f"{adata.n_obs} cells x {adata.n_vars} genes, metadata in memory"
    )
    return adata


def _load_zarr(
    store: str | Path,
    *,
    max_cells: int | None,
    seed: int | None,
    backed: bool,
    expression_layer: str | None = None,
):
    """Load an AnnData/MuData ``.zarr`` store (local directory or remote URI).

    Two modes:

    - ``backed=True`` -> *backed / cloud-native*: opened with
      :func:`anndata.experimental.read_lazy`; ``X``/``layers`` stay on disk/cloud
      and only the requested gene columns are read on demand (see
      :func:`_load_zarr_backed`). MuData has no lazy reader, so it falls back to an
      in-memory read with a clear warning.
    - ``backed=False`` -> the store is opened lazily only to bound peak memory:
      the kept cells (a random ``max_cells`` subset) are materialized via
      ``.to_memory()`` and the rest never touches RAM, mirroring the ``.h5ad``
      backed-first path. Falls back to an in-memory ``read_zarr`` if ``dask`` is
      unavailable.
    """
    rng = np.random.default_rng(seed)

    # Open the store's root group once and reuse it for MuData detection, backed
    # loading and the lazy read. Opening a remote consolidated store is a network
    # round-trip, so repeating it per step doubles/triples startup latency.
    try:
        group = _open_zarr_group(store)
    except Exception as exc:
        print(f"[guanaco] could not open zarr store {store} ({exc}); deferring to in-memory read.")
        group = None

    is_mudata = _group_is_mudata(group) if group is not None else _zarr_is_mudata(store)

    if is_mudata:
        # MuData has no lazy/group reader, so it always re-reads the store in memory.
        if backed:
            print(f"[guanaco] cloud-backed lazy mode is not available for MuData .zarr ({store}); loaded in memory.")
        return _to_gene_major(_downsample(mu.read_zarr(store), max_cells, rng))

    if backed:
        return _load_zarr_backed(store, group=group, expression_layer=expression_layer)

    try:
        from anndata.experimental import read_lazy

        if group is None:
            raise RuntimeError("zarr store could not be opened")
        lazy = read_lazy(group)
    except Exception as exc:
        # dask/fsspec missing, or the store can't be opened lazily: degrade to a
        # full in-memory read so loading never hard-fails.
        print(f"[guanaco] lazy zarr read unavailable for {store} ({exc}); using in-memory read_zarr.")
        return _to_gene_major(_downsample(ad.read_zarr(store), max_cells, rng))

    n_obs = lazy.n_obs
    if max_cells is not None and n_obs > max_cells:
        idx = _random_row_indices(n_obs, max_cells, rng)
        adata = lazy[idx, :].to_memory()
        print(f"Loaded {store}: sampled {max_cells} of {n_obs} cells into memory (lazy zarr)")
    else:
        adata = lazy.to_memory()
        print(f"Loaded {store} into memory ({n_obs} cells, lazy zarr)")
    return _to_gene_major(adata)


def _to_gene_major(data):
    """Convert CSR expression matrices to CSC (gene-major) for fast per-gene reads.

    The app's hot path is single-gene extraction (``X[:, j]``). On CSR that is
    O(nnz); on CSC it is O(nnz in that gene), typically 100-1000x faster on large
    datasets, which is what makes adding/removing genes responsive. Conversion is
    memory-neutral (CSR is replaced by CSC, not duplicated) and only matters for
    in-memory mode -- backed datasets stay on disk and should be re-encoded as CSC
    offline instead. Row-subsetting (cell selection) is the cheap direction for
    CSR, but the extraction path reads whole columns and subsets the dense result,
    so CSC does not penalize the common path.
    """
    from scipy.sparse import issparse

    def _conv(m):
        return m.tocsc() if (issparse(m) and getattr(m, "format", None) == "csr") else m

    def _convert(a):
        try:
            a.X = _conv(a.X)
            for key in list(a.layers.keys()):
                a.layers[key] = _conv(a.layers[key])
        except Exception as exc:  # never let an optimization break loading
            print(f"[guanaco] CSC conversion skipped: {exc}")

    if hasattr(data, "mod"):  # MuData: convert each modality
        for mod in data.mod:
            _convert(data.mod[mod])
    else:
        _convert(data)
    return data


def _close_backed_view(view) -> None:
    """Release a backed AnnData/MuData's on-disk file handle(s); ignore errors.

    Once the cells we keep are copied into RAM, the backing file is no longer
    needed. AnnData exposes a single ``.file`` handle; MuData holds one on the
    container plus one per modality, so close them all to avoid leaking handles
    (each backed load otherwise keeps the ``.h5mu`` open for the process lifetime).
    """
    def _close(obj):
        try:
            f = getattr(obj, "file", None)
            if f is not None:
                f.close()
        except Exception:
            pass

    _close(view)
    for mod in getattr(view, "mod", {}).values():
        _close(mod)


def load_adata(
    file: str | Path,
    *,
    max_cells: int | None = 10_000,
    seed: int | None = None,
    backed: bool = False,
    expression_layer: str | None = None,
) -> ad.AnnData | mu.MuData:
    """
    Load a single .h5ad or .h5mu file, optionally down-sampling cells.

    Memory behaviour:
        - ``backed=True``  -> keep the expression matrix on disk and serve all
          cells (lowest steady-state RAM, slower per-gene reads).
        - ``backed=False`` -> the file is still opened *backed first* so the full
          matrix is never held in RAM; only the kept cells (a random subset when
          ``max_cells`` applies) are materialized into memory. This bounds peak
          memory by the subset and works for files larger than RAM. Selecting a
          subset of rows (cells) is the cheap direction for on-disk CSR data.

    Sources:
        - Local ``.h5ad`` / ``.h5mu`` files (absolute paths).
        - Local ``.zarr`` AnnData/MuData stores (absolute directory paths).
        - Cloud ``.zarr`` stores via remote URIs (``s3://``, ``gs://``,
          ``https://``, ...). Remote URIs are passed through untouched -- never
          converted to ``pathlib.Path`` -- and require the matching ``fsspec``
          backend to be installed at read time.

    Args:
        file: Absolute local path or remote URI of the data source.
        max_cells: Maximum number of cells to keep (random down-sample if larger)
        seed: Random seed for down-sampling
        backed: If True, keep the matrix on disk/cloud (backed mode) and serve all
            cells, reading gene columns on demand.
        expression_layer: For backed ``.zarr``, name of a gene-major (CSC) layer to
            serve expression from instead of ``X`` (faster per-gene reads). If None,
            a CSC layer is auto-used when ``X`` is CSR.

    Returns:
        AnnData (for .h5ad) or MuData (for .h5mu)
    """
    suffix = _source_suffix(file)

    if _is_remote_uri(file):
        # Keep the original URI string; do not Path()-mangle the scheme.
        if suffix != ".zarr":
            raise ValueError(
                f"Remote sc_data must be a .zarr store, got '{suffix or '<none>'}' for {file}. "
                "Cloud .h5ad/.h5mu reading is not supported; convert to .zarr."
            )
        return _load_zarr(
            str(file), max_cells=max_cells, seed=seed, backed=backed, expression_layer=expression_layer
        )

    path = Path(file)
    if not path.is_absolute():
        raise ValueError(f"sc_data path must be absolute: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")
    if suffix not in (".h5ad", ".h5mu", ".zarr"):
        raise ValueError(f"Unsupported file extension: {suffix}")

    if suffix == ".zarr":
        return _load_zarr(
            path, max_cells=max_cells, seed=seed, backed=backed, expression_layer=expression_layer
        )

    is_mudata = suffix == ".h5mu"

    # Backed mode requested: keep everything on disk and serve all cells.
    if backed:
        adata = mu.read_h5mu(path, backed=True) if is_mudata else ad.read_h5ad(path, backed="r")
        print(f"Loaded {path} in backed mode (disk-based)")
        return adata

    rng = np.random.default_rng(seed)

    # In-memory mode: read backed first (metadata only; matrix stays on disk),
    # then materialize just the cells we keep -- no full-matrix RAM spike.
    try:
        if is_mudata:
            view = mu.read_h5mu(path, backed=True)
            n_obs = view.n_obs
            if max_cells is not None and n_obs > max_cells:
                idx = _random_row_indices(n_obs, max_cells, rng)
                adata = view[idx].copy()
                print(f"Loaded {path}: sampled {max_cells} of {n_obs} cells (MuData) into memory")
            else:
                adata = view.copy()
                print(f"Loaded {path} (MuData) into memory ({n_obs} cells)")
            _close_backed_view(view)  # release the backed .h5mu handle(s); data is in RAM
            return _to_gene_major(adata)

        view = ad.read_h5ad(path, backed="r")
        n_obs = view.n_obs
        if max_cells is not None and n_obs > max_cells:
            idx = _random_row_indices(n_obs, max_cells, rng)
            adata = view[idx, :].to_memory()
            print(f"Loaded {path}: sampled {max_cells} of {n_obs} cells into memory")
        else:
            adata = view.to_memory()
            print(f"Loaded {path} into memory ({n_obs} cells)")
        _close_backed_view(view)  # release the on-disk handle; data is now in RAM
        return _to_gene_major(adata)

    except Exception as exc:
        # Backed reading isn't supported for every file (notably some MuData);
        # fall back to a full in-memory read + down-sample so loading never fails.
        print(f"Backed read unavailable for {path} ({exc}); using full in-memory read.")
        adata = mu.read_h5mu(path) if is_mudata else ad.read_h5ad(path)
        if max_cells is None:
            return _to_gene_major(adata)
        if is_mudata:
            primary = next(iter(adata.mod))
            if adata.mod[primary].n_obs > max_cells:
                idx = _random_row_indices(adata.mod[primary].n_obs, max_cells, rng)
                for mod in adata.mod:
                    adata.mod[mod] = adata.mod[mod][idx, :].copy()
        elif adata.n_obs > max_cells:
            idx = _random_row_indices(adata.n_obs, max_cells, rng)
            adata = adata[idx, :].copy()
        return _to_gene_major(adata)

# ----------------------------------------------------------------------------
# Discrete label helpers
# ----------------------------------------------------------------------------

def get_discrete_labels(adata: ad.AnnData, *, max_unique: int = 50) -> list[str]:
    obs = adata.obs
    if obs is None or obs.empty:
        return []

    # Cardinality per column, vectorized. Categorical columns expose their
    # cardinality as metadata (``len(cat.categories)``, O(1)); all other columns
    # use pandas' C-level ``nunique`` in one batched call. Both avoid a Python-level
    # scan of every value, which on large (100k+ cell) backed datasets cost ~10s at
    # startup -- the bottleneck this replaces.
    counts: dict[str, int] = {}
    plain_cols: list[str] = []
    for col in obs.columns:
        dtype = obs[col].dtype
        if isinstance(dtype, pd.CategoricalDtype):
            counts[col] = len(dtype.categories)
        else:
            plain_cols.append(col)

    if plain_cols:
        plain_nunique = obs[plain_cols].nunique(dropna=True)
        for col in plain_cols:
            counts[col] = int(plain_nunique[col])

    selected = [(col, n) for col, n in counts.items() if n < max_unique]
    selected.sort(key=lambda item: item[1])
    return [col for col, _ in selected]

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
    *,
    max_cells: int | None = None,
    lazy_load: bool = True,
    backed_mode: bool = False,
    cfg: dict[str, Any] | None = None,
) -> dict[str, DatasetBundle]:
    # Use provided paths or get from environment/defaults
    if json_path is None:
        json_path = JSON_PATH

    if cfg is None:
        cfg = load_config(json_path)
    global_colors = cfg.get("color", DEFAULT_COLORS)
    genome = cfg.get("genome", "hg38")
    config_base_dir = Path(json_path).expanduser().resolve().parent
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
                
        sc_data_source = None
        if "sc_data" in dataset_cfg and dataset_cfg["sc_data"]:
            adata_file = dataset_cfg["sc_data"]
            if _is_remote_uri(adata_file):
                # Cloud URI: keep the original string, skip local path validation.
                sc_data_source = str(adata_file)
            else:
                adata_path = Path(adata_file)
                if not adata_path.is_absolute():
                    raise ValueError(
                        f"Dataset '{dataset_key}': sc_data path must be absolute, got: {adata_file}"
                    )
                sc_data_source = str(adata_path)

            if lazy_load:
                # Don't load data yet, just store the path
                adata = None
                gene_markers = dataset_cfg.get("markers", None)
                label_list = None
            else:
                adata = load_adata(
                    sc_data_source,
                    max_cells=max_cells,
                    backed=backed_mode,
                    expression_layer=dataset_cfg.get("expression_layer"),
                )
                # Use provided markers or default to first 6 genes for RNA only
                gene_markers = dataset_cfg.get("markers", None)
                label_list = get_discrete_labels(adata) if adata else None


        # Per-dataset color palette; fall back to the global/default palette.
        dataset_colors = dataset_cfg.get("color", global_colors)
        gene_annotation_path = _resolve_gene_annotation_config(
            dataset_cfg.get("gene_annotation")
            or dataset_cfg.get("gene_annotation_path")
            or dataset_cfg.get("gtf_path")
            or dataset_cfg.get("gtf"),
            config_base_dir,
        )
        # Fallback: when no explicit gene_annotation is set, the genome field (which
        # is almost always a known genome id) drives both IGV and the built-in ATAC
        # browser's gene models.  Dataset-level genome trumps the global genome.
        if gene_annotation_path is None:
            dataset_genome = dataset_cfg.get("genome", genome)
            gene_annotation_path = _resolve_gene_annotation_config(dataset_genome, config_base_dir)

        # --- per-modality overrides (modalities block) -----------------------
        modality_configs: dict[str, dict] = {}
        modalities_cfg = dataset_cfg.get("modalities", {})
        if isinstance(modalities_cfg, dict):
            for mod_name, mod_cfg in modalities_cfg.items():
                if not isinstance(mod_cfg, dict):
                    continue
                # Resolve per-modality gene_annotation:
                #  1. explicit per-modality value
                #  2. dataset-level gene_annotation_path (already resolved above)
                #  3. dataset-level / global genome (already handled above via fallback)
                mod_ga = _resolve_gene_annotation_config(
                    mod_cfg.get("gene_annotation"),
                    config_base_dir,
                )
                if mod_ga is None:
                    mod_ga = gene_annotation_path  # dataset-level (may be genome fallback)

                # Per-modality markers: override the dataset-level list.
                mod_markers = mod_cfg.get("markers")
                if mod_markers is None:
                    mod_markers = gene_markers  # dataset-level fallback

                # Per-modality scatter defaults.
                mod_scatter = {
                    "embedding_left": mod_cfg.get("default_embedding_left"),
                    "embedding_right": mod_cfg.get("default_embedding_right"),
                    "color_left": mod_cfg.get("default_color_left"),
                    "color_right": mod_cfg.get("default_color_right"),
                }
                # Fallback: dataset-level value for any key not set per-modality.
                ds_scatter = {
                    "embedding_left": dataset_cfg.get("default_embedding_left"),
                    "embedding_right": dataset_cfg.get("default_embedding_right"),
                    "color_left": dataset_cfg.get("default_color_left"),
                    "color_right": dataset_cfg.get("default_color_right"),
                }
                for k in mod_scatter:
                    if mod_scatter[k] is None:
                        mod_scatter[k] = ds_scatter[k]

                mod_opts = mod_cfg.get("optional_plot_components")
                if mod_opts is None:
                    mod_opts = dataset_cfg.get("optional_plot_components")

                modality_configs[mod_name] = {
                    "gene_markers": mod_markers,
                    "scatter_defaults": mod_scatter,
                    "optional_plot_components": mod_opts,
                    "gene_annotation_path": mod_ga,
                }

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
                dataset_colors,
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
                color_config=dataset_colors,
                adata_path=sc_data_source,
                lazy_load=lazy_load,
                backed_mode=backed_mode,
                optional_plot_components=dataset_cfg.get("optional_plot_components"),
                scatter_defaults={
                    "embedding_left": dataset_cfg.get("default_embedding_left"),
                    "embedding_right": dataset_cfg.get("default_embedding_right"),
                    "color_left": dataset_cfg.get("default_color_left"),
                    "color_right": dataset_cfg.get("default_color_right"),
                },
                max_cells=max_cells,
                expression_layer=dataset_cfg.get("expression_layer"),
                gene_annotation_path=gene_annotation_path,
                modality_configs=modality_configs or None,
            )
            datasets[dataset_key] = dataset_bundle
        else:
            print(f"Warning: Dataset '{dataset_key}' has neither AnnData nor genome browser data. Skipping.")

    return datasets
