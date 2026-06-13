# Data Loading — Current Approach

How GUANACO loads single-cell data for visualization. Grounded in `data/loader.py`,
`data/registry.py`, and the consumption layer `utils/gene_extraction_utils.py`.

## 1. Overview

The data subsystem is three files:

- **`data/loader.py`** — the engine: source detection, the loaders, sparse/CSC
  optimization, lazy materialization helpers, and S3 genome-track loading.
- **`data/registry.py`** — reads `settings` from the active JSON config and builds
  the runtime `datasets` dict by calling `initialize_data()`.
- **`data/__init__.py`** — re-exports the public API; lazily exposes the
  `registry.datasets` / `registry.color_config` singletons.

The high-level flow:

```
registry.py (reads settings)
  -> initialize_data()                # loader.py:870
       -> DatasetBundle per dataset   # loader.py:53
            -> .adata property        # loader.py:90  (deferred read)
                 -> load_adata()      # loader.py:477 (dispatches by source + mode)
```

## 2. Two notions of "lazy"

These are independent and both in play:

1. **Bundle-level deferral** — `DatasetBundle.lazy_load` (default `True`, set from
   `settings.lazy_load` in `registry.py:36`). Nothing is read at startup; the file
   is opened only when `.adata` is first accessed (`loader.py:93`). Labels are
   computed at that point via `get_discrete_labels` (`loader.py:106`).
2. **Element-level laziness** — the cloud-native zarr path keeps `X`/`layers` as
   dask arrays on disk/cloud and reads gene columns on demand. This is what
   `backed_mode` controls (`settings.backed_mode`, default `False`).

## 3. Sources and modes

`load_adata()` (`loader.py:477`) recognizes three source kinds and two memory modes.

| Pathway | Source | Entry point | `backed` |
|---|---|---|---|
| In-memory (h5ad) | local `.h5ad` | `loader.py:569` (backed-first → `to_memory`) | `False` |
| In-memory (h5mu) | local `.h5mu` | `loader.py:558` | `False` |
| In-memory (zarr) | local/remote `.zarr` | `_load_zarr` → `loader.py:436` | `False` |
| Backed HDF5 | local `.h5ad` / `.h5mu` | `loader.py:548` | `True` |
| Cloud-native lazy (zarr) | local/remote `.zarr` | `_load_zarr_backed` (`loader.py:326`) | `True` |
| Full-read fallback | any h5ad/h5mu | `loader.py:588` | n/a |

Source detection:

- `_is_remote_uri` (`loader.py:148`) — matches `s3://`, `gs://`, `gcs://`, `az://`,
  `abfs(s)://`, `http(s)://`. Remote URIs are passed through as **strings** and never
  converted to `pathlib.Path` (which would mangle `scheme://`).
- `_source_suffix` (`loader.py:153`) — extension, tolerant of trailing slashes
  (`.zarr` stores are directories) and query strings (`...dataset.zarr?versionId=`).
- Remote sources **must be `.zarr`**; remote `.h5ad`/`.h5mu` is rejected
  (`loader.py:522-526`).

## 4. How each mode works

### In-memory h5ad (`backed=False`) — `loader.py:556-582`

The file is opened `backed="r"` **first** (metadata only; matrix stays on disk),
then only the kept cells are materialized:

```
view = ad.read_h5ad(path, backed="r")     # no full-matrix RAM spike
idx  = random subset of max_cells rows
adata = view[idx, :].to_memory()          # pull just those rows into RAM
view.file.close()                         # release the on-disk handle
_to_gene_major(adata)                      # CSR -> CSC
```

Peak RAM is bounded by the subset, so files larger than RAM still load. Result is a
fully in-memory AnnData: pandas `obs`/`var`, numpy `obsm`, CSC `X`.

### In-memory h5mu (`backed=False`) — `loader.py:557-567`

muon has no sliced `to_memory`, so it opens `read_h5mu(backed=True)`, then
`view[idx].copy()` materializes the subset, then `_to_gene_major` per modality.

### In-memory zarr (`backed=False`) — `loader.py:423-441`

Opens with `read_lazy`, slices `[idx]`, `.to_memory()`. Falls back to in-memory
`ad.read_zarr` if dask/fsspec is unavailable.

### Backed HDF5 (`backed=True`) — `loader.py:548`

```
ad.read_h5ad(path, backed="r")   # AnnData
mu.read_h5mu(path, backed=True)  # MuData
```

Returns immediately. `obs`/`var`/`obsm`/`uns` are in RAM (real pandas/numpy);
`X`/`layers` stay on disk as anndata's `BackedSparseMatrix`. **All cells are served**;
gene columns are read on demand. `adata.isbacked` is `True` and `adata.filename` is set
— both are used downstream (see §6).

### Cloud-native lazy zarr (`backed=True`) — `_load_zarr_backed`, `loader.py:326-386`

The cloud-native form of "backed":

1. `_open_zarr_group` (`loader.py:163`) opens the store from its **consolidated
   metadata** (one small fetch, not a full listing). Remote URIs are wrapped in
   `zarr.storage.FsspecStore.from_url(..., read_only=True)`.
2. `read_lazy(group)` — `X`/`layers` become dask arrays; the matrix is **never
   downloaded up front**.
3. `_backed_expression_layer` (`loader.py:295`) picks a gene-major (CSC) source to
   serve expression from: the configured `expression_layer`, else an auto-detected
   CSC layer when `X` is CSR.
4. If that source is CSC, it is re-opened one **gene per chunk** via
   `read_elem_lazy(elem, chunks=(n_obs, 1))` (`loader.py:363`) so a single-gene read
   fetches exactly one column (`read_lazy`'s default chunking batches ~1000 genes).
5. `_eager_load_annotations` (`loader.py:261`) materializes the small `obs`/`var`/
   `obsm` into pandas/numpy while `X` stays remote. `uns` is read eagerly by
   `read_lazy` into a plain dict (nested sparse matrices intact).

MuData has no lazy reader, so a MuData `.zarr` in backed mode falls back to an
in-memory `mu.read_zarr` with a warning (`loader.py:415-418`).

### Full-read fallback — `loader.py:584-600`

If backed reading raises (notably some MuData), it degrades to a full in-memory
`read_h5ad` / `read_h5mu` plus down-sample, so loading never hard-fails.

## 5. Special handling

- **AnnData vs MuData (zarr):** `_zarr_is_mudata` (`loader.py:206`) checks for a
  `mod` group / `MuData` encoding marker.
- **Sparse / gene-major:** `_to_gene_major` (`loader.py:444`) converts CSR→CSC in
  memory because the hot path is single-gene `X[:, j]` extraction — O(nnz) on CSR,
  O(nnz in that gene) on CSC. Memory-neutral (replaces, not duplicates). Backed
  datasets stay on disk and should be re-encoded as CSC offline instead.
- **Per-gene chunking (zarr):** `read_elem_lazy(..., chunks=(n_obs, 1))` so a
  single-column read does not pull a ~1000-gene chunk off remote storage.
- **Down-sampling:** `_random_row_indices` (`loader.py:128`) returns a **sorted**
  random subset (sorted for efficient backed reads); `_downsample` (`loader.py:223`)
  applies it to in-memory AnnData/MuData.
- **Lazy → concrete annotations:** `_eager_load_annotations` + `_eager_read_elem`
  (`loader.py:242-292`) read `obs`/`var`/`obsm` with one batched native
  `anndata.io.read_elem` per element (fast over remote stores), falling back to the
  lazy `Dataset2D.to_dataframe()` / dask `.compute()` conversion.
- **Discrete labels:** `get_discrete_labels` (`loader.py:606`) picks obs columns with
  `< max_unique` (50) categories, using categorical metadata (O(1)) or a batched
  `nunique`, avoiding a Python-level scan on large datasets.
- **Compat shim:** `_register_anndata_null_reader_compat` (`loader.py:20`) registers a
  reader for the `null` encoding so older/edge stores don't fail.

## 6. How loaded data is consumed

**Gene expression** is always read through `utils/gene_extraction_utils.py`:

- `extract_gene_expression` / `extract_multiple_genes` →
  `_compute_gene_vector` slices `X[:, j]` → `_densify` (`gene_extraction_utils.py:142`).
- `_densify` handles all three matrix kinds uniformly: `.compute()` on dask blocks,
  `.toarray()` on scipy-sparse, numpy passes through. This is why the same plotting
  code works on in-memory, backed-HDF5, and cloud-lazy datasets unchanged.
- Results are memoized in `GeneExpressionCache` (LRU + TTL, bounded by item count and
  bytes), keyed by `_adata_id` (backed `filename`, else `id(adata)`) and `adata.shape`
  (`gene_extraction_utils.py:39-47`). Config marker genes can be **pinned** so the
  first view is instant (`pin_genes`, used from `main.py`).

**Annotations** (`obs`/`var`/`obsm`/`uns`) are consumed **directly as pandas/numpy**
by the plots — e.g. `adata.obsm[key].shape` (`embedding_layout.py:21`),
`adata.obs[groupby].isin(...)` / `adata.obs.iloc[...]` (`violin1.py`, `heatmap.py`),
`adata.uns['paga' | 'volcano' | 'grn' | 'spatial']` (`paga.py`, `volcano.py`,
`grn_demo.py`, `embedding.py`). This is why every loader produces concrete
pandas/numpy annotations even when `X` stays lazy.

Several modules branch on `adata.isbacked` / `adata.filename` to choose
row-slice-first vs extract-all strategies and to build their own caches
(`violin1.py:65,91,107,…`; `heatmap.py:68,318`; `dotplot_callbacks.py:80`;
`heatmap_callbacks.py:85`).

## 7. Configuration

Set in the dataset JSON (see `registry.py`):

- `settings.lazy_load` (default `True`) — defer reads until first access.
- `settings.backed_mode` (default `False`) — keep `X` on disk/cloud (backed / lazy).
- `settings.max_cells` (default `10000`, `null` to disable) — random down-sample cap.
- `settings.embedding_render_backend` — `scattergl` or `datashader`.
- Per-dataset: `sc_data` (path or remote `.zarr` URI), `expression_layer` (gene-major
  CSC layer for backed reads), `markers`, default embeddings/colors.

## 8. Mode selection cheat-sheet

| Situation | Mode used | Why |
|---|---|---|
| Local `.h5ad`/`.h5mu`, fits in RAM after sampling | In-memory (backed-first) | Bounded peak RAM, CSC in memory = fast reads |
| Local `.h5ad`/`.h5mu`, serve all cells | Backed HDF5 (`backed="r"`) | Lowest steady RAM; real pandas obs + `isbacked` |
| Local `.zarr` | In-memory or cloud-native lazy | Lazy `read_lazy`; CSC per-gene chunking when backed |
| Remote `.zarr` (`s3://`/`gs://`/`https://`) | Cloud-native lazy | Consolidated metadata + on-demand column reads |
| MuData (any) | In-memory or muon-backed | No lazy reader for MuData |
