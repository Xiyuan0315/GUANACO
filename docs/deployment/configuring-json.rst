Step-by-step Guidance to Run GUANACO
====================================

Command-line Usage
------------------

Basic command-line options:

.. code-block:: text

   options:
     -h, --help            Show this help message and exit
     -c CONFIG, --config CONFIG
                           Path to the configuration JSON file
                           (default: guanaco.json)
     --config-wizard, --generate-config
                           Open a GUI wizard to create a GUANACO config file

.. tip::

   **Quick start:**  
   Use an absolute config path when launching from a GUI or when the data is not in the current directory:

   .. code-block:: console

      guanaco -c /path/to/guanaco.json

   To create the config file with a GUI, run:

   .. code-block:: console

      guanaco --config-wizard



Config File Structure
---------------------

The GUANACO configuration file is written in JSON format and has **two main sections**:

1. **Studies**

   Each study is defined by a dictionary entry, where the *key* is the study name.  
   A study can include either **matrix-based data**, **track-based data**, or both:

   - **description** *(optional but recommended)*  

     A short text description of the study. Displayed in the interface to help distinguish between datasets.

   - **Matrix-based data**

     - ``sc_data`` (required): Path to the expression data. Supported sources:

       - A local ``.h5ad`` or ``.h5mu`` file (absolute path; relative paths are resolved against the config file directory).
       - A local ``.zarr`` AnnData/MuData store (a directory).
       - A **remote** ``.zarr`` store on cloud storage, given as a URL such as ``s3://…``, ``gs://…`` or ``https://…``. The matrix stays on the cloud and is read on demand — see :ref:`Additional Information 3 <additional-info-remote>`. Remote ``.h5ad``/``.h5mu`` is *not* supported; convert to ``.zarr`` first.
     - ``markers`` (optional): List of marker genes for plots (heatmaps, dot plots, violin plots, stacked bar plots, pseudotime plots). For remote stores these genes are pre-fetched at startup so the first plots render immediately.
     - ``expression_layer`` (optional): Name of a gene-major (CSC) layer to read expression from instead of ``X`` (e.g. ``"X_csc"``). This makes single-gene reads cheap over the network — see :ref:`Additional Information 3 <additional-info-remote>`.

   - **Track-based data**  

     - ``bucket_urls`` (required): List of URLs pointing to S3 buckets containing BigWig or other supported files. See example and instruction :ref:`Additional Information 1 <additional-info-tracks>`.  
     - ``genome`` (required): Reference genome build (e.g. ``hg38``, ``mm10``). All supported genome are :ref:`Additional Information 2 <additional-info-genomes>`.  
     - ``ATAC_name`` (optional): List of custom display names for ATAC tracks, aligned with ``bucket_urls``. Defaults to ``Track1``, ``Track2`` etc.  
     - ``max_height`` (optional): List of maximum display heights for tracks, aligned with ``bucket_urls``. Defaults to autoscale.  

2. **General customization**

   Fields outside of individual studies apply globally to the entire GUANACO session:

   - ``title`` (optional): Title of the GUANACO instance.  
   - ``color`` (optional): List of custom HEX colors applied to clusters and tracks.
   - ``settings`` (optional): Runtime settings shared by all datasets in the config:

     - ``host``: Host to run the Dash server on. Default: ``"0.0.0.0"``.
     - ``port``: Port to run the Dash server on. Default: ``4399``.
     - ``max_cells``: Maximum number of cells to load per dataset. Use ``null`` to disable downsampling. Default: ``10000``.
     - ``lazy_load``: Load AnnData only when first opened. Default: ``true``.
     - ``backed_mode``: Keep the expression matrix on disk/cloud instead of loading it into memory, for large datasets. Values: ``false``, ``true``, or ``"r+"``. Default: ``false``. Set to ``true`` for remote ``.zarr`` stores — see :ref:`Additional Information 3 <additional-info-remote>`.
     - ``embedding_render_backend``: Embedding scatter rendering backend. Values: ``"scattergl"`` or ``"datashader"``. Default: ``"scattergl"``.


Example Config File
~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
     "Study1": {
       "description": "PBMC single-cell RNA and ATAC integration study",
       "sc_data": "/Users/example/data/pbmc.h5ad",
       "markers": ["MS4A1", "LYZ"],
       "genome": "hg38",
       "bucket_urls": ["https://atac-bucket-1/"],
       "ATAC_name": ["PBMC-ATAC"],
       "max_height": [20]
     },
     "Study2": {
       "description": "Mouse skin study",
       "sc_data": "/Users/example/data/skin.h5mu"
     },
     "title": "GUANACO",
     "color": ["#1f77b4", "#ff7f0e", "#2ca02c"],
     "settings": {
       "host": "0.0.0.0",
       "port": 4399,
       "max_cells": 10000,
       "lazy_load": true,
       "backed_mode": false,
       "embedding_render_backend": "scattergl"
     }
   }

Once the config file is set and the AnnData and BigWig files are prepared, GUANACO can be accessed locally or deployed on a server.  
The data will be displayed in the designated areas of the platform, as shown in the figures below.

.. figure:: ../assets/Figure22.png
   :width: 600

.. figure:: ../assets/Figure23.png
   :width: 600



Additional Information
----------------------

.. _additional-info-tracks:

1. **Track Data Configuration (Formats, URLs, Policies, CORS)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users should upload BigWig or other supported files to a cloud storage service that provides public URLs.  
AWS S3 is recommended for simplicity and scalability.

**Supported URL formats**

- AWS-style: ``https://bucket.s3.region.amazonaws.com``  
- Path-style: ``https://host/bucket[/prefix...]``  

**Supported file formats**

- **BigWig**: ``.bigwig``, ``.bw`` (ATAC peaks)  
If you need to convert an ``.h5ad`` or ``.h5mu`` file to BigWig, refer to the Jupyter notebook example provided with **scCAMEL**:  
`TACoWig_Template.ipynb <https://github.com/Systems-Immunometabolism-Lab/guanaco-viz/blob/main/examples/notebooks/TACoWig_Template.ipynb>`_

- **Interaction files**: ``.bedpe``  
- **Annotation tracks** (motifs, SNPs, etc.): ``.bed``, ``.bigBed`` (``.bb``)  

  You can generate motif tracks in BigBed format in two ways:  

  1. Using the JASPAR TFBS extraction tool:  
     https://jaspar.elixir.no/tfbs_extraction/  

  2. Using our provided script:  
     https://github.com/Systems-Immunometabolism-Lab/guanaco-viz/blob/main/scripts/motif_extraction

**Bucket policy for public read access**

Your S3 bucket must allow public read access so that GUANACO can fetch the files.  
Use the following bucket policy as an example (replace ``name-of-your-bucket`` with your bucket name):

.. code-block:: json

   {
     "Version": "2012-10-17",
     "Statement": [
       {
         "Sid": "PublicListBucket",
         "Effect": "Allow",
         "Principal": "*",
         "Action": "s3:ListBucket",
         "Resource": "arn:aws:s3:::name-of-your-bucket"
       },
       {
         "Sid": "PublicReadGetObject",
         "Effect": "Allow",
         "Principal": "*",
         "Action": "s3:GetObject",
         "Resource": "arn:aws:s3:::name-of-your-bucket/*"
       }
     ]
   }

**CORS configuration**

To enable cross-origin requests (web browser access), configure CORS on your bucket.  
Example configuration:

.. code-block:: json

   [
     {
       "AllowedHeaders": ["*"],
       "AllowedMethods": ["GET", "HEAD"],
       "AllowedOrigins": ["*"],
       "ExposeHeaders": ["ETag"],
       "MaxAgeSeconds": 3000
     }
   ]

⚠️ **Security note:** For better security, replace ``"*"``
in ``AllowedOrigins`` with a specific domain (e.g. ``"http://example.com"``) to restrict access.



.. _additional-info-genomes:

2. **Supported genome builds for track-based data**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   GUANACO supports the following reference genomes (UCSC 2bit format):

   - Human: ``hg38``, ``hg19``, ``hg18``  
   - Mouse: ``mm39``, ``mm10``, ``mm9``  
   - Rat: ``rn6``, ``rn5``  
   - Zebrafish: ``danRer11``, ``danRer10``  
   - Fruit fly: ``dm6``, ``dm3``  
   - Worm (C. elegans): ``ce11``, ``ce10``  
   - Yeast: ``sacCer3``  
   - Chicken: ``galGal6``  
   - Xenopus: ``xenTro9``  
   - Dog: ``canFam3``  
   - Cow: ``bosTau9``  
   - Pig: ``susScr11``  
   - Macaque: ``rheMac10``  

   Each build corresponds to a UCSC-hosted 2bit file, e.g.:
   ``https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit``


.. _additional-info-remote:

3. **Cloud-backed (remote) matrix data**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GUANACO can read a single-cell matrix **directly from cloud storage without
downloading it**. This is meant for datasets that are too large to fit in memory
(or to copy onto every machine): the expression matrix stays on the remote store
and only the parts you actually look at are fetched, on demand.

**How it works**

When ``sc_data`` is a remote ``.zarr`` URL and ``backed_mode`` is ``true``, GUANACO:

1. **Opens the store from its consolidated metadata.** The remote URL is wrapped in
   an ``fsspec`` store and opened with ``anndata.experimental.read_lazy``. Opening
   costs a single small metadata fetch — GUANACO does *not* list or download the
   whole store.
2. **Keeps the expression matrix on the cloud.** ``X`` (and any ``layers``) stay as
   lazy arrays on the remote store; they are never pulled into memory up front.
3. **Loads only the small annotations into memory.** Per-cell metadata (``obs``),
   per-gene metadata (``var``) and embeddings (``obsm``, e.g. UMAP/t-SNE) are read
   once, in a single batched pass, so the app's tables, colour-by menus and scatter
   plots work normally. Startup time is therefore proportional to the size of this
   metadata, not to the size of the matrix.
4. **Fetches gene columns on demand and caches them.** When you add a gene to a
   plot, only that gene's column is read from the cloud and then cached. Adding more
   genes fetches only the new columns; genes you have already viewed are served from
   the cache. Genes listed in ``markers`` are pre-fetched at startup.

**Make per-gene reads fast: provide a gene-major (CSC) layer**

A plot reads one gene at a time, i.e. a single *column* of the matrix. On a
cell-major (CSR) matrix that touches the whole dataset; on a **gene-major (CSC)**
matrix it reads just that one column — far cheaper over the network. For responsive
remote browsing, store a CSC copy of the matrix as a layer (commonly ``X_csc``) and
point ``expression_layer`` at it:

- If ``expression_layer`` is set, GUANACO serves expression from that layer.
- Otherwise, if ``X`` is CSR and a CSC layer exists, it is detected and used
  automatically.

The chosen matrix is opened with **one gene per chunk**, so a single-gene read
fetches exactly that one column rather than a block of ~1000 genes.

**Requirements**

- The store must be a **consolidated** ``.zarr`` (written with consolidated
  metadata). A gene-major CSC layer is strongly recommended for speed.
- The matching ``fsspec`` backend must be installed for the URL scheme:
  ``s3fs`` for ``s3://``, ``gcsfs`` for ``gs://``, ``aiohttp`` for ``http(s)://``.
  These ship with GUANACO's dependencies.
- For ``https://`` access the bucket must allow public read and have **CORS**
  configured (see the policy and CORS examples under *Additional Information 1*
  above).

**Preparing a remote store**

Write your AnnData to ``.zarr`` with a CSC layer and consolidated metadata, then
upload the resulting ``.zarr`` directory to your bucket:

.. code-block:: python

   import zarr
   from scipy.sparse import csc_matrix

   adata.layers["X_csc"] = csc_matrix(adata.X)   # gene-major copy for fast column reads
   adata.write_zarr("dataset.zarr")
   zarr.consolidate_metadata("dataset.zarr")      # single-fetch opening

**Example config**

.. code-block:: json

   {
     "BALF_COVID": {
       "description": "BALF COVID atlas streamed from cloud storage",
       "sc_data": "https://vitessce-demo-data.storage.googleapis.com/anndata-demos/BALF_VIB-UGent_processed_cleaned.zarr",
       "expression_layer": "X_csc",
       "markers": ["CD3D", "CD8A"]
     },
     "title": "GUANACO cloud-backed demo",
     "settings": {
       "host": "0.0.0.0",
       "port": 4399,
       "max_cells": null,
       "lazy_load": true,
       "backed_mode": true,
       "embedding_render_backend": "scattergl"
     }
   }

.. note::

   Set ``max_cells`` to ``null`` so all cells stay available (the matrix is not held
   in memory, so there is no need to down-sample), ``lazy_load`` to ``true``, and
   ``backed_mode`` to ``true``. The first startup spends its time downloading
   ``obs``/``var``/``obsm`` and pre-fetching the ``markers``; after that, viewing a
   new gene reads just that gene's column from the cloud.


.. raw:: html
