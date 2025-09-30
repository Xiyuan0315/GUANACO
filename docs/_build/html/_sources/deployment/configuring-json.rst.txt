Step-by-step Guidance to Run GUANACO
====================================

Command-line Usage
------------------

Basic command-line options:

.. code-block:: text

   options:
     -h, --help            Show this help message and exit
     -c CONFIG, --config CONFIG
                           Name of configuration JSON file (relative to --data-dir)
                           (default: guanaco.json)
     -d DATA_DIR, --data-dir DATA_DIR
                           Directory containing AnnData files referenced in config
                           (default: current directory)
     -p PORT, --port PORT  Port to run the Dash server on (default: 4399)
     --host HOST           Host to run the Dash server on (default: 0.0.0.0)
     --debug               Run server in debug mode (default: False)
     --max-cells MAX_CELLS Maximum number of cells to load per dataset (default: 10000)
     --backed-mode         Enable backed mode for memory-efficient loading of
                           large datasets (default: False)

.. tip::

   **Quick start:**  
   If the config file and data files are in the current directory, you can simply run:

   .. code-block:: console

      guanaco -c config.json



Config File Structure
---------------------

The GUANACO configuration file is written in JSON format and has **two main sections**:

1. **Studies**

   Each study is defined by a dictionary entry, where the *key* is the study name.  
   A study can include either **matrix-based data**, **track-based data**, or both:

   - **description** *(optional but recommended)*  

     A short text description of the study. Displayed in the interface to help distinguish between datasets.

   - **Matrix-based data**  

     - ``sc_data`` (required): Name of a ``.h5ad`` or ``.h5mu`` file.  
     - ``markers`` (optional): List of marker genes for plots (heatmaps, dot plots, violin plots, stacked bar plots, pseudotime plots).  

   - **Track-based data**  

     - ``bucket_urls`` (required): List of URLs pointing to S3 buckets containing BigWig or other supported files. See example and instruction :ref:`Additional Information 1 <additional-info-tracks>`.  
     - ``genome`` (required): Reference genome build (e.g. ``hg38``, ``mm10``). All supported genome are :ref:`Additional Information 2 <additional-info-genomes>`.  
     - ``ATAC_name`` (optional): List of custom display names for ATAC tracks, aligned with ``bucket_urls``. Defaults to ``Track1``, ``Track2`` etc.  
     - ``max_height`` (optional): List of maximum display heights for tracks, aligned with ``bucket_urls``. Defaults to autoscale.  

2. **General customization**

   Fields outside of individual studies apply globally to the entire GUANACO session:

   - ``title`` (optional): Title of the GUANACO instance.  
   - ``color`` (optional): List of custom HEX colors applied to clusters and tracks.


Example Config File
~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
     "Study1": {
       "description": "PBMC single-cell RNA and ATAC integration study",
       "sc_data": "pbmc.h5ad",
       "markers": ["MS4A1", "LYZ"],
       "genome": "hg38",
       "bucket_urls": ["https://atac-bucket-1/"],
       "ATAC_name": ["PBMC-ATAC"],
       "max_height": [20]
     },
     "Study2": {
       "description": "Mouse skin study",
       "sc_data": "skin.h5mu"
     },
     "title": "GUANACO",
     "color": ["#1f77b4", "#ff7f0e", "#2ca02c"]
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
`TACoWig_Template.ipynb <https://github.com/Systems-Immunometabolism-Lab/guanaco-viz/blob/main/TACoWig_Template.ipynb>`_

- **Interaction files**: ``.bedpe``  
- **Annotation tracks** (motifs, SNPs, etc.): ``.bed``, ``.bigBed`` (``.bb``)  

  You can generate motif tracks in BigBed format in two ways:  

  1. Using the JASPAR TFBS extraction tool:  
     https://jaspar.elixir.no/tfbs_extraction/  

  2. Using our provided script:  
     https://github.com/Systems-Immunometabolism-Lab/guanaco-viz/blob/main/motif_extraction

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


.. raw:: html
