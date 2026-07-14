Private S3 Data
===============

GUANACO can read ``.zarr`` datasets from a private AWS S3 bucket.
The server fetches data directly from S3 using your credentials;
users only ever receive rendered visualizations — the bucket URL and
raw data are never exposed to the browser.

Architecture
------------

.. code-block:: text

   User browser  ──── HTTP ────▶  GUANACO server
                                        │
                                        │  s3:// (authenticated)
                                        ▼
                                  Private S3 bucket
                                  (zarr data, not publicly accessible)

The GUANACO server holds the AWS credentials. Users connect to the
web interface and see plots, but cannot download or access the
underlying data files.

Step 1 — Make the S3 bucket private
------------------------------------

In the AWS Console:

1. Go to **S3** → select your bucket → **Permissions** tab.
2. Under **Block public access**, click **Edit**, tick all four
   checkboxes, and save.
3. Under **Bucket policy**, remove any statement that contains
   ``"Principal": "*"`` (public access). Save.

After this, the ``https://`` URL for your bucket will return 403.

Step 2 — Configure the GUANACO config file
-------------------------------------------

Replace the ``https://`` URL with an ``s3://`` URI:

.. code-block:: json

   {
     "my_dataset": {
       "sc_data": "s3://your-bucket-name/path/to/dataset.zarr",
       "description": "My private dataset"
     },
     "settings": {
       "lazy_load": true,
       "backed_mode": true
     }
   }

The URI format is ``s3://bucket-name/key`` — no region or domain in
the URL. The region is inferred from your credentials (see below).

Step 3 — Provide AWS credentials to the server
-----------------------------------------------

The server running GUANACO needs AWS credentials with read access to
the bucket. There are two ways to provide them:

**Option A — ``aws configure`` (recommended for local / single-user deployments)**

Run once in the terminal where you launch GUANACO:

.. code-block:: console

   aws configure

Enter your **Access Key ID**, **Secret Access Key**, and the bucket's
**region** (e.g. ``eu-north-1``). Credentials are saved to
``~/.aws/credentials`` and used automatically on every subsequent
``guanaco`` invocation.

**Option B — Environment variables (recommended for servers / Docker)**

Set these before starting GUANACO:

.. code-block:: console

   export AWS_ACCESS_KEY_ID=AKIA...
   export AWS_SECRET_ACCESS_KEY=...
   export AWS_DEFAULT_REGION=eu-north-1
   guanaco -c config.json

For cloud deployments (EC2, ECS, Fargate), attach an **IAM role** to
the instance/task instead. No credentials need to be set manually —
the ``s3fs`` library picks up the role automatically.

Step 4 — Launch GUANACO
------------------------

.. code-block:: console

   guanaco -c config.json

If the credentials are correct and the IAM user/role has
``s3:GetObject`` permission on the bucket, GUANACO starts normally.

Troubleshooting
---------------

``NoCredentialsError``
   Credentials were not found. Check that ``~/.aws/credentials``
   exists (``cat ~/.aws/credentials``) or that the environment
   variables are set in the same shell session.

``AccessDenied``
   Credentials are valid but the IAM user does not have read
   permission on the bucket. Add an S3 read policy to the IAM user
   or role in the AWS Console.

``EndpointResolutionError`` / wrong region
   Set ``AWS_DEFAULT_REGION`` to the region where your bucket lives
   (e.g. ``us-east-1``, ``eu-north-1``).
