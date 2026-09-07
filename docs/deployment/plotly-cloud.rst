Publish from the configuration wizard
====================================

The wizard can prepare a portable dashboard and publish it using Plotly's official
Cloud CLI. Preparation is local; uploading requires a separate confirmation.
Your original configuration and datasets are not changed.

Install the publishing tools
----------------------------

For a published GUANACO release, install the optional publishing extra in the same
Python environment that runs the wizard:

.. code-block:: bash

   python -m pip install "guanaco-viz[publish]"
   guanaco --config-builder

When developing GUANACO **before its first PyPI release**, use the source checkout:

.. code-block:: bash

   python -m pip install -e ".[publish]"
   guanaco --config-builder

From a source checkout, the wizard builds the current code into a wheel and
includes that wheel in the export. You do not need to publish to PyPI first.
The build uses the active environment's build tools without build isolation;
install the project's build requirements if the build reports they are missing.
An installed release instead generates a requirement pinned to its installed
GUANACO version. The optional wheel selector overrides either default.

Prepare, review and publish
--------------------------

1. Configure your datasets and plots, then select **Publish to Plotly Cloud…**.
2. Enter an app name. Leave the team blank for your default team, or enter a team
   name from **Account / teams**. Leave **Update existing app** blank for a new app.
3. Select **Prepare export…** and choose a parent folder. A new
   ``guanaco-cloud-*`` folder is created; no existing export is overwritten.
4. Review the file list and size. Whole selected data files are copied, including
   their observations, annotations and embedded tissue images—not just visible
   plots or selected cells. Only upload data you are permitted to host.
5. Select **Sign in…** to authenticate through Plotly's browser-based login, if
   needed. Credentials are managed by Plotly's CLI, not stored in GUANACO JSON.
6. Select **Publish…** and confirm the upload. The log shows upload/build output.
   Use **Status** and **Build logs** to check deployment; a completed command alone
   does not guarantee a running app.
7. Use **Open dashboard**, **Copy link** or **Manage access / compute**. Set sharing
   permissions and suitable compute resources in Plotly Cloud. Account limits
   and charges remain subject to your Plotly plan.

The export contains ``app.py``, ``guanaco.json``, ``local-paths.json``,
``requirements.txt``, selected local resources, and a GUANACO wheel when applicable.
The entrypoint is the Python module ``app``. The cloud process resolves portable
data references to absolute paths before importing GUANACO. It does not launch
the desktop window or a temporary sharing tunnel.

Re-publish an update
-------------------

After publication, keep the export folder's ``plotly-cloud.toml``. It contains the
app ID and URL, allowing updates to target the same dashboard.

To update an app, configure the desired dashboard in the wizard, select the
previous ``plotly-cloud.toml`` in **Update existing app**, and prepare a new export.
The review and confirmation identify the app being updated. App name and team
fields apply to new apps only; updates retain the existing target and permissions.

After an upload attempt, the dialog offers status/log checks and requires a fresh
export before another upload. It selects the returned configuration for the next
preparation when available. If an upload failed or was cancelled before an app ID
was saved, check your Plotly account before retrying: an app may already exist.
Cancelling a local process cannot undo a remote upload.

Data size and access
-------------------

* The exporter enforces a conservative **200 MiB uncompressed bundle limit**,
  including the wheel. It never bypasses Plotly's standard size check or
  subsamples data automatically. For larger data, use a remote Zarr store or
  temporary sharing from your computer.
* Local Zarr stores are archived and restored to temporary storage at app startup.
  This preserves directories such as ``var`` that Plotly's default uploader would
  otherwise exclude. Allow for the additional disk space and startup time.
* Remote references are kept as references, not downloaded during export.
  Configure any required private-data credentials in Plotly Cloud, not in the
  bundle. Signed URLs, embedded URL credentials and secret configuration fields
  are rejected. Test remote access from the cloud environment.
* Local sharing passwords are removed. Plotly Cloud access control is separate
  from temporary password-protected sharing on your own computer.
* Symlinks are rejected and only referenced datasets/annotation files are copied.
  Files are hashed after preparation and checked again before upload. If you edit
  an export, prepare and review a new one before publishing from the wizard.
* Memory measurements for a particular plot are not the total memory of a hosted
  app. Check startup, loaded modalities, concurrent users and peak callback memory
  before choosing a compute tier. The wizard does not upgrade paid resources.

See the official `Plotly Cloud CLI documentation <https://dash.plotly.com/plotly-cloud/cli>`_
and `sharing documentation <https://dash.plotly.com/plotly-cloud/share>`_ for account
and hosting details. This integration targets ``plotly-cloud>=0.4.3,<0.5``.
