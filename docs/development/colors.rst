Color System
============

GUANACO matrix plots use one color registry:
``guanaco.utils.colors``.

Palette Sources
---------------

Discrete palettes are used for categorical annotations, groups, clusters, and
composition plots. They come from three sources:

* ``utils/cvd_color.json``: custom colorblind-friendly palettes.
  These keep their configured names, such as ``okabe_ito``,
  ``tol_vibrant``, and ``krzywinski_24``.
* Plotly qualitative palettes: added automatically from ``plotly.express`` as
  ``plotly/<name>``.
* Colorcet Glasbey palettes: added automatically as ``cc/<name>`` when
  ``colorcet`` is installed.
* Crameri categorical palettes: added automatically as ``cmc/<name>`` for the
  curated ``S`` palettes such as ``cmc/batlowS`` and ``cmc/grayCS``.
* Dynamic Glasbey palettes: shown as ``glasbey/default`` and ``glasbey/safe``.
  These generate the requested number of categorical colors on demand and cache
  the result for the current Python session. ``glasbey/safe`` uses
  colorblind-safe Glasbey generation.

Continuous colormaps are used for expression values, continuous annotations,
heatmaps, dot matrices, and gene-colored PAGA views. They come from:

* Approved perceptually uniform Plotly named color scales from
  ``PLOTLY_PERCEPTUALLY_UNIFORM`` in ``pages/matrix/colors.py``. These are
  shown as ``plotly/linear_<name>`` or ``plotly/diverging_<name>`` while the
  stored value is the Plotly colorscale name.
* Colorcet linear and diverging continuous palettes, stored as ``cc:<name>``
  values and shown as ``cc/<name>``. Names come from
  ``cc.all_original_names(group="linear")`` plus
  ``cc.all_original_names(group="diverging")``. Cyclic maps and any names
  containing ``rainbow`` are not shown.
* Crameri scientific colormaps, stored as ``cmc:<name>`` values such as
  ``cmc:batlow``. They are shown as ``cmc/linear_<name>`` or
  ``cmc/diverging_<name>`` using the curated Crameri linear and diverging lists
  in ``pages/matrix/colors.py``.

The continuous dropdown label format is always ``source/type_name``:
``plotly/linear_viridis``, ``cc/linear_bgy_10_95_c74``, and
``cmc/diverging_vik``.

Dropdown options include reversed variants for every continuous colormap using
the ``_r`` suffix, for example ``plotly/linear_viridis_r`` and
``cmc/diverging_vik_r``. The overview functions show original colormaps only,
so users can inspect the canonical collection without duplicate reversed rows.

Users and developers can inspect the continuous colormap collection from Python:

.. code-block:: python

   from guanaco.utils.colors import (
       continuous_colormap_overview,
       plot_continuous_colormap_overview,
   )

   overview = continuous_colormap_overview()
   linear_options = overview["linear"]
   diverging_options = overview["diverging"]

   fig = plot_continuous_colormap_overview()
   fig = plot_continuous_colormap_overview(colormap_type="diverging")
   fig = plot_continuous_colormap_overview(source="cc")

Each overview entry contains ``source``, ``type``, ``name``, ``label``, and
``value``. The ``label`` is what appears in the dropdown; the ``value`` is what
GUANACO stores and passes to the plotting code.

Where Colors Are Used
---------------------

Layout files build dropdown options from the registry:

* ``pages/matrix/layouts/embedding_layout.py``
* ``pages/matrix/layouts/paga_layout.py``

Embedding, heatmap, and dotplot share the embedding/scatter continuous colormap
control. Heatmap and dotplot do not define separate continuous colormap
dropdowns.

Callback files receive selected palette names from Dash controls and pass
actual palettes to plot functions:

* ``pages/matrix/callbacks/scatter_callbacks.py``
* ``pages/matrix/callbacks/heatmap_callbacks.py``
* ``pages/matrix/callbacks/violin_callbacks.py``
* ``pages/matrix/callbacks/stacked_bar_callbacks.py``
* ``pages/matrix/callbacks/pseudotime_callbacks.py``
* ``pages/matrix/callbacks/paga_callbacks.py``

Plot files resolve external continuous colormap values through
``resolve_continuous_colorscale`` before passing them to Plotly:

* ``pages/matrix/plots/embedding.py``
* ``pages/matrix/plots/heatmap.py``
* ``pages/matrix/plots/dotmatrix.py``
* ``pages/matrix/plots/paga.py``

How To Add More Colors
----------------------

To add a custom discrete palette, edit ``utils/cvd_color.json``:

.. code-block:: json

   {
     "color_palettes": {
       "my_palette": ["#1B9E77", "#D95F02", "#7570B3"]
     }
   }

To add a new programmatic source, extend
``guanaco.utils.colors``. Do not add new palette-loading logic in layout,
callback, or plot files.

To add a Plotly continuous colormap, add it to
``PLOTLY_PERCEPTUALLY_UNIFORM``. Do not expose arbitrary Plotly continuous
colormaps, because some are not perceptually uniform and can distort data
interpretation.
