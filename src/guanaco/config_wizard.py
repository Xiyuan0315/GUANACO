"""Tkinter config builder for creating GUANACO JSON config files."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from tkinter import BooleanVar, Canvas, StringVar, Text, Tk, Toplevel, filedialog, messagebox
from tkinter import ttk

try:
    from PIL import Image, ImageTk
except Exception: 
    Image = ImageTk = None

LOGO_PATH = Path(__file__).parent / "assets" / "configguanaco.png"


OPTIONAL_PLOTS = [
    ("Heatmap", "heatmap"),
    ("Violin", "violin"),
    ("Comparative Violin", "split-violin"),
    ("Dotplot", "dotplot"),
    ("Stacked Bar", "stacked-bar"),
    ("Expression Trend", "expression-trend"),
    ("PAGA", "paga"),
    ("Volcano", "volcano"),
    ("GRN", "grn"),
    ("Peak Browser", "peak-browser"),
]

DEFAULT_PLOTS = {"heatmap", "violin", "split-violin", "dotplot", "stacked-bar"}

DEFAULT_COLORS = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]

# Top-level keys reserved for global config; a dataset cannot use these names.
RESERVED_NAMES = {"title", "color", "genome", "settings"}

# Comfortable max width for the form content; wider windows just add margins.
CONTENT_MAX_WIDTH = 620

# The ttk Notebook draws its pane border a few px outside its content frame.
# Subtract this from the measured content inset to land on the visible pane edge.
NOTEBOOK_PANE_MARGIN = 7


class CollapsibleSection(ttk.Frame):
    def __init__(self, parent, title: str, *, open_by_default: bool = False, header_style: str = "Toggle.TLabel"):
        super().__init__(parent)
        self._title = title
        self._open = open_by_default
        # Plain clickable text with a disclosure chevron (no button chrome).
        self.header = ttk.Label(self, style=header_style, cursor="hand2")
        self.header.grid(row=0, column=0, sticky="w", pady=(12, 4))
        self.header.bind("<Button-1>", lambda _e: self._toggle())
        self.body = ttk.Frame(self)
        self.columnconfigure(0, weight=1)
        self._render_header()
        if open_by_default:
            self.body.grid(row=1, column=0, sticky="ew")

    def _render_header(self):
        arrow = "⌄" if self._open else "›"  # expanded / collapsed chevron
        self.header.configure(text=f"{arrow}  {self._title}")

    def _toggle(self):
        self._open = not self._open
        if self._open:
            self.body.grid(row=1, column=0, sticky="ew")
        else:
            self.body.grid_remove()
        self._render_header()


def _entry_row(parent, row: int, label: str, variable: StringVar, *, browse=None, example=None):
    if example:
        # Field name in normal text, with a small grey example beside it.
        label_cell = ttk.Frame(parent)
        label_cell.grid(row=row, column=0, sticky="w", pady=4)
        ttk.Label(label_cell, text=label).pack(side="left")
        ttk.Label(label_cell, text=f"  {example}", style="Hint.TLabel").pack(side="left")
    else:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
    entry = ttk.Entry(parent, textvariable=variable)
    entry.grid(row=row, column=1, sticky="ew", pady=4)
    if browse is not None:
        ttk.Button(parent, text="Browse", command=browse).grid(row=row, column=2, sticky="e", padx=(8, 0))
    return entry


def _hinted_label(parent, text: str):
    """Build a label whose trailing (parenthetical) is greyed and smaller."""
    frame = ttk.Frame(parent)
    head, sep, tail = text.partition("(")
    ttk.Label(frame, text=head.rstrip()).pack(side="left")
    if sep:
        ttk.Label(frame, text=f"  {sep}{tail}", style="Hint.TLabel").pack(side="left")
    return frame


def _text_box(parent, height: int):
    """A multi-line Text widget with a thin border matching the ttk entries."""
    return Text(
        parent,
        height=height,
        wrap="word",
        relief="flat",
        borderwidth=0,
        highlightthickness=1,
        highlightbackground="#c9c9c9",
        highlightcolor="#a3a3a3",
    )


def _split_values(raw: str) -> list[str]:
    values: list[str] = []
    for chunk in raw.replace(",", "\n").splitlines():
        value = chunk.strip()
        if value:
            values.append(value)
    return values


# Cloud/remote URI schemes accepted for sc_data (mirrors guanaco.data.loader).
_REMOTE_SCHEMES = ("s3://", "gs://", "gcs://", "az://", "abfs://", "abfss://", "http://", "https://")


def _is_remote_uri(value: str) -> bool:
    """True if value is a cloud/remote URI rather than a local filesystem path."""
    return value.lower().startswith(_REMOTE_SCHEMES)


def _zarr_suffix(value: str) -> bool:
    """True if the path/URL (ignoring trailing slash and query/fragment) ends in .zarr."""
    base = value.split("?", 1)[0].split("#", 1)[0].rstrip("/")
    return base.lower().endswith(".zarr")


def _detect_modalities(path: str) -> list[str] | None:
    """Modality names for a local ``.h5mu`` file, or ``None`` for single AnnData/unknown.

    Reads only the ``mod`` group keys from the HDF5 container (no cell data), so it
    stays cheap even for large files. Remote stores and ``.h5ad`` are treated as a
    single (unnamed) modality.
    """
    if not path:
        return None
    text = path.strip()
    if _is_remote_uri(text) or not text.lower().endswith(".h5mu"):
        return None
    try:
        import h5py

        with h5py.File(Path(text).expanduser(), "r") as f:
            if "mod" in f:
                return list(f["mod"].keys())
    except Exception:
        return None
    return None


class ViewBlock:
    """One set of default-view fields: embeddings, colors, markers, optional gene annotation.

    Used both for a single-modality dataset (``modality_name=None`` -> writes
    dataset-level ``default_*``/``markers`` keys) and for each modality of a MuData
    dataset (writes those keys plus ``gene_annotation`` under ``modalities[name]``).
    """

    def __init__(self, parent, *, modality_name: str | None = None):
        self.modality_name = modality_name
        self.embedding_left = StringVar()
        self.color_left = StringVar()
        self.embedding_right = StringVar()
        self.color_right = StringVar()
        self.gene_annotation = StringVar()

        if modality_name is None:
            self.frame = ttk.Frame(parent)
        else:
            self.frame = ttk.LabelFrame(parent, text=modality_name, padding=8)
        self.frame.columnconfigure(1, weight=1)

        row = 0
        _entry_row(self.frame, row, "Left embedding", self.embedding_left, example="e.g. X_umap"); row += 1
        _entry_row(self.frame, row, "Left color by", self.color_left, example="e.g. cell_type"); row += 1
        _entry_row(self.frame, row, "Right embedding", self.embedding_right, example="e.g. X_umap"); row += 1
        _entry_row(self.frame, row, "Right color by", self.color_right, example="e.g. CD3D"); row += 1
        if modality_name is not None:
            _entry_row(self.frame, row, "Gene annotation", self.gene_annotation, example="hg38 / URL / path"); row += 1
        _hinted_label(self.frame, "Marker genes (comma or new line)").grid(
            row=row, column=0, columnspan=3, sticky="w", pady=(6, 0)
        )
        row += 1
        self.markers_text = _text_box(self.frame, height=2)
        self.markers_text.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(2, 0))

    def to_dict(self) -> dict:
        """Config fragment for this block (empty keys omitted)."""
        d: dict = {}
        for key, var in (
            ("default_embedding_left", self.embedding_left),
            ("default_color_left", self.color_left),
            ("default_embedding_right", self.embedding_right),
            ("default_color_right", self.color_right),
        ):
            val = var.get().strip()
            if val:
                d[key] = val
        if self.modality_name is not None:
            ga = self.gene_annotation.get().strip()
            if ga:
                d["gene_annotation"] = ga
        markers = _split_values(self.markers_text.get("1.0", "end"))
        if markers:
            d["markers"] = markers
        return d


class DatasetTab:
    """Widgets and variables describing a single dataset (one Notebook tab)."""

    def __init__(self, notebook: ttk.Notebook, default_name: str):
        self.notebook = notebook
        self.frame = ttk.Frame(notebook, padding=14)
        self.frame.columnconfigure(0, weight=1)

        self.name = StringVar(value=default_name)
        self.sc_data = StringVar()
        self.description = StringVar()

        self.genome = StringVar(value="hg38")
        self.atac_names = StringVar()

        self.plot_vars = {value: BooleanVar(value=value in DEFAULT_PLOTS) for _, value in OPTIONAL_PLOTS}

        # Default-view blocks: one for single AnnData, one per modality for MuData.
        # Rebuilt whenever the selected data file's modality layout is detected.
        self._view_blocks: list[ViewBlock] = []
        self._detected_mods: list[str] | None = None
        self._last_detected_path: str | None = None
        self._has_peak_modality = False

        self._build()
        # Keep the tab label in sync with the dataset name.
        self.name.trace_add("write", lambda *_: self._sync_tab_text())

    # -- construction -------------------------------------------------------
    def _build(self):
        # Core fields ( * marks what is required).
        basics = ttk.Frame(self.frame)
        basics.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        basics.columnconfigure(1, weight=1)
        _entry_row(basics, 0, "Dataset name *", self.name)
        self.sc_data_entry = _entry_row(basics, 1, "Data file/URL *", self.sc_data, browse=self._browse_sc_data)
        _entry_row(basics, 2, "Description", self.description)
        # Reference genome is a dataset-wide property (species/assembly), but only
        # peak data (ATAC/ChIP) needs it -- so the row is hidden until something that
        # uses it is enabled (Peak Browser, IGV tracks, or a detected peak modality).
        self._genome_row = ttk.Frame(basics)
        self._genome_row.grid(row=3, column=0, columnspan=3, sticky="ew")
        self._genome_row.columnconfigure(1, weight=1)
        _entry_row(self._genome_row, 0, "Reference genome", self.genome, example="e.g. hg38 — for ATAC/peak data")
        self._genome_row.grid_remove()
        ttk.Label(basics, text="Plot tabs to show").grid(row=4, column=0, columnspan=3, sticky="w", pady=(6, 0))
        plot_box = ttk.Frame(basics)
        plot_box.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(4, 0))
        # 3 columns so longer labels (Comparative Violin, Expression Trend) fit.
        for col in range(3):
            plot_box.columnconfigure(col, weight=1, uniform="plots")
        for i, (label, value) in enumerate(OPTIONAL_PLOTS):
            ttk.Checkbutton(plot_box, text=label, variable=self.plot_vars[value]).grid(
                row=i // 3, column=i % 3, sticky="w", padx=(0, 12), pady=4
            )
        # Re-read the modality layout when the data file changes (MuData -> one
        # view block per modality; AnnData -> a single block).
        self.sc_data_entry.bind("<FocusOut>", self._detect_and_rebuild)

        # Advanced settings.
        ttk.Label(self.frame, text="Advanced settings", style="Section.TLabel").grid(
            row=1, column=0, sticky="w", pady=(14, 0)
        )

        # Default Views: reference genome + initial scatter state (per modality for
        # MuData) + marker genes + categorical palette.
        views = CollapsibleSection(self.frame, "Default Views", open_by_default=False)
        views.grid(row=2, column=0, sticky="ew")
        views.body.columnconfigure(0, weight=1)

        self._blocks_frame = ttk.Frame(views.body)
        self._blocks_frame.grid(row=0, column=0, sticky="ew")
        self._blocks_frame.columnconfigure(0, weight=1)
        self._rebuild_view_blocks(None)

        _hinted_label(views.body, "Categorical color palette (comma or new line)").grid(
            row=1, column=0, sticky="w", pady=(8, 0)
        )
        self.colors_text = _text_box(views.body, height=4)
        self.colors_text.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        self.colors_text.insert("1.0", "\n".join(DEFAULT_COLORS))

        # IGV(optional)
        igv = CollapsibleSection(self.frame, "IGV", open_by_default=False)
        igv.grid(row=3, column=0, sticky="ew")
        igv.body.columnconfigure(1, weight=1)
        ttk.Label(
            igv.body,
            text="Add bucket URLs to enable the IGV genome browser; leave empty to disable it.",
            style="Hint.TLabel",
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))
        ttk.Label(igv.body, text="Bucket URLs, one per line").grid(row=1, column=0, sticky="nw", pady=4)
        self.bucket_text = _text_box(igv.body, height=3)
        self.bucket_text.grid(row=1, column=1, sticky="ew", pady=4)
        _entry_row(igv.body, 2, "Track names", self.atac_names)

        # Reveal the Reference genome row only when the dataset needs it: the Peak
        # Browser is enabled, IGV tracks are added, or a peak modality is detected.
        self.plot_vars["peak-browser"].trace_add("write", lambda *_: self._sync_genome_visibility())
        self.bucket_text.bind("<FocusOut>", lambda _e: self._sync_genome_visibility())
        self._sync_genome_visibility()

    @staticmethod
    def _is_peak_modality(name: str) -> bool:
        low = name.lower()
        return any(tag in low for tag in ("atac", "peak", "chip", "cut"))

    def _genome_needed(self) -> bool:
        """True if the dataset uses the reference genome (Peak Browser / IGV / peak modality)."""
        return (
            self.plot_vars["peak-browser"].get()
            or self._has_peak_modality
            or bool(_split_values(self.bucket_text.get("1.0", "end")))
        )

    def _sync_genome_visibility(self, *_):
        """Show the Reference genome row iff some peak/genome feature is in use."""
        if self._genome_needed():
            self._genome_row.grid()
        else:
            self._genome_row.grid_remove()

    def _rebuild_view_blocks(self, modalities: list[str] | None):
        """(Re)build the Default-Views blocks: one per modality, or a single block."""
        for vb in self._view_blocks:
            vb.frame.destroy()
        self._view_blocks = []

        if modalities:
            for i, mod in enumerate(modalities):
                vb = ViewBlock(self._blocks_frame, modality_name=mod)
                vb.frame.grid(row=i, column=0, sticky="ew", pady=(0, 8))
                self._view_blocks.append(vb)
        else:
            vb = ViewBlock(self._blocks_frame, modality_name=None)
            vb.frame.grid(row=0, column=0, sticky="ew")
            self._view_blocks.append(vb)

    def _detect_and_rebuild(self, *_):
        """Detect the data file's modalities and rebuild the view blocks if changed."""
        path = self.sc_data.get().strip()
        if path == self._last_detected_path:
            return
        self._last_detected_path = path
        mods = _detect_modalities(path)
        if mods == self._detected_mods:
            return
        self._detected_mods = mods
        self._has_peak_modality = bool(mods) and any(self._is_peak_modality(m) for m in mods)
        self._rebuild_view_blocks(mods)
        self._sync_genome_visibility()

    def _browse_sc_data(self):
        path = filedialog.askopenfilename(
            title="Select AnnData or MuData file",
            filetypes=[
                ("AnnData / MuData", "*.h5ad *.h5mu"),
                ("AnnData", "*.h5ad"),
                ("MuData", "*.h5mu"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.sc_data.set(path)
            self._detect_and_rebuild()

    def _sync_tab_text(self):
        try:
            self.notebook.tab(self.frame, text=self.name.get().strip() or "Untitled")
        except Exception:
            pass

    # -- serialization ------------------------------------------------------
    def build_dataset(self) -> tuple[str, dict]:
        """Return (dataset_name, dataset_dict); raise ValueError if invalid."""
        name = self.name.get().strip()
        sc_data = self.sc_data.get().strip()
        if not name:
            raise ValueError("Every dataset needs a name.")
        if not sc_data:
            raise ValueError(f"Dataset '{name}': a data file or URL is required.")

        if _is_remote_uri(sc_data):
            # Remote store: keep the URL verbatim (never resolve to a local Path).
            if not _zarr_suffix(sc_data):
                raise ValueError(
                    f"Dataset '{name}': remote data must be a .zarr store URL, got: {sc_data}"
                )
            dataset: dict = {"sc_data": sc_data}
        else:
            sc_data_path = Path(sc_data).expanduser().resolve()
            if not sc_data_path.exists():
                raise ValueError(f"Dataset '{name}': file does not exist: {sc_data_path}")
            dataset = {"sc_data": str(sc_data_path)}

        description = self.description.get().strip()
        if description:
            dataset["description"] = description

        # Reference genome — saved only when the dataset actually uses it (Peak
        # Browser, IGV tracks, or a peak modality). Drives IGV's ref track and the
        # Peak Browser's gene models.
        genome_value = self.genome.get().strip()
        if genome_value and self._genome_needed():
            dataset["genome"] = genome_value

        dataset["optional_plot_components"] = [
            value for _, value in OPTIONAL_PLOTS if self.plot_vars[value].get()
        ]

        # Default views. A single block writes dataset-level keys (default_*, markers);
        # modality blocks write a `modalities` map (each with default_*, markers,
        # gene_annotation).
        if self._detected_mods:
            modalities: dict = {}
            for vb in self._view_blocks:
                block = vb.to_dict()
                if block:
                    modalities[vb.modality_name] = block
            if modalities:
                dataset["modalities"] = modalities
        elif self._view_blocks:
            dataset.update(self._view_blocks[0].to_dict())

        colors = _split_values(self.colors_text.get("1.0", "end"))
        if colors:
            dataset["color"] = colors

        # Genome browser tracks (IGV) — enabled by providing bucket URLs.
        bucket_urls = _split_values(self.bucket_text.get("1.0", "end"))
        if bucket_urls:
            dataset["bucket_urls"] = bucket_urls

            names = _split_values(self.atac_names.get())
            if names:
                if len(names) != len(bucket_urls):
                    raise ValueError(f"Dataset '{name}': track names must match the number of bucket URLs.")
                dataset["ATAC_name"] = names

        return name, dataset


class ConfigWizard:
    def __init__(self, default_config_path: Path):
        self.root = Tk()
        self.root.title("GUANACO Config Builder")
        self.root.geometry("680x800")
        self.result: tuple[Path, bool] | None = None

        style = ttk.Style()
        style.configure("Title.TLabel", font=("", 18, "bold"))
        # Top-level section headers (Dataset Settings, App Settings).
        style.configure("Header.TLabel", font=("", 15, "bold"))
        # Sub-section header 'Advanced settings'.
        style.configure("Section.TLabel", font=("", 14, "bold"))
        # Collapsible sub-section toggles — same size, unbold, so they sit below 'Advanced settings'.
        style.configure("Toggle.TLabel", font=("", 14))
        style.configure("Hint.TLabel", foreground="#5f6368", font=("", 11))
        # Highlighted primary action button (bold accent text, theme-safe).
        style.configure("Action.TButton", font=("", 13))
        style.configure("Accent.TButton", font=("", 13, "bold"), foreground="#0a7d2c")
        style.map("Accent.TButton", foreground=[("active", "#0a7d2c"), ("disabled", "#9aa0a6")])

        self._logo_img = self._load_logo()

        self.config_path = StringVar(value=str(default_config_path))
        self.title = StringVar(value="GUANACO")

        self.host = StringVar(value="0.0.0.0")
        self.port = StringVar(value="4399")
        self.max_cells = StringVar(value="10000")
        self.backed_mode = StringVar(value="true")
        self.embedding_backend = StringVar(value="scattergl")

        self.dataset_tabs: list[DatasetTab] = []

        self._build()
        self._add_dataset()

    def _load_logo(self, height: int = 86):
        """Load and scale the GUANACO logo; returns None if unavailable."""
        if ImageTk is None or not LOGO_PATH.exists():
            return None
        try:
            img = Image.open(LOGO_PATH)
            width = max(1, round(img.width * height / img.height))
            img = img.resize((width, height), Image.LANCZOS)
            return ImageTk.PhotoImage(img)
        except Exception:
            return None

    # -- construction -------------------------------------------------------
    def _build(self):
        container = ttk.Frame(self.root, padding=18)
        container.pack(fill="both", expand=True)

        scroll_canvas = Canvas(container, highlightthickness=0)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=scroll_canvas.yview)
        canvas = ttk.Frame(scroll_canvas)
        content_window = scroll_canvas.create_window((0, 0), window=canvas, anchor="nw")
        scroll_canvas.configure(yscrollcommand=scrollbar.set)
        scroll_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        def update_scroll_region(_event=None):
            bbox = scroll_canvas.bbox("all")
            if not bbox:
                return
            # Clamp the top to 0 so the view can never scroll above the content.
            scroll_canvas.configure(scrollregion=(bbox[0], 0, bbox[2], bbox[3]))
            # When everything fits, pin to the top so there's no empty scroll space.
            if (bbox[3] - bbox[1]) <= scroll_canvas.winfo_height():
                scroll_canvas.yview_moveto(0.0)

        self._update_scroll_region = update_scroll_region

        def fit_content_width(event):
            # Keep the form at a comfortable reading width and center it,
            # so fields don't stretch edge-to-edge on a wide window.
            width = min(event.width, CONTENT_MAX_WIDTH)
            x = max(0, (event.width - width) // 2)
            scroll_canvas.coords(content_window, x, 0)
            scroll_canvas.itemconfigure(content_window, width=width)

        def on_mousewheel(event):
            # Do nothing when the whole form already fits (avoids overscroll wobble).
            first, last = scroll_canvas.yview()
            if first <= 0.0 and last >= 1.0:
                return
            # Cross-platform wheel handling: Linux/X11 fires Button-4/Button-5,
            # Windows sends delta in multiples of 120, macOS sends small deltas.
            if event.num == 4:
                step = -1
            elif event.num == 5:
                step = 1
            elif sys.platform == "darwin":
                step = -int(event.delta)
            else:
                step = int(-event.delta / 120)
            if step:
                scroll_canvas.yview_scroll(step, "units")

        canvas.bind("<Configure>", update_scroll_region)
        scroll_canvas.bind("<Configure>", fit_content_width)
        scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)
        scroll_canvas.bind_all("<Button-4>", on_mousewheel)
        scroll_canvas.bind_all("<Button-5>", on_mousewheel)
        canvas.columnconfigure(0, weight=1)

        # Header: logo + title + subtitle.
        header = ttk.Frame(canvas)
        self._header = header
        header.grid(row=0, column=0, sticky="ew", pady=(0, 12))
        header.columnconfigure(1, weight=1)
        if self._logo_img is not None:
            ttk.Label(header, image=self._logo_img).grid(
                row=0, column=0, rowspan=2, sticky="n", padx=(0, 16)
            )
        ttk.Label(header, text="GUANACO Config Builder", style="Title.TLabel").grid(
            row=0, column=1, sticky="sw"
        )
        ttk.Label(
            header,
            text="Create a configuration file for one or more AnnData / MuData datasets.",
            style="Hint.TLabel",
        ).grid(row=1, column=1, sticky="nw", pady=(2, 0))

        # Dataset Settings: header + toolbar + notebook of dataset tabs.
        self._dataset_title = ttk.Label(canvas, text="Dataset Settings", style="Header.TLabel")
        self._dataset_title.grid(row=2, column=0, sticky="w", pady=(4, 6))
        datasets_frame = ttk.Frame(canvas)
        datasets_frame.grid(row=3, column=0, sticky="ew", pady=(0, 12))
        datasets_frame.columnconfigure(0, weight=1)

        self._dataset_toolbar = ttk.Frame(datasets_frame)
        self._dataset_toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Button(self._dataset_toolbar, text="+ Add Dataset", command=self._add_dataset).pack(side="left")
        ttk.Button(
            self._dataset_toolbar, text="Remove Current Dataset", command=self._remove_current_dataset
        ).pack(side="left", padx=(8, 0))

        self.notebook = ttk.Notebook(datasets_frame)
        self.notebook.grid(row=1, column=0, sticky="ew")

        # App Settings (global): output + title + runtime options, ordered by importance.
        runtime = CollapsibleSection(canvas, "App Settings", open_by_default=False, header_style="Header.TLabel")
        self._app_section = runtime
        runtime.grid(row=4, column=0, sticky="ew")
        runtime.body.columnconfigure(0, weight=1)
        # Bordered card so the expanded content reads as a contained block.
        card = ttk.LabelFrame(runtime.body, padding=12)
        card.grid(row=0, column=0, sticky="ew")
        card.columnconfigure(1, weight=1)
        _entry_row(card, 0, "Config output *", self.config_path, browse=self._browse_config)
        _entry_row(card, 1, "App title", self.title)
        ttk.Label(card, text="Backed mode").grid(row=2, column=0, sticky="w", pady=4)
        ttk.Combobox(
            card,
            textvariable=self.backed_mode,
            values=("false", "true"),
            state="readonly",
            width=14,
        ).grid(row=2, column=1, sticky="w", pady=4)
        ttk.Label(card, text="Embedding backend").grid(row=3, column=0, sticky="w", pady=4)
        ttk.Combobox(
            card,
            textvariable=self.embedding_backend,
            values=("scattergl", "datashader"),
            state="readonly",
            width=18,
        ).grid(row=3, column=1, sticky="w", pady=4)
        self.max_cells_entry = _entry_row(card, 4, "Max cells", self.max_cells, example="(blank = all)")
        self.max_cells_hint = ttk.Label(
            card, text="Ignored in backed mode (all cells served from disk)", style="Hint.TLabel"
        )
        self.max_cells_hint.grid(row=5, column=1, sticky="w")
        _entry_row(card, 6, "Port", self.port)
        _entry_row(card, 7, "Host", self.host)

        # Down-sampling is skipped in backed mode, so disable Max cells there.
        def _sync_max_cells_state(*_):
            backed = self.backed_mode.get().strip().lower() != "false"
            self.max_cells_entry.configure(state="disabled" if backed else "normal")
            if backed:
                self.max_cells_hint.grid()
            else:
                self.max_cells_hint.grid_remove()

        self.backed_mode.trace_add("write", _sync_max_cells_state)
        _sync_max_cells_state()

        actions = ttk.Frame(canvas)
        self._actions = actions
        actions.grid(row=5, column=0, sticky="ew", pady=(18, 0))
        ttk.Button(actions, text="Preview JSON", command=self._preview, style="Action.TButton").pack(side="left")
        ttk.Button(
            actions,
            text="Save And Launch GUANACO",
            command=lambda: self._save(launch=True),
            style="Accent.TButton",
        ).pack(side="right")
        ttk.Button(
            actions, text="Save Config", command=lambda: self._save(launch=False), style="Action.TButton"
        ).pack(side="right", padx=(0, 8))

        # Align all major blocks to the notebook's visible pane edge once rendered.
        self.root.after_idle(self._align_main_blocks)

    def _align_main_blocks(self):
        """The ttk Notebook insets its grey pane, so it's narrower than the
        App Settings card. Inset the title, toolbar, App Settings section, and
        action buttons to that same pane edge so all major blocks share one
        left/right boundary and both grey cards are the same width."""
        if not self.dataset_tabs:
            return
        self.root.update_idletasks()
        nb = self.notebook
        tab = self.dataset_tabs[0].frame
        if not tab.winfo_ismapped():
            self.root.after(60, self._align_main_blocks)
            return
        left = max(0, (tab.winfo_rootx() - nb.winfo_rootx()) - NOTEBOOK_PANE_MARGIN)
        right = max(
            0,
            (nb.winfo_rootx() + nb.winfo_width())
            - (tab.winfo_rootx() + tab.winfo_width())
            - NOTEBOOK_PANE_MARGIN,
        )
        if left == 0 and right == 0:
            return
        self._header.grid_configure(padx=(left, 0))
        self._dataset_title.grid_configure(padx=(left, 0))
        self._dataset_toolbar.grid_configure(padx=(left, right))
        self._app_section.grid_configure(padx=(left, right))
        self._actions.grid_configure(padx=(left, right))
        # Layout shifted; refresh the scroll region and pin to top if it now fits.
        self.root.after_idle(self._update_scroll_region)

    # -- dataset management -------------------------------------------------
    def _add_dataset(self):
        default_name = f"Dataset{len(self.dataset_tabs) + 1}"
        tab = DatasetTab(self.notebook, default_name)
        self.dataset_tabs.append(tab)
        self.notebook.add(tab.frame, text=default_name)
        self.notebook.select(tab.frame)

    def _remove_current_dataset(self):
        if len(self.dataset_tabs) <= 1:
            messagebox.showinfo("Keep One Dataset", "At least one dataset is required.", parent=self.root)
            return
        current = self.notebook.select()
        for tab in self.dataset_tabs:
            if str(tab.frame) == current:
                self.notebook.forget(tab.frame)
                self.dataset_tabs.remove(tab)
                break

    # -- serialization ------------------------------------------------------
    @staticmethod
    def _parse_backed_mode(value: str) -> bool:
        return value.strip().lower() == "true"

    def _build_config(self) -> dict:
        if not self.dataset_tabs:
            raise ValueError("At least one dataset is required.")

        datasets: dict = {}
        for tab in self.dataset_tabs:
            name, dataset = tab.build_dataset()
            if name in RESERVED_NAMES:
                raise ValueError(f"Dataset name '{name}' is reserved; choose another name.")
            if name in datasets:
                raise ValueError(f"Duplicate dataset name: '{name}'.")
            datasets[name] = dataset

        settings = {
            "host": self.host.get().strip() or "0.0.0.0",
            "port": int(self.port.get().strip()),
            "max_cells": None
            if self.max_cells.get().strip().lower() in {"", "none", "null"}
            else int(self.max_cells.get()),
            "lazy_load": True,  # always lazy-load; not user-configurable
            "backed_mode": self._parse_backed_mode(self.backed_mode.get()),
            "embedding_render_backend": self.embedding_backend.get(),
        }

        config = {
            **datasets,
            "title": self.title.get().strip() or "GUANACO",
            "settings": settings,
        }
        return config

    # -- actions ------------------------------------------------------------
    def _browse_config(self):
        path = filedialog.asksaveasfilename(
            title="Save GUANACO config",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if path:
            self.config_path.set(path)

    def _preview(self):
        try:
            config = self._build_config()
        except Exception as exc:
            messagebox.showerror("Invalid Config", str(exc), parent=self.root)
            return

        preview = Toplevel(self.root)
        preview.title("GUANACO Config Preview")
        preview.geometry("720x520")
        text = Text(preview, wrap="none")
        text.pack(fill="both", expand=True)
        text.insert("1.0", json.dumps(config, indent=2))
        text.configure(state="disabled")

    def _save(self, *, launch: bool):
        try:
            config = self._build_config()
            path = Path(self.config_path.get()).expanduser()
            if not path.is_absolute():
                path = path.resolve()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(config, indent=2) + "\n")
        except Exception as exc:
            messagebox.showerror("Could Not Save Config", str(exc), parent=self.root)
            return

        if launch:
            self.result = (path, True)
            self.root.destroy()
        else:
            messagebox.showinfo("Config Saved", f"Saved config to:\n{path}", parent=self.root)
            self.result = (path, False)

    def run(self) -> tuple[Path, bool] | None:
        self.root.mainloop()
        return self.result


def launch_config_wizard(default_config_path: Path) -> tuple[Path, bool] | None:
    wizard = ConfigWizard(default_config_path)
    return wizard.run()
