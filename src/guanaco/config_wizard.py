"""Tkinter config wizard for creating GUANACO JSON config files."""

from __future__ import annotations

import json
from pathlib import Path
from tkinter import BooleanVar, Canvas, StringVar, Text, Tk, Toplevel, filedialog, messagebox
from tkinter import ttk


OPTIONAL_PLOTS = [
    ("Heatmap", "heatmap"),
    ("Violin", "violin"),
    ("Dotplot", "dotplot"),
    ("Stacked Bar", "stacked-bar"),
    ("Expression Trend", "expression-trend"),
    ("PAGA", "paga"),
    ("Volcano", "volcano"),
    ("GRN", "grn"),
]

DEFAULT_COLORS = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]


class CollapsibleSection(ttk.Frame):
    def __init__(self, parent, title: str, *, open_by_default: bool = False):
        super().__init__(parent)
        self._open = BooleanVar(value=open_by_default)
        self.header = ttk.Checkbutton(
            self,
            text=title,
            variable=self._open,
            command=self._toggle,
            style="Toolbutton",
        )
        self.header.grid(row=0, column=0, sticky="ew", pady=(12, 4))
        self.body = ttk.Frame(self)
        self.columnconfigure(0, weight=1)
        if open_by_default:
            self.body.grid(row=1, column=0, sticky="ew")

    def _toggle(self):
        if self._open.get():
            self.body.grid(row=1, column=0, sticky="ew")
        else:
            self.body.grid_remove()


class ConfigWizard:
    def __init__(self, default_config_path: Path):
        self.root = Tk()
        self.root.title("GUANACO Config Wizard")
        self.root.geometry("880x760")
        self.result: tuple[Path, bool] | None = None

        style = ttk.Style()
        style.configure("Title.TLabel", font=("", 18, "bold"))
        style.configure("Section.TLabel", font=("", 12, "bold"))
        style.configure("Hint.TLabel", foreground="#5f6368")

        self.config_path = StringVar(value=str(default_config_path))
        self.title = StringVar(value="GUANACO")
        self.dataset_name = StringVar(value="Dataset1")
        self.sc_data = StringVar()
        self.description = StringVar()

        self.host = StringVar(value="0.0.0.0")
        self.port = StringVar(value="4399")
        self.max_cells = StringVar(value="10000")
        self.lazy_load = BooleanVar(value=True)
        self.backed_mode = StringVar(value="false")
        self.embedding_backend = StringVar(value="scattergl")

        self.enable_igv = BooleanVar(value=False)
        self.genome = StringVar(value="hg38")
        self.atac_names = StringVar()
        self.max_heights = StringVar()

        self.plot_vars = {
            value: BooleanVar(value=value in {"heatmap", "violin", "dotplot", "stacked-bar"})
            for _, value in OPTIONAL_PLOTS
        }

        self._build()

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
            scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

        def fit_content_width(event):
            scroll_canvas.itemconfigure(content_window, width=event.width)

        def on_mousewheel(event):
            scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind("<Configure>", update_scroll_region)
        scroll_canvas.bind("<Configure>", fit_content_width)
        scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)

        ttk.Label(canvas, text="GUANACO Config Wizard", style="Title.TLabel").grid(row=0, column=0, sticky="w")

        row = 2
        important = ttk.LabelFrame(canvas, text="Project Overview", padding=14)
        important.grid(row=row, column=0, sticky="ew", pady=(0, 12))
        important.columnconfigure(1, weight=1)
        ttk.Label(
            important,
            text="Only the AnnData / MuData file is required. Other fields can keep their defaults.",
            style="Hint.TLabel",
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._entry_row(important, 1, "AnnData / MuData file *", self.sc_data, browse=self._browse_sc_data)
        self._entry_row(important, 2, "Config output", self.config_path, browse=self._browse_config)
        self._entry_row(important, 3, "Dataset name", self.dataset_name)
        self._entry_row(important, 4, "App title", self.title)
        self._entry_row(important, 5, "Dataset description", self.description)
        row += 1

        marker_frame = ttk.LabelFrame(canvas, text="Default Genes And Plot Tabs", padding=14)
        marker_frame.grid(row=row, column=0, sticky="ew", pady=(0, 12))
        marker_frame.columnconfigure(0, weight=1)
        ttk.Label(marker_frame, text="Default gene list, separated by comma or new line").grid(row=0, column=0, sticky="w")
        self.markers_text = Text(marker_frame, height=4, wrap="word")
        self.markers_text.grid(row=1, column=0, sticky="ew", pady=(4, 10))

        plot_box = ttk.Frame(marker_frame)
        plot_box.grid(row=2, column=0, sticky="ew")
        for i, (label, value) in enumerate(OPTIONAL_PLOTS):
            ttk.Checkbutton(plot_box, text=label, variable=self.plot_vars[value]).grid(
                row=i // 4,
                column=i % 4,
                sticky="w",
                padx=(0, 18),
                pady=2,
            )
        row += 1

        igv = CollapsibleSection(canvas, "Genome Browser Tracks", open_by_default=False)
        igv.grid(row=row, column=0, sticky="ew")
        igv.body.columnconfigure(1, weight=1)
        ttk.Checkbutton(igv.body, text="Enable genome browser", variable=self.enable_igv).grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 8)
        )
        self._entry_row(igv.body, 1, "Genome", self.genome)
        ttk.Label(igv.body, text="Bucket URLs, one per line").grid(row=2, column=0, sticky="nw", pady=4)
        self.bucket_text = Text(igv.body, height=4, wrap="word")
        self.bucket_text.grid(row=2, column=1, sticky="ew", pady=4)
        self._entry_row(igv.body, 3, "Track names", self.atac_names)
        self._entry_row(igv.body, 4, "Track max heights", self.max_heights)
        row += 1

        visual = CollapsibleSection(canvas, "Colors", open_by_default=False)
        visual.grid(row=row, column=0, sticky="ew")
        visual.body.columnconfigure(0, weight=1)
        ttk.Label(visual.body, text="Global color palette, separated by comma or new line").grid(
            row=0, column=0, sticky="w"
        )
        self.colors_text = Text(visual.body, height=4, wrap="word")
        self.colors_text.grid(row=1, column=0, sticky="ew", pady=(4, 0))
        self.colors_text.insert("1.0", "\n".join(DEFAULT_COLORS))
        row += 1

        runtime = CollapsibleSection(canvas, "Runtime Settings", open_by_default=False)
        runtime.grid(row=row, column=0, sticky="ew")
        runtime.body.columnconfigure(1, weight=1)
        self._entry_row(runtime.body, 0, "Host", self.host)
        self._entry_row(runtime.body, 1, "Port", self.port)
        self._entry_row(runtime.body, 2, "Max cells", self.max_cells)
        ttk.Checkbutton(runtime.body, text="Lazy load data", variable=self.lazy_load).grid(
            row=3, column=1, sticky="w", pady=4
        )
        ttk.Label(runtime.body, text="Backed mode").grid(row=4, column=0, sticky="w", pady=4)
        ttk.Combobox(
            runtime.body,
            textvariable=self.backed_mode,
            values=("false", "true", "r+"),
            state="readonly",
            width=14,
        ).grid(row=4, column=1, sticky="w", pady=4)
        ttk.Label(runtime.body, text="Embedding backend").grid(row=5, column=0, sticky="w", pady=4)
        ttk.Combobox(
            runtime.body,
            textvariable=self.embedding_backend,
            values=("scattergl", "datashader"),
            state="readonly",
            width=18,
        ).grid(row=5, column=1, sticky="w", pady=4)
        row += 1

        actions = ttk.Frame(canvas)
        actions.grid(row=row, column=0, sticky="ew", pady=(18, 0))
        ttk.Button(actions, text="Preview JSON", command=self._preview).pack(side="left")
        ttk.Button(actions, text="Save Config", command=lambda: self._save(launch=False)).pack(side="right")
        ttk.Button(actions, text="Save And Launch GUANACO", command=lambda: self._save(launch=True)).pack(
            side="right", padx=(0, 8)
        )
        canvas.columnconfigure(0, weight=1)

    def _entry_row(self, parent, row: int, label: str, variable: StringVar, *, browse=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        entry = ttk.Entry(parent, textvariable=variable)
        entry.grid(row=row, column=1, sticky="ew", pady=4)
        if browse is not None:
            ttk.Button(parent, text="Browse", command=browse).grid(row=row, column=2, sticky="e", padx=(8, 0))

    def _browse_config(self):
        path = filedialog.asksaveasfilename(
            title="Save GUANACO config",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if path:
            self.config_path.set(path)

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

    @staticmethod
    def _split_values(raw: str) -> list[str]:
        values: list[str] = []
        for chunk in raw.replace(",", "\n").splitlines():
            value = chunk.strip()
            if value:
                values.append(value)
        return values

    @staticmethod
    def _parse_bool_or_backed(value: str):
        normalized = value.strip().lower()
        if normalized == "r+":
            return "r+"
        return normalized == "true"

    def _build_config(self) -> dict:
        dataset_name = self.dataset_name.get().strip()
        sc_data = self.sc_data.get().strip()
        if not dataset_name:
            raise ValueError("Dataset name is required.")
        if not sc_data:
            raise ValueError("AnnData / MuData file is required.")
        sc_data_path = Path(sc_data).expanduser().resolve()
        if not sc_data_path.exists():
            raise ValueError(f"AnnData / MuData file does not exist: {sc_data_path}")

        dataset = {
            "sc_data": str(sc_data_path),
        }
        description = self.description.get().strip()
        if description:
            dataset["description"] = description

        markers = self._split_values(self.markers_text.get("1.0", "end"))
        if markers:
            dataset["markers"] = markers

        optional_plots = [value for _, value in OPTIONAL_PLOTS if self.plot_vars[value].get()]
        dataset["optional_plot_components"] = optional_plots

        if self.enable_igv.get():
            bucket_urls = self._split_values(self.bucket_text.get("1.0", "end"))
            if not bucket_urls:
                raise ValueError("Bucket URLs are required when genome browser is enabled.")
            dataset["genome"] = self.genome.get().strip() or "hg38"
            dataset["bucket_urls"] = bucket_urls

            names = self._split_values(self.atac_names.get())
            if names:
                if len(names) != len(bucket_urls):
                    raise ValueError("Track names must have the same count as bucket URLs.")
                dataset["ATAC_name"] = names
            heights = self._split_values(self.max_heights.get())
            if heights:
                if len(heights) != len(bucket_urls):
                    raise ValueError("Track max heights must have the same count as bucket URLs.")
                dataset["max_height"] = [int(value) for value in heights]

        settings = {
            "host": self.host.get().strip() or "0.0.0.0",
            "port": int(self.port.get().strip()),
            "max_cells": None if self.max_cells.get().strip().lower() in {"", "none", "null"} else int(self.max_cells.get()),
            "lazy_load": bool(self.lazy_load.get()),
            "backed_mode": self._parse_bool_or_backed(self.backed_mode.get()),
            "embedding_render_backend": self.embedding_backend.get(),
        }

        config = {
            dataset_name: dataset,
            "title": self.title.get().strip() or "GUANACO",
            "settings": settings,
        }
        colors = self._split_values(self.colors_text.get("1.0", "end"))
        if colors:
            config["color"] = colors
        return config

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
