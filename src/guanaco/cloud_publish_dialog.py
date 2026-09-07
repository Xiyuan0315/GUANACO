"""Tk front end for the local export / explicit cloud publish workflow."""

from functools import partial
from pathlib import Path
import queue
import threading
from tkinter import StringVar, Text, Toplevel, filedialog, messagebox, ttk
import webbrowser

from guanaco.cloud_publish import (
    deployment_links,
    prepare_bundle,
    publish_bundle,
    run_plotly,
    source_checkout,
)


class CloudPublishDialog:
    def __init__(self, parent, config, base_dir):
        self.window = window = Toplevel(parent)
        window.title("GUANACO · Publish to Plotly Cloud")
        window.geometry("780x720")
        window.minsize(780, 650)
        window.transient(parent)
        window.grab_set()
        window.protocol("WM_DELETE_WINDOW", self._close)
        self.config = config
        self.base_dir = Path(base_dir)
        self.bundle = None
        self.busy = False
        self.attempted = False
        self.preparation_inputs = None
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.name = StringVar(window, value=config.get("title") or "GUANACO dashboard")
        self.team = StringVar(window)
        self.existing = StringVar(window)
        self.wheel = StringVar(window)
        self.status = StringVar(
            window, value="Prepare an export first. Nothing is uploaded automatically."
        )

        frame = ttk.Frame(window, padding=16)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)
        ttk.Label(
            frame, text="Publish your dashboard", font=("TkDefaultFont", 16, "bold")
        ).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        ttk.Label(
            frame,
            text="1. Prepare and review files   →   2. Sign in   →   3. Confirm upload",
            wraplength=680,
        ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 12))
        self.controls = []
        for row, label, variable in [
            (2, "App name", self.name),
            (3, "Team (blank = default)", self.team),
            (4, "Update existing app", self.existing),
            (5, "GUANACO wheel (optional)", self.wheel),
        ]:
            ttk.Label(frame, text=label).grid(
                row=row, column=0, sticky="w", padx=(0, 10), pady=4
            )
            entry = ttk.Entry(frame, textvariable=variable)
            entry.grid(row=row, column=1, sticky="ew", pady=4)
            self.controls.append(entry)
        for row, variable, title, patterns in [
            (
                4,
                self.existing,
                "Select the previous export's plotly-cloud.toml",
                [("Plotly config", "*.toml")],
            ),
            (
                5,
                self.wheel,
                "Select a built GUANACO wheel",
                [("Python wheel", "*.whl")],
            ),
        ]:
            button = ttk.Button(
                frame,
                text="Browse…",
                command=partial(self._browse, variable, title, patterns),
            )
            button.grid(row=row, column=2, padx=(8, 0))
            self.controls.append(button)
        hint = (
            "Source checkout detected: preparation builds your current code into a wheel."
            if source_checkout()
            else "Installed release: preparation pins your GUANACO version on PyPI."
        )
        ttk.Label(
            frame,
            text=hint + " Leave Update blank to create a new app.",
            wraplength=710,
        ).grid(row=6, column=0, columnspan=3, sticky="w", pady=(6, 12))
        toolbar = ttk.Frame(frame)
        toolbar.grid(row=7, column=0, columnspan=3, sticky="ew")
        for label, command in [
            ("Prepare export…", self._prepare),
            (
                "Sign in…",
                lambda: self._command(
                    ["user", "login"], "Complete sign-in in your browser…"
                ),
            ),
            ("Account / teams", self._account),
        ]:
            button = ttk.Button(toolbar, text=label, command=command)
            button.pack(side="left", padx=(0, 6))
            self.controls.append(button)
        self.publish = ttk.Button(
            toolbar, text="Publish…", command=self._publish, state="disabled"
        )
        self.publish.pack(side="right")

        self.output = Text(frame, wrap="word", height=16, state="disabled")
        self.output.grid(row=8, column=0, columnspan=2, sticky="nsew", pady=12)
        scroll = ttk.Scrollbar(frame, command=self.output.yview)
        scroll.grid(row=8, column=2, sticky="ns", pady=12)
        self.output.configure(yscrollcommand=scroll.set)
        frame.rowconfigure(8, weight=1)
        ttk.Label(frame, textvariable=self.status, wraplength=710).grid(
            row=9, column=0, columnspan=3, sticky="w"
        )
        footer = ttk.Frame(frame)
        footer.grid(row=10, column=0, columnspan=3, sticky="ew", pady=(12, 0))
        self.links = []
        for label, command in [
            ("Open dashboard", self._open),
            ("Copy link", self._copy),
            ("Manage access / compute", self._manage),
            ("Status", lambda: self._app_command("status")),
            ("Build logs", lambda: self._app_command("logs")),
        ]:
            button = ttk.Button(footer, text=label, command=command, state="disabled")
            button.pack(side="left", padx=(0, 4))
            self.links.append(button)
        self.cancel_button = ttk.Button(
            frame, text="Cancel operation", command=self.cancel.set, state="disabled"
        )
        self.cancel_button.grid(row=11, column=0, sticky="w", pady=(10, 0))
        ttk.Label(
            frame,
            text="Only upload data you are permitted to host. Access and resource limits are managed in Plotly Cloud.",
            wraplength=710,
        ).grid(row=12, column=0, columnspan=3, sticky="w", pady=(10, 0))
        # Tk reads and writes happen only in this main-thread queue consumer.
        self.timer = window.after(100, self._poll)
        for variable in (self.wheel, self.existing):
            variable.trace_add("write", lambda *_: self._refresh())

    def _browse(self, variable, title, patterns):
        path = filedialog.askopenfilename(
            parent=self.window, title=title, filetypes=patterns
        )
        if path:
            variable.set(path)
            self.bundle = None
            self._refresh()

    def _log(self, text):
        self.output.configure(state="normal")
        self.output.insert("end", text + "\n")
        self.output.see("end")
        self.output.configure(state="disabled")

    def _refresh(self):
        for control in self.controls:
            control.configure(state="disabled" if self.busy else "normal")
        self.cancel_button.configure(state="normal" if self.busy else "disabled")
        unchanged = self.preparation_inputs == (self.wheel.get(), self.existing.get())
        self.publish.configure(
            state="normal"
            if self.bundle and unchanged and not self.busy and not self.attempted
            else "disabled"
        )
        try:
            dashboard, _ = (
                deployment_links(self.bundle.directory) if self.bundle else ("", "")
            )
        except (ValueError, OSError):
            dashboard = ""
        deployed = bool(
            self.bundle and (self.bundle.directory / "plotly-cloud.toml").exists()
        )
        for index, button in enumerate(self.links):
            available = bool(dashboard) if index < 2 else deployed
            button.configure(
                state="normal" if available and not self.busy else "disabled"
            )

    def _start(self, operation, label, kind="command"):
        if self.busy:
            return
        self.busy = True
        self.cancel.clear()
        self.status.set(label)
        self._refresh()
        # Capture only plain Python objects: workers never own or touch Tk variables.
        events, cancel = self.events, self.cancel

        def work():
            try:
                result = operation(
                    cancel=cancel, emit=lambda line: events.put(("log", line))
                )
                events.put((kind, result))
            except Exception as exc:
                events.put(("error", str(exc)))

        threading.Thread(target=work, daemon=True).start()

    def _poll(self):
        try:
            while True:
                kind, result = self.events.get_nowait()
                if kind == "log":
                    self._log(result)
                    continue
                self.busy = False
                if kind == "prepared":
                    self.bundle = result
                    self.attempted = False
                    self._log("\n" + result.review())
                    self.status.set(
                        "Ready for review. Nothing has been uploaded. Keep this export folder for future updates."
                    )
                elif kind == "published":
                    self.status.set(
                        "Publishing command finished. Check status below or in Plotly Cloud; completion is not a guarantee the app is running."
                    )
                elif kind == "error":
                    self._log(result)
                    self.status.set(
                        "Operation did not complete. Details are shown above."
                    )
                    messagebox.showerror("Plotly Cloud", result, parent=self.window)
                else:
                    self.status.set("Command finished. See details above.")
                if self.attempted and self.bundle and kind in {"published", "error"}:
                    target = self.bundle.directory / "plotly-cloud.toml"
                    if target.exists():
                        self.existing.set(str(target))
                    self._log(
                        "Before another upload, prepare a fresh export. Preserve plotly-cloud.toml to update this app instead of creating another."
                    )
                self._refresh()
        except queue.Empty:
            pass
        self.timer = self.window.after(100, self._poll)

    def _prepare(self):
        parent = filedialog.askdirectory(
            parent=self.window, title="Choose a parent folder for a NEW export"
        )
        if not parent:
            return
        self.bundle = None
        self.preparation_inputs = (self.wheel.get(), self.existing.get())
        self._start(
            partial(
                prepare_bundle,
                self.config,
                Path(parent),
                base_dir=self.base_dir,
                wheel=Path(self.wheel.get()) if self.wheel.get().strip() else None,
                existing=Path(self.existing.get())
                if self.existing.get().strip()
                else None,
            ),
            "Preparing locally; no upload…",
            "prepared",
        )

    def _command(self, args, label):
        self._start(partial(run_plotly, args), label)

    def _account(self):
        def account(**kwargs):
            run_plotly(["user", "whoami"], **kwargs)
            return run_plotly(["user", "teams"], **kwargs)

        self._start(account, "Checking account and available teams…")

    def _publish(self):
        if not self.bundle or self.busy or self.attempted:
            return
        if not self.name.get().strip():
            messagebox.showerror(
                "App name", "Enter an app name first.", parent=self.window
            )
            return
        action = (
            f"Update app {self.bundle.app_id}"
            if self.bundle.app_id
            else f"Create '{self.name.get().strip()}' (team: {self.team.get().strip() or 'default'})"
        )
        if not messagebox.askyesno(
            "Confirm upload",
            f"{action}?\n\nUpload the reviewed {len(self.bundle.files)} files "
            f"({self.bundle.size / 1024**2:.2f} MiB) to Plotly Cloud, including whole selected datasets.\n\n"
            "Confirm you have permission to host these data. New apps use Plotly's default access settings; "
            "updates retain existing access. Review sharing settings in Plotly Cloud.",
            parent=self.window,
        ):
            return
        self.attempted = True
        self._start(
            partial(publish_bundle, self.bundle, self.name.get(), self.team.get()),
            "Uploading and checking build status…",
            "published",
        )

    def _app_command(self, command):
        args = ["app", command, "--project-path", str(self.bundle.directory)]
        if command == "logs":
            args += ["--type", "build"]
        self._command(args, f"Retrieving {command}…")

    def _open(self):
        webbrowser.open(deployment_links(self.bundle.directory)[0])

    def _copy(self):
        self.window.clipboard_clear()
        self.window.clipboard_append(deployment_links(self.bundle.directory)[0])

    def _manage(self):
        webbrowser.open(deployment_links(self.bundle.directory)[1])

    def _close(self):
        if self.busy:
            self.cancel.set()
            self.status.set(
                "Cancelling… Wait for the operation to stop, then close. An upload cannot be rolled back locally."
            )
            return
        self.window.after_cancel(self.timer)
        self.window.grab_release()
        self.window.destroy()
