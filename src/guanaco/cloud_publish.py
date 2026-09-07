"""Portable, reviewed Plotly Cloud exports. No GUI or cloud SDK at import time.

Preparing an export is entirely local. Only ``run_plotly`` can contact Plotly;
callers must obtain upload confirmation before invoking ``app publish``.
"""

from __future__ import annotations

import copy
from contextlib import ExitStack
from dataclasses import dataclass
from email.parser import BytesParser
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import tomllib
from urllib.parse import urlsplit
import zipfile


MAX_BUNDLE_BYTES = 200 * 1024 * 1024
MAX_BUNDLE_FILES = 100_000
_GLOBAL_KEYS = {"title", "color", "genome", "settings"}
_ANNOTATIONS = {"gene_annotation", "gene_annotation_path", "gtf_path", "gtf", "genome"}
_SECRET_KEY = re.compile(
    r"password|secret|token|credential|api[_-]?key|access[_-]?key", re.I
)

# Self-contained so a bundle can use a released GUANACO wheel, too. Resolve only
# recorded local paths; the existing loader still requires absolute sc_data paths.
APP_SOURCE = '''"""GUANACO dashboard exported by the configuration wizard."""
import json
import os
from pathlib import Path
import tempfile
import zipfile

_root = Path(__file__).resolve().parent
_runtime = tempfile.TemporaryDirectory(prefix="guanaco-runtime-")
_config = json.loads((_root / "guanaco.json").read_text(encoding="utf-8"))
for _keys in json.loads((_root / "local-paths.json").read_text(encoding="utf-8")):
    _node = _config
    for _key in _keys[:-1]:
        _node = _node[_key]
    _path = (_root / _node[_keys[-1]]).resolve()
    if not _path.is_relative_to(_root) or not _path.exists():
        raise ValueError("Missing or unsafe bundled data path")
    if _path.name.endswith(".zarr.zip"):
        _unpacked = Path(_runtime.name) / _path.relative_to(_root).with_suffix("")
        if not _unpacked.exists():
            with zipfile.ZipFile(_path) as _store:
                for _member in _store.namelist():
                    if not (_unpacked / _member).resolve().is_relative_to(_unpacked.resolve()):
                        raise ValueError("Unsafe Zarr archive member")
                _store.extractall(_unpacked)
        _path = _unpacked
    _node[_keys[-1]] = str(_path)
# Keep the runtime config alive for this worker, without altering the export.
_resolved = Path(_runtime.name) / "guanaco.json"
_resolved.write_text(json.dumps(_config), encoding="utf-8")
os.environ["GUANACO_CONFIG"] = str(_resolved)

from guanaco.main import app

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", "8050")), debug=False)
'''


class Cancelled(RuntimeError):
    """A local operation was cancelled; remote work may already have started."""


def _check_cancel(cancel):
    if cancel is not None and cancel.is_set():
        raise Cancelled(
            "Cancelled. If uploading had started, check Plotly Cloud before retrying."
        )


@dataclass(frozen=True)
class BundleFile:
    name: str
    size: int
    sha256: str


@dataclass(frozen=True)
class PreparedBundle:
    directory: Path
    files: tuple[BundleFile, ...]
    remote_sources: tuple[str, ...]
    app_id: str = ""

    @property
    def size(self):
        return sum(file.size for file in self.files)

    def review(self):
        action = (
            f"UPDATE existing app {self.app_id}" if self.app_id else "CREATE a new app"
        )
        lines = [
            action,
            f"Export: {self.directory}",
            f"Upload: {len(self.files)} files, {self.size / 1024**2:.2f} MiB",
            "",
        ]
        lines.extend(
            f"{file.size / 1024**2:8.3f} MiB  {file.name}" for file in self.files
        )
        if self.remote_sources:
            lines.extend(
                [
                    "",
                    "Remote references (not copied; cloud server needs access):",
                    *self.remote_sources,
                ]
            )
        lines.extend(
            [
                "",
                "Whole selected data files are included, including annotations and images.",
                "Local sharing passwords and credentials are not exported.",
                "Review consent and data-governance requirements before uploading.",
            ]
        )
        return "\n".join(lines)


def _files(root: Path):
    """Walk without following symlinks or accidentally including special files."""
    count = 0
    for directory, dirs, names in os.walk(root):
        for name in dirs + names:
            path = Path(directory) / name
            if path.is_symlink():
                raise ValueError(f"Symlinks cannot be exported: {path.name}")
        for name in sorted(names):
            path = Path(directory) / name
            if not path.is_file():
                raise ValueError(f"Not a regular file: {path.name}")
            count += 1
            if count > MAX_BUNDLE_FILES:
                raise ValueError(
                    "Too many files for this export. Use a remote Zarr store instead."
                )
            yield path


def _snapshot(root: Path, cancel=None):
    result = []
    total = 0
    for path in sorted(_files(root)):
        _check_cancel(cancel)
        size = path.stat().st_size
        total += size
        if total > MAX_BUNDLE_BYTES:
            raise ValueError(
                "Export exceeds the 200 MiB publishing limit. Use a remote Zarr "
                "store or temporary local sharing; data will not be subsampled."
            )
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                _check_cancel(cancel)
                digest.update(block)
        result.append(
            BundleFile(path.relative_to(root).as_posix(), size, digest.hexdigest())
        )
    return tuple(result)


def verify_bundle(bundle: PreparedBundle, cancel=None):
    if _snapshot(bundle.directory, cancel) != bundle.files:
        raise ValueError(
            "Export changed after review. Prepare and review a fresh export before publishing."
        )
    if (
        deployment_info(bundle.directory / "plotly-cloud.toml").get("app_id", "")
        != bundle.app_id
    ):
        raise ValueError("Deployment target changed after review.")


def deployment_info(path: Path):
    """Read only official app identity fields, never arbitrary TOML instructions."""
    if not path.exists():
        return {}
    if path.is_symlink() or path.stat().st_size > 64 * 1024:
        raise ValueError("Invalid Plotly Cloud configuration file.")
    with path.open("rb") as stream:
        config = tomllib.load(stream)
    result = {
        key: str(config[key])
        for key in ("name", "app_id", "app_url")
        if config.get(key)
    }
    for key in ("app_id", "app_url"):
        if key in result and not re.fullmatch(r"[A-Za-z0-9_-]+", result[key]):
            raise ValueError(f"Invalid {key} in Plotly Cloud configuration.")
    return result


def deployment_links(directory: Path):
    info = deployment_info(directory / "plotly-cloud.toml")
    dashboard = f"https://{info['app_url']}.plotly.app" if info.get("app_url") else ""
    settings = (
        f"https://cloud.plotly.com/app/{info['app_id']}/settings"
        if info.get("app_id")
        else "https://cloud.plotly.com"
    )
    return dashboard, settings


def _check_config(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if _SECRET_KEY.search(key) and item:
                raise ValueError(
                    f"Remove credential field '{key}' before exporting. "
                    "Configure private-data access in Plotly Cloud instead."
                )
            _check_config(item)
    elif isinstance(value, list):
        for item in value:
            _check_config(item)
    elif isinstance(value, str) and "://" in value:
        parsed = urlsplit(value)
        if parsed.username or parsed.password or parsed.query or parsed.fragment:
            raise ValueError(
                "URLs containing credentials, query parameters or fragments cannot "
                "be exported. Use a clean URL and configure cloud-side credentials."
            )


def source_checkout():
    root = Path(__file__).resolve().parents[2]
    project = root / "pyproject.toml"
    if project.is_file() and (root / "src/guanaco/cloud_publish.py").is_file():
        with project.open("rb") as stream:
            if tomllib.load(stream).get("project", {}).get("name") == "guanaco-viz":
                return root
    return None


def _wheel_version(path: Path):
    with zipfile.ZipFile(path) as wheel:
        manifests = [
            name for name in wheel.namelist() if name.endswith(".dist-info/METADATA")
        ]
        if len(manifests) != 1:
            raise ValueError("Select a GUANACO wheel with one package metadata record.")
        info = BytesParser().parsebytes(wheel.read(manifests[0]))
    if re.sub(r"[-_.]+", "-", info.get("Name", "")).lower() != "guanaco-viz":
        raise ValueError("The selected wheel is not guanaco-viz.")
    if not re.fullmatch(r"[A-Za-z0-9.!+_-]+", info.get("Version", "")):
        raise ValueError("Invalid GUANACO wheel version.")
    # GUANACO itself is pure Python; reject wheels specific to the developer's OS.
    if not path.name.endswith("-none-any.whl"):
        raise ValueError("Use a platform-independent GUANACO wheel (*-none-any.whl).")
    if not re.fullmatch(r"guanaco_viz-[A-Za-z0-9.+_!-]+-none-any\.whl", path.name):
        raise ValueError("Invalid GUANACO wheel filename.")
    return info["Version"]


def prepare_bundle(
    config: dict,
    parent: Path,
    *,
    base_dir: Path,
    wheel: Path | None = None,
    existing: Path | None = None,
    cancel=None,
    emit=lambda text: None,
):
    """Copy only referenced data/resources into a NEW folder; never upload.

    A source checkout builds its current wheel. An installed release pins the
    installed version on PyPI, unless an explicit wheel was selected.
    """
    from guanaco.pages.matrix.plots.gene_annotation import is_known_genome_id

    cfg = copy.deepcopy(config)
    settings = cfg.setdefault("settings", {})
    for key in ("share", "share_username", "share_password", "host", "port"):
        settings.pop(key, None)
    _check_config(cfg)
    identity = deployment_info(existing) if existing else {}
    if existing and not identity.get("app_id"):
        raise ValueError("The selected configuration has no existing app_id.")
    sources = {}
    local_paths = []
    remote_sources = set()
    extras = set()

    def reference(node, key, keys, *, data=False):
        value = node.get(key)
        if not value:
            return
        if not isinstance(value, str):
            raise ValueError(f"{key} must be a path or URL.")
        if "://" in value:
            remote_sources.add(value)
            extras.add("cloud")
            return
        if not data and is_known_genome_id(value):
            return
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = base_dir / path
        if path.is_symlink():
            raise ValueError(f"Symlinks cannot be exported: {path.name}")
        path = path.resolve(strict=True)
        if data and path.suffix.lower() not in {".h5ad", ".h5mu", ".zarr"}:
            raise ValueError("Only .h5ad, .h5mu and .zarr data can be bundled.")
        if not data and not (
            path.name.lower().endswith(
                (".gtf", ".gtf.gz", ".gff", ".gff3", ".gff.gz", ".gff3.gz")
            )
        ):
            raise ValueError("Local gene annotations must be GTF or GFF files.")
        if path.is_dir():
            if path.suffix.lower() != ".zarr" or not any(
                (path / marker).is_file() for marker in (".zgroup", "zarr.json")
            ):
                raise ValueError("Only validated Zarr directories can be bundled.")
            extras.add("cloud")
        if path not in sources:
            sources[path] = f"data/{len(sources):03d}/{path.name}" + (
                ".zip" if path.is_dir() else ""
            )
        node[key] = sources[path]
        local_paths.append(keys + [key])

    if cfg.get("genome"):
        reference(cfg, "genome", [])
    for name, dataset in cfg.items():
        if name in _GLOBAL_KEYS:
            continue
        if not isinstance(dataset, dict) or not dataset.get("sc_data"):
            raise ValueError(f"Dataset '{name}' requires sc_data for cloud export.")
        reference(dataset, "sc_data", [name], data=True)
        for key in _ANNOTATIONS:
            reference(dataset, key, [name])
        for modality, block in dataset.get("modalities", {}).items():
            for key in _ANNOTATIONS:
                reference(block, key, [name, "modalities", modality])
        if dataset.get("bucket_urls"):
            extras.add("tracks")
            remote_sources.update(dataset["bucket_urls"])
    if settings.get("embedding_render_backend") == "datashader":
        extras.add("datashader")

    # Inventory and bound the copy before building or creating the export folder.
    copies = []
    total = 0
    for source, target in sources.items():
        for file in _files(source) if source.is_dir() else [source]:
            _check_cancel(cancel)
            relative = file.relative_to(source) if source.is_dir() else Path()
            if any(
                part.lower() in {".env", ".git", ".aws", "credentials"}
                or part.lower().startswith(".env.")
                for part in relative.parts
            ):
                raise ValueError(
                    "A data directory contains credential or project files; clean it first."
                )
            total += file.stat().st_size
            if total > MAX_BUNDLE_BYTES or len(copies) >= MAX_BUNDLE_FILES:
                raise ValueError(
                    "Selected data exceed the 200 MiB / 100,000-file export limit. "
                    "Use a remote Zarr store or temporary local sharing."
                )
            copies.append((file, Path(target), relative, source.is_dir()))

    parent = Path(parent).expanduser().resolve(strict=True)
    destination = Path(tempfile.mkdtemp(prefix="guanaco-cloud-", dir=parent))
    emit(f"Preparing local export: {destination}")
    try:
        with ExitStack() as stack:
            archives = {}
            for source, target, relative, archive in copies:
                _check_cancel(cancel)
                target = destination / target
                target.parent.mkdir(parents=True, exist_ok=True)
                if archive:
                    # Plotly's default excludes drop Zarr's `var/` directory.
                    # Archive once per store, preserving every chunk for startup.
                    if target not in archives:
                        archives[target] = stack.enter_context(
                            zipfile.ZipFile(target, "w", compression=zipfile.ZIP_STORED)
                        )
                    archives[target].write(source, relative.as_posix())
                else:
                    shutil.copyfile(source, target)
        suffix = f"[{','.join(sorted(extras))}]" if extras else ""
        checkout = source_checkout()
        if wheel is None and checkout is not None:
            emit("Building the current GUANACO wheel (no PyPI publication required)…")
            run_process(
                [
                    sys.executable,
                    "-m",
                    "build",
                    "--wheel",
                    "--no-isolation",
                    "--outdir",
                    str(destination),
                    str(checkout),
                ],
                cancel=cancel,
                emit=emit,
            )
            wheels = list(destination.glob("*.whl"))
            if len(wheels) != 1:
                raise ValueError("Wheel build did not produce exactly one wheel.")
            wheel = wheels[0]
        if wheel is not None:
            wheel = Path(wheel).expanduser().resolve(strict=True)
            _wheel_version(wheel)
            if wheel.parent != destination:
                shutil.copyfile(wheel, destination / wheel.name)
            requirement = f"./{wheel.name}{suffix}"
        else:
            version = metadata.version("guanaco-viz")
            requirement = f"guanaco-viz{suffix}=={version}"
        (destination / "requirements.txt").write_text(
            requirement + "\n", encoding="utf-8"
        )
        (destination / "app.py").write_text(APP_SOURCE, encoding="utf-8")
        (destination / "guanaco.json").write_text(
            json.dumps(cfg, indent=2), encoding="utf-8"
        )
        (destination / "local-paths.json").write_text(
            json.dumps(local_paths), encoding="utf-8"
        )
        if identity:
            (destination / "plotly-cloud.toml").write_text(
                "\n".join(
                    f"{key} = {json.dumps(value, ensure_ascii=False)}"
                    for key, value in identity.items()
                )
                + "\n",
                encoding="utf-8",
            )
        return PreparedBundle(
            destination,
            _snapshot(destination, cancel),
            tuple(sorted(remote_sources)),
            identity.get("app_id", ""),
        )
    except Exception as exc:
        # Keep the exact new directory recoverable; never remove the chosen parent.
        raise RuntimeError(
            f"{exc}\nIncomplete export retained at: {destination}"
        ) from exc


def run_process(command, *, cancel=None, emit=lambda text: None, timeout=600):
    """Drain output off-thread so silent subprocesses remain cancellable."""
    _check_cancel(cancel)
    events = queue.Queue()
    env = {**os.environ, "NO_COLOR": "1", "TERM": "dumb", "PYTHONUNBUFFERED": "1"}
    with subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    ) as process:

        def read():
            for line in process.stdout:
                events.put(line.rstrip())
            events.put(None)

        reader = threading.Thread(target=read, daemon=True)
        reader.start()
        deadline = time.monotonic() + timeout
        tail = []
        try:
            while True:
                _check_cancel(cancel)
                if time.monotonic() > deadline:
                    raise TimeoutError(
                        "Operation timed out. Check Plotly Cloud before retrying an upload."
                    )
                try:
                    line = events.get(timeout=0.1)
                except queue.Empty:
                    continue
                if line is None:
                    break
                emit(line)
                tail.append(line)
                tail = tail[-80:]
            code = process.wait(timeout=5)
            if code:
                raise RuntimeError(
                    f"Command failed (exit {code}).\n" + "\n".join(tail[-15:])
                )
            return "\n".join(tail)
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            reader.join(timeout=1)


def run_plotly(arguments, **kwargs):
    try:
        version = metadata.version("plotly-cloud")
    except metadata.PackageNotFoundError as exc:
        raise RuntimeError(
            "Publishing tools are not installed. In GUANACO's Python environment run:\n"
            'python -m pip install "plotly-cloud>=0.4.3,<0.5"\n'
            'Future installs can use: pip install "guanaco-viz[publish]"'
        ) from exc
    if not re.fullmatch(r"0\.4\.(?:[3-9]|[1-9][0-9]+)", version):
        raise RuntimeError(
            f"Unsupported plotly-cloud {version}; install plotly-cloud>=0.4.3,<0.5."
        )
    # Use this interpreter, not an unrelated `plotly` executable on PATH.
    # Run this stdlib-only file directly to avoid importing GUANACO's plotting
    # stack in the CLI subprocess. The entrypoint also preserves real CLI errors.
    return run_process(
        [sys.executable, str(Path(__file__).resolve()), *arguments],
        **kwargs,
    )


def _plotly_cli_main():
    """Work around 0.4.3's error handler reparsing command flags with no schema.

    The initial parse succeeds, but a later API/build/etc. exception invokes
    ``parse_args([])`` and masks the real error as an unknown command flag. Reuse
    the successfully parsed arguments for that error handler only. Authentication,
    publishing and error rendering remain owned by Plotly's CLI.
    """
    from plotly_cloud import cli

    original_parse = cli.parse_args
    parsed = None

    def parse_once(arguments, args_index=3):
        nonlocal parsed
        if not arguments and parsed is not None:
            return parsed
        parsed = original_parse(arguments, args_index)
        return parsed

    cli.parse_args = parse_once
    try:
        cli.main()
    finally:
        cli.parse_args = original_parse


def publish_bundle(bundle, name, team="", **kwargs):
    if not name.strip():
        raise ValueError("An application name is required.")
    verify_bundle(bundle, kwargs.get("cancel"))
    args = [
        "app",
        "publish",
        "--project-path",
        str(bundle.directory),
        "--poll-timeout",
        "180",
        "--poll-interval",
        "3",
    ]
    if not bundle.app_id:
        args += ["--name", name.strip(), "--entrypoint-module", "app"]
        if team.strip():
            args += ["--team", team.strip()]
    return run_plotly(args, **kwargs)


if __name__ == "__main__":
    _plotly_cli_main()
