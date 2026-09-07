"""Publishing tests never authenticate, contact Plotly, or upload data."""

import copy
import json
from pathlib import Path
import runpy
import sys
import threading
import types
import zipfile

import pytest

from guanaco import cloud_publish as cloud


@pytest.fixture
def sample(tmp_path, monkeypatch):
    source = tmp_path / "cells.h5ad"
    source.write_bytes(b"test data")
    config = {
        "Cells": {"sc_data": str(source)},
        "title": "Test dashboard",
        "settings": {
            "share": True,
            "share_username": "private-user",
            "share_password": "secret!",
            "backed_mode": True,
            "max_cells": None,
        },
    }
    monkeypatch.setattr(cloud, "source_checkout", lambda: None)
    original_version = cloud.metadata.version
    monkeypatch.setattr(
        cloud.metadata,
        "version",
        lambda name: "1.2.3" if name == "guanaco-viz" else original_version(name),
    )
    return config


def prepare(tmp_path, config, **kwargs):
    return cloud.prepare_bundle(config, tmp_path, base_dir=tmp_path, **kwargs)


def test_portable_allowlisted_export_strips_sharing_credentials(tmp_path, sample):
    before = copy.deepcopy(sample)
    (tmp_path / ".env").write_text("DO_NOT_UPLOAD=secret")
    (tmp_path / "unrelated.txt").write_text("DO_NOT_UPLOAD")
    bundle = prepare(tmp_path, sample)
    assert sample == before
    assert bundle.directory != tmp_path
    names = {file.name for file in bundle.files}
    assert names == {
        "app.py",
        "guanaco.json",
        "local-paths.json",
        "requirements.txt",
        "data/000/cells.h5ad",
    }
    config = json.loads((bundle.directory / "guanaco.json").read_text())
    assert config["Cells"]["sc_data"] == "data/000/cells.h5ad"
    assert config["settings"] == {"backed_mode": True, "max_cells": None}
    assert (bundle.directory / "requirements.txt").read_text() == "guanaco-viz==1.2.3\n"
    assert "secret!" not in bundle.review()
    cloud.verify_bundle(bundle)


def test_sources_are_deduplicated_and_annotations_relocated(tmp_path, sample):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("test genes")
    sample["Cells"]["modalities"] = {"rna": {"gene_annotation": "genes.gtf"}}
    sample["Other"] = {
        "sc_data": sample["Cells"]["sc_data"],
        "gene_annotation_path": str(gtf),
    }
    sample["genome"] = "hg38"
    bundle = prepare(tmp_path, sample)
    config = json.loads((bundle.directory / "guanaco.json").read_text())
    assert config["Cells"]["sc_data"] == config["Other"]["sc_data"]
    assert (
        config["Cells"]["modalities"]["rna"]["gene_annotation"]
        == config["Other"]["gene_annotation_path"]
    )
    assert config["genome"] == "hg38"
    assert sum(file.name.endswith(".h5ad") for file in bundle.files) == 1


@pytest.mark.parametrize(
    "url",
    [
        "https://example.org/cells.zarr",
        "s3://bucket/cells.zarr",
        "https://example.org/cells.h5ad",
    ],
)
def test_remote_references_are_not_downloaded(tmp_path, sample, url):
    sample["Cells"]["sc_data"] = url
    sample["Cells"]["bucket_urls"] = ["s3://tracks/bigwigs/"]
    sample["settings"]["embedding_render_backend"] = "datashader"
    bundle = prepare(tmp_path, sample)
    assert not (bundle.directory / "data").exists()
    assert url in bundle.remote_sources
    assert (
        bundle.directory / "requirements.txt"
    ).read_text() == "guanaco-viz[cloud,datashader,tracks]==1.2.3\n"


@pytest.mark.parametrize(
    "url",
    [
        "https://user:secret@example.org/data.zarr",
        "https://example.org/data.zarr?token=secret",
        "s3://bucket/data.zarr#secret",
    ],
)
def test_credential_urls_rejected_without_disclosing_them(tmp_path, sample, url):
    sample["Cells"]["sc_data"] = url
    with pytest.raises(ValueError, match="URLs containing") as error:
        prepare(tmp_path, sample)
    assert url not in str(error.value)


def test_arbitrary_credentials_are_not_silently_exported(tmp_path, sample):
    sample["settings"]["aws_secret_access_key"] = "sensitive"
    with pytest.raises(ValueError, match="credential field"):
        prepare(tmp_path, sample)


@pytest.mark.parametrize("mutation", ["add", "modify", "symlink", "target"])
def test_export_changes_block_publishing(tmp_path, sample, mutation, monkeypatch):
    bundle = prepare(tmp_path, sample)
    if mutation == "add":
        (bundle.directory / ".env").write_text("secret")
    elif mutation == "modify":
        (bundle.directory / "data/000/cells.h5ad").write_bytes(b"different")
    elif mutation == "symlink":
        (bundle.directory / "extra").symlink_to(tmp_path / "cells.h5ad")
    else:
        (bundle.directory / "plotly-cloud.toml").write_text('app_id = "different"\n')
    monkeypatch.setattr(
        cloud, "run_plotly", lambda *a, **kw: pytest.fail("Unreviewed bundle uploaded")
    )
    with pytest.raises(ValueError):
        cloud.publish_bundle(bundle, "App")


def test_size_limit_checked_before_copy(tmp_path, sample, monkeypatch):
    monkeypatch.setattr(cloud, "MAX_BUNDLE_BYTES", 4)
    with pytest.raises(ValueError, match="200 MiB"):
        prepare(tmp_path, sample)
    assert not list(tmp_path.glob("guanaco-cloud-*"))


def test_source_symlink_rejected(tmp_path, sample):
    link = tmp_path / "link.h5ad"
    link.symlink_to(tmp_path / "cells.h5ad")
    sample["Cells"]["sc_data"] = str(link)
    with pytest.raises(ValueError, match="Symlinks"):
        prepare(tmp_path, sample)


def test_non_data_directory_rejected(tmp_path, sample):
    fake = tmp_path / "project.zarr"
    fake.mkdir()
    sample["Cells"]["sc_data"] = str(fake)
    with pytest.raises(ValueError, match="validated Zarr"):
        prepare(tmp_path, sample)


def make_wheel(tmp_path, name="guanaco-viz"):
    wheel = tmp_path / "guanaco_viz-1.2.3-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            "guanaco_viz-1.2.3.dist-info/METADATA", f"Name: {name}\nVersion: 1.2.3\n"
        )
    return wheel


def test_unpublished_wheel_can_be_exported(tmp_path, sample):
    wheel = make_wheel(tmp_path)
    bundle = prepare(tmp_path, sample, wheel=wheel)
    assert (bundle.directory / wheel.name).read_bytes() == wheel.read_bytes()
    assert (bundle.directory / "requirements.txt").read_text() == f"./{wheel.name}\n"


def test_wrong_package_wheel_rejected(tmp_path, sample):
    with pytest.raises(RuntimeError, match="not guanaco-viz"):
        prepare(tmp_path, sample, wheel=make_wheel(tmp_path, name="unrelated"))


def test_checkout_builds_current_code_without_pypi(tmp_path, sample, monkeypatch):
    monkeypatch.setattr(cloud, "source_checkout", lambda: tmp_path)
    commands = []

    def build(command, **kwargs):
        commands.append(command)
        make_wheel(Path(command[command.index("--outdir") + 1]))

    monkeypatch.setattr(cloud, "run_process", build)
    bundle = prepare(tmp_path, sample)
    assert commands[0][:5] == [
        sys.executable,
        "-m",
        "build",
        "--wheel",
        "--no-isolation",
    ]
    assert (
        (bundle.directory / "requirements.txt").read_text().startswith("./guanaco_viz-")
    )


def test_new_publish_uses_dotted_entrypoint_and_explicit_team(
    tmp_path, sample, monkeypatch
):
    bundle = prepare(tmp_path, sample)
    commands = []
    monkeypatch.setattr(cloud, "run_plotly", lambda args, **kw: commands.append(args))
    cloud.publish_bundle(bundle, " My app ", "My team")
    assert commands[0][-6:] == [
        "--name",
        "My app",
        "--entrypoint-module",
        "app",
        "--team",
        "My team",
    ]
    assert "--skip-size-check" not in commands[0]


def test_update_preserves_identity_without_overwriting_previous_export(
    tmp_path, sample, monkeypatch
):
    previous = tmp_path / "plotly-cloud.toml"
    original = 'name = "Original 🦙"\napp_id = "app-123"\napp_url = "my-app"\n'
    previous.write_text(original)
    bundle = prepare(tmp_path, sample, existing=previous)
    assert bundle.app_id == "app-123"
    assert "UPDATE" in bundle.review()
    commands = []
    monkeypatch.setattr(cloud, "run_plotly", lambda args, **kw: commands.append(args))
    cloud.publish_bundle(bundle, "Ignored new name", "Ignored new team")
    assert "--name" not in commands[0]
    assert "--team" not in commands[0]
    assert previous.read_text() == original
    assert cloud.deployment_links(bundle.directory) == (
        "https://my-app.plotly.app",
        "https://cloud.plotly.com/app/app-123/settings",
    )


def test_update_requires_app_id(tmp_path, sample):
    previous = tmp_path / "plotly-cloud.toml"
    previous.write_text('name = "No identity"\n')
    with pytest.raises(ValueError, match="no existing app_id"):
        prepare(tmp_path, sample, existing=previous)


def test_deployment_url_cannot_escape_plotly_domain(tmp_path):
    previous = tmp_path / "plotly-cloud.toml"
    previous.write_text('app_url = "evil.example/path?"\n')
    with pytest.raises(ValueError, match="Invalid app_url"):
        cloud.deployment_links(tmp_path)


@pytest.mark.parametrize("zarr", [False, True])
def test_entrypoint_resolves_paths_before_import_from_other_cwd(
    tmp_path, sample, monkeypatch, zarr
):
    if zarr:
        store = tmp_path / "cells.zarr"
        (store / "var").mkdir(parents=True)
        (store / ".zgroup").write_text('{"zarr_format": 2}')
        (store / "var/.zattrs").write_text("{}")
        sample["Cells"]["sc_data"] = str(store)
    bundle = prepare(tmp_path, sample)
    portable = (bundle.directory / "guanaco.json").read_bytes()
    fake = types.ModuleType("guanaco.main")
    fake.app = object()
    monkeypatch.setitem(sys.modules, "guanaco.main", fake)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("GUANACO_CONFIG", "old")
    namespace = runpy.run_path(str(bundle.directory / "app.py"))
    import os

    resolved = json.loads(Path(os.environ["GUANACO_CONFIG"]).read_text())
    path = Path(resolved["Cells"]["sc_data"])
    assert path.is_absolute() and path.exists()
    if zarr:
        assert (path / "var/.zattrs").exists()
    assert namespace["app"] is fake.app
    assert (bundle.directory / "guanaco.json").read_bytes() == portable
    namespace["_runtime"].cleanup()


@pytest.mark.parametrize("format", ["h5ad", "h5mu", "zarr"])
def test_runtime_uses_real_loader(tmp_path, sample, monkeypatch, format):
    import anndata as ad
    import muon as mu
    import numpy as np
    from guanaco.data.loader import initialize_data

    data = ad.AnnData(np.zeros((4, 3)))
    path = tmp_path / f"cells.{format}"
    if format == "h5mu":
        mu.MuData({"rna": data}).write_h5mu(path)
    elif format == "zarr":
        data.write_zarr(path)
    else:
        data.write_h5ad(path)
    sample["Cells"]["sc_data"] = str(path)
    bundle = prepare(tmp_path, sample)
    fake = types.ModuleType("guanaco.main")
    fake.app = object()
    monkeypatch.setitem(sys.modules, "guanaco.main", fake)
    monkeypatch.setenv("GUANACO_CONFIG", "old")
    namespace = runpy.run_path(str(bundle.directory / "app.py"))
    registry = initialize_data(
        json_path=namespace["_resolved"], lazy_load=True, backed_mode=True
    )
    assert registry["Cells"].adata.shape == (4, 3)
    file = getattr(registry["Cells"].adata, "file", None)
    if file is not None:
        file.close()
    namespace["_runtime"].cleanup()


def test_process_drains_output_and_reports_failures():
    lines = []
    assert (
        cloud.run_process([sys.executable, "-c", "print('hello')"], emit=lines.append)
        == "hello"
    )
    assert lines == ["hello"]
    with pytest.raises(RuntimeError, match="exit 2"):
        cloud.run_process([sys.executable, "-c", "raise SystemExit(2)"])


def test_silent_process_is_cancellable():
    cancel = threading.Event()
    timer = threading.Timer(0.2, cancel.set)
    timer.start()
    try:
        with pytest.raises(cloud.Cancelled):
            cloud.run_process(
                [sys.executable, "-c", "import time; time.sleep(60)"], cancel=cancel
            )
    finally:
        timer.cancel()


def test_silent_process_times_out():
    with pytest.raises(TimeoutError):
        cloud.run_process(
            [sys.executable, "-c", "import time; time.sleep(60)"], timeout=0.1
        )


def test_optional_cli_has_actionable_install_error(monkeypatch):
    def missing(name):
        raise cloud.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(cloud.metadata, "version", missing)
    with pytest.raises(RuntimeError, match="not installed"):
        cloud.run_plotly(["--help"])


def test_cli_uses_current_interpreter_without_shell(monkeypatch):
    monkeypatch.setattr(cloud.metadata, "version", lambda name: "0.4.3")
    monkeypatch.setattr(cloud, "run_process", lambda args, **kw: args)
    command = cloud.run_plotly(["user", "teams"])
    assert command == [
        sys.executable,
        str(Path(cloud.__file__).resolve()),
        "user",
        "teams",
    ]


@pytest.mark.parametrize("command", ["publish", "status", "logs"])
@pytest.mark.parametrize("known_error", [True, False])
def test_cli_preserves_original_error_after_successful_parse(
    monkeypatch, capsys, command, known_error
):
    cli = pytest.importorskip("plotly_cloud.cli")
    from plotly_cloud.exceptions import AppCreationError

    handler = cli.CommandRegistry.commands["app"]["commands"][command]
    original_parse = cli.parse_args
    calls = []

    async def fail(args):
        calls.append(args.project_path)
        error = AppCreationError if known_error else RuntimeError
        raise error("SIMULATED original failure")

    # Replace execution before any authentication, file access or network call.
    monkeypatch.setattr(handler, "execute", fail)
    monkeypatch.setattr(
        sys, "argv", ["plotly", "app", command, "--project-path", "/local-only-test"]
    )
    with pytest.raises(SystemExit) as exit_info:
        cloud._plotly_cli_main()
    output = capsys.readouterr()
    assert exit_info.value.code == 1
    assert "SIMULATED original failure" in output.out
    assert "unknown argument" not in output.err
    assert calls == ["/local-only-test"]
    assert cli.parse_args is original_parse


def test_cli_adapter_still_rejects_genuinely_unknown_arguments(monkeypatch, capsys):
    cli = pytest.importorskip("plotly_cloud.cli")
    original_parse = cli.parse_args
    monkeypatch.setattr(
        sys, "argv", ["plotly", "app", "publish", "--not-a-real-option"]
    )
    with pytest.raises(SystemExit) as exit_info:
        cloud._plotly_cli_main()
    assert exit_info.value.code == 2
    assert "unknown argument '--not-a-real-option'" in capsys.readouterr().err
    assert cli.parse_args is original_parse
