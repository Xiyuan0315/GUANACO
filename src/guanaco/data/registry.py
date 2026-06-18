"""Runtime dataset registry built from the active GUANACO config."""

from guanaco.data.loader import JSON_PATH, initialize_data, load_config


def _bool_setting(settings: dict, key: str, default: bool) -> bool:
    value = settings.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "on"}:
            return True
        if normalized in {"false", "0", "no", "off"}:
            return False
    return default


def _backed_mode_from_settings(settings: dict) -> bool:
    return _bool_setting(settings, 'backed_mode', False)


def _embedding_render_backend_from_settings(settings: dict) -> str:
    backend = str(settings.get('embedding_render_backend', 'scattergl')).lower()
    if backend in {'scattergl', 'datashader'}:
        return backend
    print(f"Warning: invalid settings.embedding_render_backend '{backend}', using 'scattergl'.")
    return 'scattergl'


cfg = load_config(JSON_PATH)
settings = cfg.get("settings", {})
if not isinstance(settings, dict):
    settings = {}

lazy_load = _bool_setting(settings, 'lazy_load', True)
backed_mode = _backed_mode_from_settings(settings)
embedding_render_backend = _embedding_render_backend_from_settings(settings)
max_cells = settings.get('max_cells', 10000)
max_cells = None if max_cells is None else int(max_cells)

datasets = initialize_data(
    json_path=JSON_PATH,
    lazy_load=lazy_load,
    backed_mode=backed_mode,
    max_cells=max_cells,
    cfg=cfg,
)
color_config = cfg.get("color", [])
