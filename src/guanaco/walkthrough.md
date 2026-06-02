# GUANACO 单细胞可视化平台 代码导读

*2026-06-02T13:21:00Z by Showboat 0.6.1*
<!-- showboat-id: 6911219b-1f89-407b-a932-5a030b23f8af -->

## 项目概览

**GUANACO**（Genome bUild ANnotAtion COmprehensive visualization）是一个基于 [Plotly Dash](https://dash.plotly.com/) 的交互式单细胞多组学可视化 Web 平台。

### 核心功能

| 模块 | 说明 |
|------|------|
| 散点图（Scatter）| 展示 UMAP/PCA/t-SNE 等降维嵌入，支持细胞注释着色与基因表达着色 |
| 矩阵图表 | Heatmap、Dotplot、Violin Plot、Stacked Bar、Pseudotime 等 |
| 基因组浏览器（IGV）| 集成 ATAC-seq BigWig/BED 轨道，支持 JASPAR 转录因子 motif 检索 |
| 多数据集切换 | 通过顶栏 Tab 切换不同实验数据集 |
| 多模态支持 | 同时可视化 RNA + ATAC 等 MuData 多模态数据 |

### 整体执行流

    guanaco (CLI)
    │
    ├── cli.py          ← 解析命令行参数，设置环境变量
    ├── app.py          ← 初始化 Dash app + 读取 guanaco.json
    ├── data/
    │   ├── loader.py   ← 加载 .h5ad/.h5mu 文件，构建 DatasetBundle
    │   └── registry.py ← 模块级单例：datasets 注册表
    ├── layouts/        ← 构建所有 HTML/Dash 组件布局
    │   └── app_layout.py
    ├── main.py         ← 组装布局 + 注册所有回调
    └── pages/
        ├── matrix/     ← 矩阵类图表（散点、热图、提琴图等）
        │   ├── plots/  ← 纯绘图函数（不含 Dash 依赖）
        │   ├── layouts/← Dash 控件布局
        │   └── callbacks/← Dash 回调函数
        └── track/      ← 基因组浏览器（IGV + JASPAR）

## 第一步：CLI 入口 `cli.py`

用户通过命令 `guanaco -c /path/to/guanaco.json` 启动服务器。CLI 做三件事：

1. **解析参数**：配置文件 `--config`
2. **写入环境变量**：将路径和参数存入 `GUANACO_DATA_DIR`、`GUANACO_CONFIG` 等 env var，供后续模块读取
3. **延迟导入 + 进度提示**：用 Spinner 包裹耗时的库导入，避免用户面对静默等待

**关键设计**：参数通过环境变量传递给下游模块，而不是直接函数调用。这样 `guanaco.main` 等模块可以在模块级初始化时读取配置，无需修改函数签名。

```bash
sed -n '76,136p' cli.py
```

```output
def main():
    """Main entry point for GUANACO CLI."""
    import time
    import sys
    from guanaco.utils.progress_utils import Spinner
    
    args = parse_args()
    
    # Set environment variables from CLI args
    os.environ['GUANACO_DATA_DIR'] = str(args.config_dir)
    os.environ['GUANACO_CONFIG'] = str(args.config)
    os.environ['GUANACO_LAZY_LOAD'] = "true" if args.settings["lazy_load"] else "false"
    os.environ['GUANACO_BACKED_MODE'] = str(args.settings["backed_mode"]).lower()
    
    print("🧬 Starting GUANACO server...")
    print(f"📄 Config file: {args.config}")
    print(f"📁 Relative data path base: {args.config_dir}")
    print(f"🌐 Server will run on {args.settings['host']}:{args.settings['port']}")
    print(f"💾 Backed mode: {args.settings['backed_mode']}")
    print("─" * 60)
    
    start_time = time.time()
    
    try:
        # Load core libraries
        with Spinner("Loading core libraries (scanpy, muon)...Please wait 5-10 seconds..."):
            import anndata 
            import muon 
        
        # Load visualization
        with Spinner("Loading visualization libraries (matplotlib, plotly)..."):
            import matplotlib
            matplotlib.use('Agg')
            import plotly
        
        
        with Spinner("Importing GUANACO modules and Loading data..."):
            from guanaco.main import app
        
        elapsed = time.time() - start_time
        print(f"\n✅ Startup completed in {elapsed:.1f} seconds")
        print("─" * 60)
        
    except KeyboardInterrupt:
        print("\n\n❌ Startup cancelled by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error during startup: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Now actually start the server
    print()  # Empty line before Flask messages
    app.run_server(
        host=args.settings["host"],
        debug=False,
        port=args.settings["port"],
    )
```

## 第二步：Dash 应用初始化 `app.py`

`app.py` 是 Dash 实例的唯一来源（**单例模式**）。它在模块级完成两件事：

- **读取配置文件**：从 `GUANACO_CONFIG` 环境变量获取 JSON 路径，调用 `load_config()` 读取 `guanaco.json`
- **创建 Dash App**：挂载 Bootstrap LUX 主题和自定义 CSS，设置 `suppress_callback_exceptions=True`（因为回调目标是动态生成的）

其他模块通过 `from guanaco.app import app` 获取同一个 `app` 对象，注册回调时都往这个对象上挂载。

```bash
sed -n '1,24p' app.py | tr -d '\r'
```

```output
import dash
import dash_bootstrap_components as dbc
from pathlib import Path
import os
import json

def load_config(json_path):
    if not json_path.exists():
        raise FileNotFoundError(f"Config file not found: {json_path}")
    return json.loads(json_path.read_text())

# Get config path from environment variable or use default
JSON_PATH = Path(os.environ.get("GUANACO_CONFIG", "guanaco.json"))

# Initialize the app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.LUX, "/assets/scientific_style.css"],
    suppress_callback_exceptions=True
)

config = load_config(JSON_PATH)
app.title = config.get("title", "GUANACO")
```

## 第三步：数据层 `data/loader.py` + `data/registry.py`

### DatasetBundle —— 数据集的统一封装

`DatasetBundle` 是贯穿整个系统的核心数据对象，封装了一个实验数据集的所有要素：

| 属性 | 类型 | 说明 |
|------|------|------|
| `adata` | AnnData / MuData | 单细胞矩阵数据（延迟加载） |
| `genome_tracks` | dict | IGV 轨道配置（来自 S3 Bucket） |
| `ref_track` | dict | 参考基因组 URL |
| `gene_markers` | list | 默认展示的标记基因 |
| `color_config` | list | 全局离散颜色配置 |

**懒加载（Lazy Load）设计**：`adata` 属性是一个 `@property`，首次访问时才从磁盘加载文件。这样服务器启动时不会因为加载大型 `.h5ad` 文件而阻塞，只有用户切换到该数据集标签页时才触发实际 IO。

### initialize_data() —— 配置驱动的数据集构建

`initialize_data()` 遍历 `guanaco.json` 中的每个 dataset key，根据 `sc_data` 字段加载 AnnData、根据 `bucket_urls` 字段从 S3 加载 BigWig/BED 轨道，返回 `dict[str, DatasetBundle]`。

### 模块级注册表 `data/registry.py`

`registry.py` 在**模块导入时**立即调用 `initialize_data()`，把结果存入模块变量 `datasets`。所有需要数据集信息的模块都通过 `from guanaco.data.registry import datasets` 获取此单例字典。

```bash
sed -n '51,97p' data/loader.py
```

```output
class DatasetBundle:
    def __init__(
        self,
        title: str,
        description: str,
        adata: ad.AnnData | None,
        gene_markers: list[str] | None,
        label_list: list[str] | None,
        genome_tracks: dict[str, list[dict[str, Any]]] | None,
        ref_track: dict[str, str] | None,
        color_config: list[str],
        adata_path: str | None = None,
        lazy_load: bool = False,
        backed_mode: bool | str = True,
        optional_plot_components: list[str] | None = None,
    ):

        self.title = title
        self.description = description
        self._adata = adata
        self.gene_markers = gene_markers  # Only for RNA, other modalities get first 10 dynamically
        self.label_list = label_list
        self.genome_tracks = genome_tracks
        self.ref_track = ref_track
        self.color_config = color_config
        self.adata_path = adata_path
        self.lazy_load = lazy_load
        self.backed_mode = backed_mode
        self.optional_plot_components = optional_plot_components

    @property
    def adata(self):
        """Lazy load AnnData when first accessed."""
        if self.lazy_load and self._adata is None and self.adata_path:
            print(f"Loading dataset {self.title}...")
            max_cells = int(os.environ.get('GUANACO_MAX_CELLS', '8000'))
            seed = int(os.environ.get('GUANACO_SEED', '42')) if 'GUANACO_SEED' in os.environ else None
            backed = self.backed_mode if self.backed_mode else False
            self._adata = load_adata(self.adata_path, max_cells=max_cells, seed=seed, backed=backed)
            # Update gene markers and labels after loading
            # if self.gene_markers is None:
            #     self.gene_markers = self._adata.var_names[:6].tolist()
            if self.label_list is None:
                self.label_list = get_discrete_labels(self._adata)
        return self._adata
```

```bash
sed -n '118,190p' data/loader.py
```

```output
def load_adata(
    file: str | Path,
    *,
    max_cells: int | None = 10_000,
    seed: int | None = None,
    base_dir: Path = BASE_DIR,
    backed: bool | str = False,
) -> ad.AnnData | mu.MuData:
    """
    Load a single .h5ad or .h5mu file and optionally down-sample cells.
    
    Args:
        file: Path to the data file
        max_cells: Maximum number of cells to load (downsampling if needed)
        seed: Random seed for downsampling
        base_dir: Base directory for relative paths
        backed: If True, use backed mode (disk-based). If 'r+', use read-write backed mode.

    Returns:
        AnnData (for .h5ad) or MuData (for .h5mu)
    """

    path = Path(file)
    if not path.is_absolute():
        path = base_dir / path
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    if path.suffix == ".h5mu":
        if backed:
            mode = 'r+' if backed == 'r+' else True
            adata = mu.read_h5mu(path, backed=mode)
            print(f"Loaded {path} in backed mode (disk-based)")
        else:
            adata = mu.read_h5mu(path)
    elif path.suffix == ".h5ad":
        if backed:
            # Use backed mode for disk-based access
            mode = 'r+' if backed == 'r+' else 'r'
            adata = ad.read_h5ad(path, backed=mode)
            print(f"Loaded {path} in backed mode (disk-based)")
        else:
            adata = ad.read_h5ad(path)
    else:
        raise ValueError(f"Unsupported file extension: {path.suffix}")

    if max_cells is not None and not backed:
        # Handle both AnnData and MuData separately
        # Note: Downsampling is not supported in backed mode
        if isinstance(adata, ad.AnnData):
            if adata.n_obs > max_cells:
                rng = np.random.default_rng(seed)
                idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
                idx.sort()
                # Copy to release references to the full matrix.
                adata = adata[idx, :].copy()
                print(f"Down-sampled {path} to {max_cells} cells")
        elif isinstance(adata, mu.MuData):
            # For MuData, down-sample the primary modality
            primary_mod = list(adata.mod.keys())[0]
            n_obs = adata.mod[primary_mod].n_obs
            if n_obs > max_cells:
                rng = np.random.default_rng(seed)
                idx = rng.choice(n_obs, size=max_cells, replace=False)
                idx.sort()
                for mod in adata.mod:
                    # Copy per modality to avoid retaining the full object via views.
                    adata.mod[mod] = adata.mod[mod][idx, :].copy()
                print(f"Down-sampled MuData {path} to {max_cells} cells")
    elif max_cells is not None and backed:
        print(f"Warning: Downsampling not supported in backed mode for {path}")

    return adata
```

## 第四步：全局布局 `layouts/app_layout.py`

布局层负责把数据对象转化为 Dash HTML 组件树。主要函数：

- **`navbar(datasets)`**：顶栏导航，把 `datasets` 字典的 key 渲染为 `dbc.Tabs` 标签页
- **`tab_content(dataset, tab)`**：每个标签页的内容，包含三个懒渲染占位 Div：
  - `description-layout-div`：数据集描述信息
  - `ann-layout-div`：AnnData 可视化区域（散点图 + 其他图表）
  - `igv-layout-div`：基因组浏览器区域
- **`anndata_layout(adata, ...)`**：上半部分是 embedding 散点图区域，下半部分是左侧控件栏 + 右侧多标签图表区域
- **`igv_layout(session_names, prefix)`**：基因组浏览器的容器布局

**模式化懒渲染**：三个 Div 的内容不在 `tab_content` 里生成，而是在 `main.py` 注册的 `@app.callback` 中按需填充。这使得切换标签页时只触发当前数据集的渲染。

```bash
sed -n '267,285p' layouts/app_layout.py
```

```output
# Entire tab content
def tab_content(dataset, tab):
    return html.Div([
        html.Div(id={"type": "description-layout-div", "index": tab}),
        html.Div([
            create_modality_tabs(dataset, tab),
            html.Div(id={"type": "ann-layout-div", "index": tab})
        ]),
        # Wrap IGV in loading component with delay
        dcc.Loading(
            id=f"loading-igv-{tab}",
            type="circle",
            children=[
                html.Div(id={"type": "igv-layout-div", "index": tab}),
            ],
            style={"marginTop": "20px"},
            delay_show=1000,  # Wait 1 second before showing loading spinner
        )
    ])
```

## 第五步：主文件组装 `main.py`

`main.py` 是整个系统的组装点，做三件事：

### 1. 设置应用布局

将所有静态组件（导航栏、页脚、懒渲染占位 Div）挂载到 `app.layout`。

### 2. 逐数据集注册回调

对每个 DatasetBundle，根据数据类型注册对应的回调模块：

- 如果有 `adata`（AnnData 或 MuData），调用 `matrix_callbacks(app, adata, prefix, ...)`
  - MuData 会拆分为各子模态，分别以 `{dataset}-{mod}` 为 prefix 注册
- 如果有 `genome_tracks`，调用 `gene_browser_callbacks(app, ...)`

**prefix 机制**是 Dash 多实例的关键：所有 DOM ID 都以 `{prefix}-` 开头，确保不同数据集的组件 ID 不冲突。

### 3. 注册跨数据集回调

- **`update_tab_content`**：监听顶栏 Tab 切换，填充 `tabs-content`
- **`update_anndata_layout`**：监听模态 Tab 切换，重建 `ann-layout-div`（含所有图表控件）
- **`update_description_layout`**：填充数据集描述信息
- **`update_igv_layout`**：填充基因组浏览器容器

```bash
sed -n '30,76p' main.py
```

```output
app.layout = html.Div([
    dcc.Location(id="url", refresh=False),
    dcc.Store(id="tip-store", storage_type="session", data={"shown": False}),
    navbar(datasets),
    resize_tip_toast(),
    html.Div(id="tabs-content", style={"paddingTop": "70px"}),
    footprint,
    guanaco_footer,
])

# Register callbacks for scatter and other plots for each dataset
for name, dataset in datasets.items():
    dataset_adata = dataset.adata

    # Register AnnData callbacks if adata exists
    if dataset_adata is not None:
        if isinstance(dataset_adata, mu.MuData):
            for mod in dataset_adata.mod.keys():
                mod_adata = dataset_adata.mod[mod]
                prefix = f"{name}-{mod}"
                matrix_callbacks(
                    app,
                    mod_adata,
                    prefix,
                    embedding_render_backend=embedding_render_backend,
                )
        else:
            prefix = name
            matrix_callbacks(
                app,
                dataset_adata,
                prefix,
                embedding_render_backend=embedding_render_backend,
            )

    if dataset.genome_tracks is not None and dataset.ref_track is not None:
        gene_browser_callbacks(app, dataset.genome_tracks, dataset.ref_track, dataset.title)


@app.callback(
    Output("tabs-content", "children"),
    Input("tabs-dataset", "active_tab")
)
def update_tab_content(tab):
    dataset = datasets[tab]
    return tab_content(dataset, tab)

```

## 第六步：矩阵图表回调注册中心 `pages/matrix/callbacks/register.py`

这是整个可视化系统的**神经中枢**，`matrix_callbacks()` 函数为单个数据集注册所有交互回调。

### FigureMemoCache —— 图形缓存

`FigureMemoCache` 是一个带 TTL（300秒）的 LRU 内存缓存，避免相同参数下重复渲染同一张图：

- `max_items=24`：最多缓存 24 张图的序列化字典
- `_make_cache_key()`：将所有绘图参数序列化为 MD5 哈希作为 key
- 缓存命中时，图形以 `go.Figure(dict)` 方式恢复，再叠加当前的 zoom 状态（`relayoutData`）

### 全局过滤器（Global Filter）

用户可在散点图上方展开过滤面板，按元数据（obs 列）过滤细胞子集。过滤结果存入 `dcc.Store` 组件 `{prefix}-global-filtered-data`，下游的所有图表回调均以此为输入，实现联动过滤。

### generate_other_plots() —— 可选图表的动态组装

根据配置中的 `optional_plot_components` 列表，决定渲染哪些图表标签页。**配置驱动 + 运行时检查**：
- `paga`：检查 `adata.uns['paga']` 是否存在
- `volcano`：调用 `has_volcano_data(adata)` 检测是否有差异表达结果

```bash
sed -n '53,116p' pages/matrix/callbacks/register.py
```

```output
class FigureMemoCache:
    """Small in-process TTL+LRU cache for expensive callback figure payloads."""

    def __init__(self, max_items=24, ttl_seconds=300):
        self.max_items = max_items
        self.ttl_seconds = ttl_seconds
        self._store = OrderedDict()

    def _prune_expired(self):
        now = time.time()
        expired = [k for k, (_, ts) in self._store.items() if now - ts > self.ttl_seconds]
        for k in expired:
            self._store.pop(k, None)

    def get(self, key):
        self._prune_expired()
        item = self._store.get(key)
        if item is None:
            return None
        value, ts = item
        # Refresh LRU position
        self._store.move_to_end(key)
        return value

    def set(self, key, value):
        self._prune_expired()
        self._store[key] = (value, time.time())
        self._store.move_to_end(key)
        while len(self._store) > self.max_items:
            self._store.popitem(last=False)


_figure_cache = FigureMemoCache(max_items=24, ttl_seconds=300)


def _hash_list_signature(values):
    """Compact signature for potentially large lists."""
    if values is None:
        return None
    if not isinstance(values, (list, tuple)):
        return values
    n = len(values)
    if n == 0:
        return {"len": 0, "hash": None}
    # Hash full content for correctness while keeping key compact.
    payload = json.dumps(list(values), sort_keys=False, default=str, separators=(",", ":"))
    digest = hashlib.md5(payload.encode("utf-8")).hexdigest()
    return {"len": n, "hash": digest}


def _filtered_data_signature(filtered_data):
    if not filtered_data:
        return None
    return {
        "n_cells": filtered_data.get("n_cells"),
        "cell_indices": _hash_list_signature(filtered_data.get("cell_indices")),
    }


def _make_cache_key(kind, adata, **kwargs):
    payload = {"kind": kind, "adata_id": id(adata), **kwargs}
    payload_json = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.md5(payload_json.encode("utf-8")).hexdigest()

```

## 第七步：散点图回调 `pages/matrix/callbacks/scatter_callbacks.py`

散点图是 GUANACO 的核心交互入口，包含两个并排的嵌入图：

- **左图（Annotation Scatter）**：按细胞注释（obs 列）着色，支持框选/套索选取细胞子集
- **右图（Gene Scatter）**：按基因表达着色，支持共表达（co-expression）模式

### 两图联动机制

两图通过 `relayoutData` 传递 zoom 状态。当左图发生 zoom/reset 时，右图在下次触发时通过 `cross_relayout` 应用相同的视口范围（仅在两图使用相同 embedding 时才同步）。

### 细胞选择工作流

1. 用户在左图框选 → `selectedData` 触发 `update_left_highlighted_cells`，将选中细胞 ID 存入 `left-highlighted-cells-store`
2. 右图实时高亮这些细胞（通过 `highlighted_cell_ids` 参数）
3. 用户点击Update

```bash
sed -n '210,342p' pages/matrix/callbacks/scatter_callbacks.py
```

```output
    @app.callback(
        Output(f"{prefix}-annotation-scatter", "figure"),
        [
            Input(f"{prefix}-clustering-dropdown", "value"),
            Input(f"{prefix}-x-axis", "value"),
            Input(f"{prefix}-y-axis", "value"),
            Input(f"{prefix}-annotation-dropdown", "value"),
            Input(f"{prefix}-marker-size-slider", "value"),
            Input(f"{prefix}-opacity-slider", "value"),
            Input(f"{prefix}-scatter-legend-toggle", "value"),
            Input(f"{prefix}-axis-toggle", "value"),
            Input(f"{prefix}-discrete-color-map-dropdown", "value"),
            Input(f"{prefix}-scatter-log-or-zscore", "value"),
            Input(f"{prefix}-plot-order", "value"),
            Input(f"{prefix}-scatter-color-map-dropdown", "value"),
            Input(f"{prefix}-global-filtered-data", "data"),
            Input(f"{prefix}-spatial-imgkey-dropdown", "value"),
            Input(f"{prefix}-gene-scatter", "relayoutData"),
        ],
        [
            State(f"{prefix}-annotation-scatter", "relayoutData"),
            State(f"{prefix}-right-clustering-dropdown", "value"),
            State(f"{prefix}-right-x-axis", "value"),
            State(f"{prefix}-right-y-axis", "value"),
        ],
    )
    def update_annotation_scatter(
        clustering_method,
        x_axis,
        y_axis,
        annotation,
        marker_size,
        opacity,
        legend_show,
        axis_show,
        discrete_color_map,
        transformation,
        order,
        continuous_color_map,
        filtered_data,
        spatial_img_key,
        gene_relayout,
        annotation_relayout,
        right_clustering,
        right_x_axis,
        right_y_axis,
    ):
        if not annotation:
            raise exceptions.PreventUpdate

        triggered_prop = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        # Keep left/right plots synced when the right plot was zoomed/reset.
        same_embedding_view = (
            clustering_method == right_clustering
            and x_axis == right_x_axis
            and y_axis == right_y_axis
        )
        cross_relayout = gene_relayout if (
            triggered_prop == f"{prefix}-gene-scatter.relayoutData" and same_embedding_view
        ) else None
        # Preserve current view on style-only updates; reset on geometry/data changes.
        if triggered_prop in {
            f"{prefix}-clustering-dropdown.value",
            f"{prefix}-x-axis.value",
            f"{prefix}-y-axis.value",
            f"{prefix}-global-filtered-data.data",
            f"{prefix}-spatial-imgkey-dropdown.value",
        }:
            self_relayout = None
        else:
            self_relayout = annotation_relayout
        effective_relayout = cross_relayout if cross_relayout is not None else self_relayout

        plot_adata = resolve_plot_adata_from_filter(filtered_data)
        filtered_cell_idx = _filtered_cell_indices(filtered_data)
        filtered_sig = filtered_data_signature(filtered_data)

        is_continuous = annotation in adata.var_names or is_continuous_annotation(adata, annotation)
        render_backend = embedding_render_backend
        n_annotation_categories = adata.obs[annotation].nunique() if annotation in adata.obs.columns else 0
        discrete_palette = resolve_discrete_palette(
            discrete_color_map,
            n_annotation_categories,
            default=color_config,
        )
        use_transformation = transformation if annotation in adata.var_names else None
        cache_key = make_cache_key(
            "annotation_scatter",
            adata,
            clustering_method=clustering_method,
            x_axis=x_axis,
            y_axis=y_axis,
            annotation=annotation,
            marker_size=marker_size,
            opacity=opacity,
            render_backend=render_backend,
            legend_show=legend_show,
            axis_show=axis_show,
            discrete_color_map=discrete_color_map,
            transformation=use_transformation if is_continuous else None,
            order=order if is_continuous else None,
            continuous_color_map=continuous_color_map or "Viridis",
            filtered_data=filtered_sig,
            spatial_img_key=spatial_img_key,
            mode="continuous" if is_continuous else "categorical",
        )
        cached_fig = cached_figure_get(cache_key)
        if cached_fig is not None:
            return apply_relayout(cached_fig, effective_relayout)

        fig = plot_embedding(
            adata=plot_adata,
            adata_full=adata,
            embedding_key=clustering_method,
            color=annotation,
            x_axis=x_axis,
            y_axis=y_axis,
            mode="continuous" if is_continuous else "categorical",
            transformation=use_transformation if is_continuous else None,
            order=order if is_continuous else None,
            continuous_color_map=continuous_color_map or "Viridis",
            discrete_color_map=discrete_palette,
            marker_size=marker_size,
            opacity=opacity,
            render_backend=render_backend,
            legend_show=legend_show,
            axis_show=axis_show,
            img_key=spatial_img_key,
            source_adata=adata,
            cell_indices=filtered_cell_idx,
        )
        cached_figure_set(cache_key, fig)
        return apply_relayout(fig, effective_relayout)
```

## 第八步：嵌入图绘制函数 `pages/matrix/plots/embedding.py`

`plot_embedding()` 是实际绘图的核心函数，根据 `mode` 参数（`categorical` 或 `continuous`）分支处理：

- **Categorical 模式**：按 obs 列的类别分组，每个类别生成一条 Plotly trace（`Scattergl` 或 Datashader 栅格），颜色由离散调色板配置
- **Continuous 模式**：单条 trace，颜色映射到基因表达值或连续型 obs 列，支持 log/z-score 变换
- **Spatial 模式**：从 `adata.uns['spatial']` 读取组织切片图像，以 `go.Image` 作为背景层叠加散点

**Datashader 渲染路径**：`embedding_render_backend='datashader'` 时，服务端将百万级散点栅格化为 PNG 图像，通过 `go.Image` 返回给前端。这避免了向浏览器传输海量点的序列化开销。

### 颜色工具 `utils/colors.py`

颜色系统支持三大来源：
- **Plotly** 内置色阶（仅感知均匀子集）
- **Colorcet**（前缀 `cc:`）：线性/发散色阶
- **Crameri 科学色彩图**（前缀 `cmc:`）：e.g. `cmc:batlow`
- **Glasbey 动态离散调色板**：按实际类别数动态生成最大化色差的颜色列表

```bash
grep -n 'def plot_embedding\|def plot_coexpression\|datashader\|go.Image\|Scattergl' pages/matrix/plots/embedding.py | head -25
```

```output
212:def _build_datashader_continuous_figure(
227:        import datashader as ds
229:        return None, f"datashader unavailable: {exc}"
268:        return None, "datashader aggregation produced empty grid"
321:    datashader_err = None
323:    if render_backend == "datashader" and not has_highlight:
324:        ds_fig, datashader_err = _build_datashader_continuous_figure(
348:            go.Scattergl(
363:        go.Scattergl(
394:    if datashader_err:
400:            text=f"Datashader fallback to scattergl ({datashader_err})",
410:def plot_embedding(
535:    fig.add_trace(go.Scattergl(
551:        fig.add_trace(go.Scattergl(
608:def plot_coexpression_embedding(
681:        fig.add_trace(go.Scattergl(
701:            fig.add_trace(go.Scattergl(
```

## 第九步：其他图表模块

其他图表均遵循相同的三层结构：**layout（控件布局）→ callback（数据处理）→ plot（纯绘图函数）**。

| 图表 | 核心绘图函数 | 主要特点 |
|------|------------|---------|
| Heatmap | `plot_unified_heatmap()` | 支持基因 × 细胞注释的聚合热图，行列均可层次聚类 |
| Dotplot | `plot_dot_matrix()` | 点大小编码表达比例，点颜色编码平均表达量；支持自定义行列排序 |
| Violin Plot | `plot_violin1()` / `plot_violin2_new()` | 两种展示模式：按基因分组（violin1）vs 按注释分组（violin2） |
| Stacked Bar | `plot_stacked_bar()` | 展示不同注释类别在另一注释维度上的组成比例 |
| Pseudotime | `plot_genes_in_pseudotime()` | 基因表达沿拟时间轴的趋势曲线 |
| PAGA | `build_paga_cytoscape()` | 细胞间轨迹推断图，用 Cytoscape 渲染 |
| Volcano | 内联于 callback | 差异表达火山图，数据来自 `adata.uns` |
| GRN Demo | 内联于 callback | 基因调控网络（Gene Regulatory Network）演示图 |

**细胞子集联动**：所有图表的回调都接受 `{prefix}-selected-cells-store` 作为 Input，通过 `filter_data(adata, annotation, selected_labels, selected_cells)` 将散点图的细胞选择传递给下游图表，实现真正的联动分析。

```bash
grep -n 'def register_\|def plot_\|def build_' pages/matrix/callbacks/register.py | head -20
```

```output
```

```bash
grep -rn 'def register_\|def plot_\|def build_' pages/matrix/ | grep -v '__pycache__' | head -25
```

```output
pages/matrix//callbacks/stacked_bar_callbacks.py:9:def register_stacked_bar_callbacks(
pages/matrix//callbacks/paga_callbacks.py:20:def register_paga_callbacks(
pages/matrix//callbacks/heatmap_callbacks.py:8:def register_heatmap_callbacks(
pages/matrix//callbacks/violin_callbacks.py:10:def register_violin_callbacks(
pages/matrix//callbacks/dotplot_callbacks.py:6:def register_dotplot_callbacks(
pages/matrix//callbacks/scatter_callbacks.py:12:def register_scatter_callbacks(
pages/matrix//callbacks/volcano_callbacks.py:33:def register_volcano_callbacks(app, adata, prefix):
pages/matrix//callbacks/grn_demo_callbacks.py:11:def register_grn_demo_callbacks(app, adata, prefix):
pages/matrix//callbacks/pseudotime_callbacks.py:8:def register_pseudotime_callbacks(
pages/matrix//plots/volcano.py:412:def plot_volcano(
pages/matrix//plots/heatmap.py:197:def plot_unified_heatmap(
pages/matrix//plots/heatmap.py:367:def plot_heatmap2_continuous(
pages/matrix//plots/embedding.py:410:def plot_embedding(
pages/matrix//plots/embedding.py:608:def plot_coexpression_embedding(
pages/matrix//plots/stacked_bar.py:4:def plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=None, y_order=None, x_order=None):
pages/matrix//plots/pseudotime.py:19:def plot_genes_in_pseudotime(
pages/matrix//plots/violin2.py:418:def plot_violin2_new(adata, key, meta1, meta2, mode, transformation='log',
pages/matrix//plots/grn_demo.py:192:def build_grn_cytoscape(
pages/matrix//plots/paga.py:347:def build_paga_cytoscape(
pages/matrix//plots/violin1.py:158:def plot_violin1(
pages/matrix//plots/dotmatrix.py:12:def plot_dot_matrix(
```

## 第十步：基因组浏览器 `pages/track/`

### IGV 嵌入 `pages/track/callbacks.py`

基因组浏览器使用 `dash-bio` 的 `dashbio.Igv` 组件，核心回调 `return_igv` 根据用户选择的 ATAC-seq session，动态加载对应的轨道列表：

```python
return dashbio.Igv(
    genome=ref_track['label'],   # e.g. 'hg38'
    locus='chr1:1-10000000',
    tracks=genome_tracks[selected_atac]  # list of BigWig/BED 轨道配置
)
```

轨道数据在启动时由 `load_tracks_from_s3()` 从 AWS S3（或兼容的对象存储）批量获取，支持 BigWig、BED、BigBed、BEDPE 四种格式。

### JASPAR Motif 检索 `pages/track/callbacks.py`

用户输入 JASPAR 2024 的 matrix ID（如 `MA0001.1`），回调通过 `pyjaspar` 库查询转录因子信息，调用 `plot_motif()` 生成序列 logo（Position Weight Matrix 可视化），以 base64 内嵌图像方式返回。

```bash
sed -n '8,53p' pages/track/callbacks.py
```

```output
def gene_browser_callbacks(app, genome_tracks, ref_track, prefix):
    """
    Register genome browser callbacks for a specific dataset.
    Arguments:
        app: Dash app
        genome_tracks: dict of genome tracks (from DatasetBundle)
        ref_track: reference genome (from DatasetBundle)
        prefix: dataset name (used to ensure unique IDs)
    """
    if genome_tracks is None or ref_track is None:
        return

    @app.callback(
        Output(f'{prefix}-igv-container', 'children'),
        Input(f'{prefix}-igv-genome-select', 'value')
    )
    def return_igv(selected_atac):
        if selected_atac is None:
            return html.Div(
                [
                    html.P(
                        "Please select an IGV session from the dropdown above to view the genome browser.",
                        style={
                            'textAlign': 'center',
                            'color': '#868e96',
                            'fontSize': '16px',
                            'padding': '40px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '8px',
                            'margin': '20px 0'
                        }
                    )
                ]
            )

        if genome_tracks.get(selected_atac) is None:
            raise Exception(f"No tracks configured for ATAC {selected_atac}")

        return html.Div([
            dashbio.Igv(
                id=f'igv-{prefix}',
                genome=ref_track['label'],
                locus='chr1:1-10000000',
                tracks=genome_tracks[selected_atac]
            )
        ])
```

## 第十一步：完整数据流

    用户执行: guanaco -c /path/to/guanaco.json -p 4399
    │
    ├── cli.py: parse_args()
    │   ├── 设置 GUANACO_DATA_DIR, GUANACO_CONFIG, GUANACO_MAX_CELLS
    │   └── 延迟 import guanaco.main → 触发下列模块级初始化
    │
    ├── app.py (模块级)
    │   └── dash.Dash() + load_config(JSON_PATH) → app 单例
    │
    ├── data/registry.py (模块级)
    │   └── initialize_data()
    │       ├── 遍历 guanaco.json 每个 dataset key
    │       ├── sc_data → DatasetBundle(lazy_load=True)   [不立即加载 .h5ad]
    │       └── bucket_urls → load_tracks_from_s3() → genome_tracks
    │
    ├── main.py (模块级)
    │   ├── app.layout = navbar + 占位 Div
    │   ├── for dataset in datasets:
    │   │   ├── matrix_callbacks(app, adata, prefix)    [注册散点/图表回调]
    │   │   └── gene_browser_callbacks(app, tracks, ref, prefix)
    │   └── @app.callback update_tab_content / update_anndata_layout / ...
    │
    └── app.run_server(host, port)   ← Flask 开始监听
        │
        用户浏览器请求
        │
        ├── 切换 Tab → update_tab_content()
        │   └── tab_content(dataset) → 返回三个空 Div
        │
        ├── 模态 Tab 切换 → update_anndata_layout()
        │   ├── dataset.adata 首次访问 → 触发 lazy load (.h5ad 从磁盘读入)
        │   └── anndata_layout() → 返回散点图控件 + 图表 Tabs
        │
        ├── 用户选择标注/基因 → update_annotation_scatter() / update_gene_scatter()
        │   ├── 检查 FigureMemoCache (MD5 key)
        │   ├── 缓存未命中 → plot_embedding() → Plotly Figure
        │   └── apply_relayout(fig, relayoutData) → 返回前端
        │
        ├── 框选细胞 → update_left_highlighted_cells() → store → 右图高亮
        ├── 点击 Update Plots → store_selected_cells() → 更新下游所有图表
        │
        └── 选择 IGV session → return_igv() → dashbio.Igv() 组件
            └── 浏览器直接从 S3 URL 加载 BigWig/BED 数据流
